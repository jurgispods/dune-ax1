#ifndef DUNE_AX1_EVALGRIDSSETUP_HH
#define DUNE_AX1_EVALGRIDSSETUP_HH

#include <dune/common/array.hh>

#include <dune/pdelab/backend/istlvectorbackend.hh>
//#include<dune/pdelab/finiteelementmap/q12dfem.hh>
//#include<dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/common/clock.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/function/selectcomponent.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/backend/seqistlsolverbackend.hh>
#include <dune/pdelab/backend/novlpistlsolverbackend.hh>
#include <dune/pdelab/backend/ovlpistlsolverbackend.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/ordering/transformations.hh>
#include <dune/pdelab/ordering/entityblockedlocalordering.hh>
#include <dune/pdelab/ordering/interleavedordering.hh>
#include <dune/pdelab/ordering/permutationordering.hh>

#include <dune/ax1/common/ax1_newton.hh>
#include <dune/ax1/common/ax1_linearproblem.hh>
#include <dune/ax1/common/ax1_discretegridfunction.hh>
#include <dune/ax1/common/ax1_parallelconstraintshelper.hh>
#include <dune/ax1/common/ax1_simulationdata.hh>
#include <dune/ax1/common/ax1_simulationstate.hh>
#include <dune/ax1/common/ax1_solverbackend.hh>
#include <dune/ax1/common/chargedensitygridfunction.hh>
#include <dune/ax1/common/error_norms.hh>
#include <dune/ax1/common/ionicstrengthgridfunction.hh>
#include <dune/ax1/common/membranefluxgridfunction.hh>
#include <dune/ax1/common/poisson_boltzmann_concentrationgridfunction.hh>
#include <dune/ax1/common/poisson_boltzmann_rhs_gridfunction.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/vectordiscretegridfunctiongradient.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_boundary.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_geometrytools.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_initial.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_output.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_setup.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_solutionvectors.hh>
#include <dune/ax1/acme2_cyl/common/nernst_planck_parameters.hh>
#include <dune/ax1/acme2_cyl/common/poisson_parameters.hh>
#include <dune/ax1/acme2_cyl/common/poisson_boltzmann_parameters.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_simulation.hh>
#include <dune/ax1/acme2_cyl/operator/convectiondiffusionfem.hh>
#include <dune/ax1/acme2_cyl/operator/poisson_boltzmann_operator.hh>
#include <dune/ax1/acme2_cyl/operator/acme2_cyl_operator_fully_coupled.hh>
#include <dune/ax1/acme2_cyl/operator/acme2_cyl_toperator.hh>

template<class Grid, class GV, class PHYSICS, class SubGV>
class EvalGridsSetup
{
  public:
    
    typedef Acme2CylSetup<Grid,GV,PHYSICS,SubGV> BaseT;
    typedef typename BaseT::Traits Traits;

    // Choose domain and range field type
    typedef GV GridView;
    typedef typename GV::Grid GridType;
    typedef typename GV::Grid::ctype Coord;
    typedef double Real;

    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

    // SubGrid stuff
    typedef typename SubGV::Grid SubGrid;

    static const bool writeIfNotConverged = true;

    EvalGridsSetup(std::vector<Dune::shared_ptr<Grid> >& grids_, GV& gv_, PHYSICS& physics_,
        SubGV& elecGV_, SubGV& membGV_)
    : grids(grids_),
      gv(gv_),
      physics(physics_),
      elecGV(elecGV_),
      membGV(membGV_)
    {}

    void setup (double dtstart, double tend)
    {
      Real time = 0.0;
      Real dt = dtstart;
      const Acme2CylParameters& params = physics.getParams();

      Real tEquilibrium = physics.getParams().tEquilibrium();
      if(params.doEquilibration() && tEquilibrium > time)
      {
        dt = physics.getParams().dtEquilibrium();
        debug_info << "== Starting with initial dt = " << dt << " until tEquilibirum " << tEquilibrium
            << " ==" << std::endl;
      }

      const int dim = GV::dimension;
      const int degreePot = 1, degreeCon = 1;

			// Finite element map for CG
      typename Traits::FEM_CON femCon;
      typename Traits::FEM_POT femPot;

      // =================== GFS SETUP ================================================================
#if USE_PARALLEL == 1 && USE_OVERLAP == 0
      typename Traits::CONSTRAINTS constraints(gv);
#else
      typename Traits::CONSTRAINTS constraints;
#endif

      // Single-component GFS for one ion species
      typename Traits::GFS_SINGLE_CON gfsSingleCon(elecGV,femCon,constraints);
#if USE_PARALLEL == 1 && USE_OVERLAP == 0
      if(gv.grid().ghostSize(0) > 0)
      {
        constraints.compute_ghosts(gfsSingleCon);
      }
#endif

      // Power grid function space for all ion species
      typename Traits::GFS_CON gfsCon(gfsSingleCon);

      // Full GFS for potential on electrolyte and membrane subdomain
      typename Traits::GFS_POT gfsPot(gv,femPot,constraints);
#if USE_PARALLEL == 1 && USE_OVERLAP == 0
      if(gv.grid().ghostSize(0) > 0)
      {
        constraints.compute_ghosts(gfsPot);
      }
#endif

      // Multidomain GFS
#ifdef MULTIPLE_MEMBRANE_ELEMENTS
      if(params.nMembraneElements() < 2 || !params.refineMembrane())
      {
        debug_warn << std::endl;
        debug_warn << "========================================================================" << std::endl;
        debug_warn << "You are using the 'multiple membrane elements' version acme2_cyl_par_mme, but the number"
          << " of membrane elements in the config file is " << params.nMembraneElements() << "! It is recommended"
          << " to use the executable acme2_cyl_par instead, which is more performant as it uses a blocked matrix!"
          << std::endl;
        debug_warn << "========================================================================" << std::endl;
        debug_warn << std::endl;
      }

      // Permutation for DOFs
      std::vector<std::size_t> permutation, inv_permutation;

      // Set initial permutation to identity
      permutation.resize(gfsCon.size()+gfsPot.size());
      inv_permutation.resize(gfsCon.size()+gfsPot.size());
      for(std::size_t i = 0; i<(gfsCon.size()+gfsPot.size()); i++)
      {
        permutation[i] = i;
        inv_permutation[i] = i;
      }

      // Reorder matrix by grouping together DOFs at vertices
      // Old method (error-prone!)
//      if(params.doReorderMatrix())
//      {
//        Tools::setupDOFPermutation(physics, gfsCon.size(), gfsPot.size(), permutation, inv_permutation);
//      }

      // Create MultiGFS with permutation ordering (identity for now)
      Dune::PDELab::PermutationOrderingTag permutationOrderingTag(inv_permutation);
      typename Traits::MultiGFS multigfs(grid,permutationOrderingTag,gfsCon,gfsPot);

      if(! params.doReorderMatrix())
      {
        DUNE_THROW(Dune::Exception, "It is necessary to use matrix reordering when using more than one membrane element layer!");
      }

      if(params.doReorderMatrix())
      {
        // Use multigfs object to set up the real permutation and update in ordering
        physics.setupDOFPermutation(multigfs, permutation, inv_permutation);
        multigfs.orderingTag().updatePermutation(inv_permutation);
      }

      // Check calculated permutation
//      for(int i=0; i<permutation.size(); i++)
//      {
//        debug_jochen << "perm[" << i << "] = " << permutation[i] << " -- " << "inv_perm[" << i << "] = " << inv_permutation[i] << std::endl;
//        assert(inv_permutation[permutation[i]] == i);
//      }

      // Check calculated permutation
//      for(int i=0; i<multigfs.orderingTag().permutation().size(); i++)
//      {
//        debug_jochen << "perm[" << i << "] = " << multigfs.orderingTag().permutation()[i] << std::endl;
//      }

      //typename Traits::MultiGFS multigfs(grid,gfsCon,gfsPot);
#else
      if(params.nMembraneElements() > 1)
      {
        DUNE_THROW(Dune::Exception, "Number of membrane elements > 1, this does not work with this executable, use"
          << " acme2_cyl_par_mme instead!");
      }
      if(params.doReorderMatrix())
      {
        DUNE_THROW(Dune::Exception, "Matrix reordering makes no sense when using only one membrane element layer!");
      }

      std::vector<std::size_t> blocksizes(2);
      blocksizes[0] = NUMBER_OF_SPECIES;
      blocksizes[1] = 1;
      Dune::PDELab::InterleavedOrderingTag interleavedOrderingTag(blocksizes);
      typename Traits::MultiGFS multigfs(*grids.back(),interleavedOrderingTag,gfsCon,gfsPot);

      // Fallback multigfs creation for LexicographicOrderingTag
      //typename Traits::MultiGFS multigfs(grid,gfsCon,gfsPot);
#endif
      debug_verb << "Instantiated MultiGFS." << std::endl;

      debug_verb << "Setting up concentration Sub GFS..." << std::endl;
      // SubGridFunctionSpaces for use in gridfunctions
      typename Traits::GFS_CON_SUB gfsConSub(multigfs);
      typename Traits::GFS_POT_SUB gfsPotSub(multigfs);
      // ==============================================================================================


      // =========== Define solution vectors ==========================================================
      // Coefficient vectors on subdomains
      typename Traits::U uold(multigfs,0.0);
      typename Traits::U unew = uold;

      debug_info << "=========== #DOFs =============" << std::endl;
      int rootNode = params.general.get("rootOutputNode", 0);
      int nBlockedDofs = unew.N();
      int nFlatDofs = unew.flatsize();

      std::vector<int> vBlockedDofs(gv.comm().size());
      std::vector<int> vFlatDofs(gv.comm().size());
      gv.comm().gather(&nBlockedDofs,&vBlockedDofs[0],1,rootNode);
      gv.comm().gather(&nFlatDofs,&vFlatDofs[0],1,rootNode);

      int nTotalBlockedDofsWithOverlap = 0;
      int nTotalFlatDofsWithOverlap = 0;
      for(int i=0; i<gv.comm().size(); i++)
      {
        debug_info << "p" << i << " |  U.N() / U.flatsize() : " << vBlockedDofs[i] << " / " << vFlatDofs[i] << std::endl;
        nTotalBlockedDofsWithOverlap += vBlockedDofs[i];
        nTotalFlatDofsWithOverlap += vFlatDofs[i];
      }
      debug_info << "Total #DOFs on all processors (including overlap): " << nTotalBlockedDofsWithOverlap << " / "
          << nTotalFlatDofsWithOverlap << std::endl;

      // Total 'real' number of DOFs (substract missing concentration DOFs on inner membrane elements)
      int nTotalDofsWithoutOverlap = (params.nNodes() * (NUMBER_OF_SPECIES+1))
          - ((params.nMembraneElements()-1) * params.X().size() * NUMBER_OF_SPECIES);
      debug_info << "Total #DOFs of the sequential problem (i.e., without overlap): " << nTotalDofsWithoutOverlap << std::endl;
      debug_info << "===============================" << std::endl;


      // Now let's put 'em all into one structure!
      Acme2CylSolutionVectors<Traits> solutionVectors(uold, unew);
      // ==============================================================================================

      // ====== Create grid functions and interpolate initial values onto coefficient vectors =========
      // Initial concentrations
      typename Traits::INITIAL_GF_CON initialGFCon(gv,physics.getParams());

      // Generic wrapper class for initial grid function
      typename Traits::INITIAL_CON initialCon(gv,physics,initialGFCon);
      initialCon.setTime(time);

      //typedef Dune::PDELab::VectorDiscreteGridFunction <typename Traits::MultiGFS,typename Traits::U> DGF_ALL;
      //DGF_ALL(multigfs,unew);

      typedef typename Dune::PDELab::GridFunctionSubSpace<typename Traits::GFS_CON_SUB,
          Dune::PDELab::TypeTree::TreePath<0> > GFS_NA;
      typedef typename Dune::PDELab::GridFunctionSubSpace<typename Traits::GFS_CON_SUB,
          Dune::PDELab::TypeTree::TreePath<1> > GFS_K;
      typedef typename Dune::PDELab::GridFunctionSubSpace<typename Traits::GFS_CON_SUB,
          Dune::PDELab::TypeTree::TreePath<2> > GFS_CL;
      GFS_NA gfsNa(gfsConSub);
      GFS_K gfsK(gfsConSub);
      GFS_CL gfsCl(gfsConSub);
      typedef Dune::PDELab::DiscreteGridFunction<GFS_NA,typename Traits::U> DGF_NA;
      typedef Dune::PDELab::DiscreteGridFunction<GFS_K,typename Traits::U> DGF_K;
      typedef Dune::PDELab::DiscreteGridFunction<GFS_CL,typename Traits::U> DGF_CL;
      DGF_NA dgfNa(gfsNa,unew);
      DGF_K dgfK(gfsK,unew);
      DGF_CL dgfCl(gfsCl,unew);

      typename Traits::DGF_CON_SUB dgfConElec(gfsConSub,unew);
      typename Traits::DGF_CON_SUB dgfOldConElec(gfsConSub,uold);
      typename Traits::DGF_CON_GRAD_SUB dgfGradConElec(gfsConSub, unew);

      typename Traits::DGF_CON dgfCon(gv, dgfConElec, physics);
      typename Traits::DGF_CON dgfOldCon(gv, dgfOldConElec, physics);
      typename Traits::DGF_CON_GRAD dgfGradCon(gv, dgfGradConElec, physics);

      typename Traits::GF_CD gfChargeDensity(dgfCon, dgfOldCon, physics);
      typename Traits::GF_IS gfIonicStrength(dgfCon, physics);

      // Initial potential
      typename Traits::INITIAL_GF_POT initialGFPot(gv,physics.getParams());
      typename Traits::INITIAL_POT initialPot(gv,physics,initialGFPot);
      initialPot.setTime(time);



      typename Traits::DGF_POT dgfPot(gfsPotSub, unew);
      typename Traits::DGF_POT dgfOldPot(gfsPotSub, uold);
      typename Traits::DGF_POT_GRAD dgfGradPot(gfsPotSub, unew);

      // Flag 'true' for updating channel states in each time step
      typename Traits::GF_MEMB_FLUX gfMembFlux(dgfOldCon, dgfOldPot, physics, true, tEquilibrium);

      // Composite initial grid function for the electrolyte subdomain
      typename Traits::INITIAL_ELEC initialElec(initialCon,initialPot);
      // ==============================================================================================


      // ========== Define Parameter classes containing the model problem =============================
      typename Traits::BOUNDARY_CON boundaryCon(physics.getParams(), initialGFCon);
      typename Traits::PARAMETERS_CON parametersCon(gv,physics,boundaryCon,gfMembFlux,tEquilibrium);

      typename Traits::BOUNDARY_POT boundaryPot(physics.getParams(), initialGFPot);
      typename Traits::PARAMETERS_POT parametersPot(physics,boundaryPot);
      // ==============================================================================================


      // ==============================================================================================
      // Define and instantiate output class
      typedef Acme2CylOutput<Traits> ACME2CYL_OUTPUT;
      ACME2CYL_OUTPUT acme2_cylOutput(gv,elecGV,membGV,multigfs,gfsConSub,gfsPotSub,solutionVectors,
          gfMembFlux,physics,parametersCon,2*degreePot,2*degreeCon);
      // ==============================================================================================

      // ========== Load saved state ======================================================================
      if(physics.getParams().doLoadState())
      {
        //solutionVectors.printHostGridVectors();
        std::map<std::string,std::string> loadFileNames = acme2_cylOutput.loadState(time,dt);
        debug_info << "======================================================================================"
            << std::endl;
        debug_info << "Loaded simulation state from following files: " << std::endl;
        for(std::map<std::string,std::string>::const_iterator it = loadFileNames.begin();
        		it != loadFileNames.end(); ++it)
        {
          debug_info << "'" << it->first << "': " << it->second << std::endl;
        }
        debug_info << "======================================================================================"
            << std::endl;
      }


      // ========== Do grid convergence test =================================================================
      typedef typename Acme2CylOutput<Traits>::Traits OutputTraits;

      std::vector<double> l2ErrorsPot(grids.size());
      std::vector<double> l2ErrorsCon(grids.size());
      std::vector<double> maxErrorsPot(grids.size());
      std::vector<double> maxErrorsCon(grids.size());
      std::vector<int> dofs(grids.size());

      // Construct gridfunction spaces for each of the grids
      std::vector<int> default_levels(grids.size(),0);
      std::vector<int> levels = params.general.get("gridLevels",default_levels);

      std::vector<std::string> default_filenames(grids.size(),std::string("bla"));
      std::vector<std::string> load_filenames = params.general.get("gridStateFiles",default_filenames);
      for(int i=0; i<grids.size(); i++)
      {
        debug_info << std::endl << std::endl << std::endl;
        debug_info << " ============================================================= " << std::endl;
        debug_info << "  GRID #" << i << " (level " << levels[i] << ")" << std::endl;
        debug_info << " ============================================================= " << std::endl;

        GV gv_coarse = grids[i]->leafGridView();

        typedef typename Grid::SubDomainGrid SubDomainGrid;
        SubDomainGrid& elecGrid_coarse = grids[i]->subDomain(0);
        SubDomainGrid& membGrid_coarse = grids[i]->subDomain(1);

        SubGV elecGV_coarse = elecGrid_coarse.leafGridView();
        SubGV membGV_coarse = membGrid_coarse.leafGridView();

        // Reconstruct pot GFS
        typename Traits::FEM_POT femPot_coarse;
        typename Traits::GFS_POT gfsPot_coarse(gv_coarse,femPot_coarse);

        // Reconstruct con GFS
        typename Traits::FEM_CON femCon_coarse;
        typename Traits::GFS_SINGLE_CON gfsSingleCon_coarse(elecGV_coarse,femCon_coarse);
        typename Traits::GFS_CON gfsCon_coarse(gfsSingleCon_coarse);

      // Reconstruct MultiGFS
#ifdef MULTIPLE_MEMBRANE_ELEMENTS
        DUNE_THROW(Dune::NotImplemented, "Copy stuff from simloader::extrapolateData() here!");
#else
        // Use InterleavedOrdering
        typename Traits::MultiGFS multigfs_coarse(*grids[i], {NUMBER_OF_SPECIES,1},
            gfsCon_coarse, gfsPot_coarse);
        // Fallback call when using LexicographicOrderingTag
        //typename Acme2CylTraits::MultiGFS multigfs_Equi(md_equilibrationGrid,
        //      gfsCon_Equi_TEMP, gfsPot_Equi_TEMP);
#endif

        // Extract sub GFS's
        typename Traits::GFS_CON_SUB gfsConSub_coarse(multigfs_coarse);
        typename Traits::GFS_POT_SUB gfsPotSub_coarse(multigfs_coarse);

        GFS_NA gfsNa_coarse(gfsConSub_coarse);
        GFS_K gfsK_coarse(gfsConSub_coarse);
        GFS_CL gfsCl_coarse(gfsConSub_coarse);

        // Load solution vectors (one per loaded simulation state)
        typename Traits::U uold_coarse(multigfs_coarse,0.0);
        typename Traits::U unew_coarse(multigfs_coarse,0.0);

        // TODO load data into solution vector
        Ax1SimulationState<Real> state;
        state.loadState(load_filenames[i]);

        debug_jochen << "Ok, loaded state: nElems = " << state.getData().getNElements() << std::endl;

        // Get simulation data from first entry in states vector
        Ax1SimulationData<Real>& simulationData_coarse = state.getData();
        Acme2CylSolutionVectors<Traits> solutionVectors_coarse(uold_coarse, unew_coarse);

        solutionVectors_coarse.deserialize(simulationData_coarse);

//        // Load channel states
//        if(physics.getParams().doLoadChannelStates())
//        {
//          if(not simulationData.hasVector("channelStates"))
//            DUNE_THROW(Dune::Exception, "No vector 'channelStates' found in saved file!");
//
//          // Now load that shit!
//          std::vector<Real> stdVec = simulationData_coarse.getVector("channelStates").data;
//        }

        // One potential/concentration gridfunction for each loaded solution vector
        std::vector<Dune::shared_ptr<typename OutputTraits::DGF_POT> > dgfPot_coarse;
        std::vector<Dune::shared_ptr<typename OutputTraits::DGF_CON> > dgfCon_coarse;
        std::vector<Dune::shared_ptr<DGF_NA> > dgfNa_coarse;
        std::vector<Dune::shared_ptr<DGF_K> > dgfK_coarse;
        std::vector<Dune::shared_ptr<DGF_CL> > dgfCl_coarse;
        //std::vector<Dune::shared_ptr<DGF_ALL> > dgfAll_coarse;

        Dune::shared_ptr<typename OutputTraits::DGF_POT> dgfPot_coarse_ptr(
           new typename OutputTraits::DGF_POT(gfsPotSub_coarse, unew_coarse));
         dgfPot_coarse.push_back(dgfPot_coarse_ptr);

       Dune::shared_ptr<typename OutputTraits::DGF_CON> dgfCon_coarse_ptr(
          new typename OutputTraits::DGF_CON(gfsConSub_coarse, unew_coarse));
        dgfCon_coarse.push_back(dgfCon_coarse_ptr);

        Dune::shared_ptr<DGF_NA> dgfNa_coarse_ptr(
          new DGF_NA(gfsNa_coarse, unew_coarse));
        dgfNa_coarse.push_back(dgfNa_coarse_ptr);

        Dune::shared_ptr<DGF_K> dgfK_coarse_ptr(
          new DGF_K(gfsK_coarse, unew_coarse));
        dgfK_coarse.push_back(dgfK_coarse_ptr);

        Dune::shared_ptr<DGF_CL> dgfCl_coarse_ptr(
          new DGF_CL(gfsCl_coarse, unew_coarse));
        dgfCl_coarse.push_back(dgfCl_coarse_ptr);

        //Dune::shared_ptr<DGF_ALL> dgfAll_coarse_ptr(
        //  new DGF_ALL(multigfs_coarse, unew_coarse));
        //dgfAll_coarse.push_back(dgfAll_coarse_ptr); // TODO Use this for calculations

        // Calculate number of nodes of this coarse grid (we don't have a params class for it);
        // assume levels.back() is the finest level
        int levelDiff = levels.back() - levels[i];

        int nNodesXCoarse = state.getData().getVector("x").data.size();
        int nNodesYCoarse = state.getData().getVector("y").data.size();

        debug_jochen << "Number of nodes for coarse grid #" << i << ": " << nNodesXCoarse << " / "
            << nNodesYCoarse << std::endl;

        debug_info << "=========== #DOFs =============" << std::endl;
        int rootNode = params.general.get("rootOutputNode", 0);
        int nBlockedDofs = unew_coarse.N();
        int nFlatDofs = unew_coarse.flatsize();

        std::vector<int> vBlockedDofs(gv.comm().size());
        std::vector<int> vFlatDofs(gv.comm().size());
        gv.comm().gather(&nBlockedDofs,&vBlockedDofs[0],1,rootNode);
        gv.comm().gather(&nFlatDofs,&vFlatDofs[0],1,rootNode);

        int nTotalBlockedDofsWithOverlap = 0;
        int nTotalFlatDofsWithOverlap = 0;
        for(int i=0; i<gv.comm().size(); i++)
        {
          debug_info << "p" << i << " |  U.N() / U.flatsize() : " << vBlockedDofs[i] << " / " << vFlatDofs[i] << std::endl;
          nTotalBlockedDofsWithOverlap += vBlockedDofs[i];
          nTotalFlatDofsWithOverlap += vFlatDofs[i];
        }
        debug_info << "Total #DOFs on all processors (including overlap): " << nTotalBlockedDofsWithOverlap << " / "
            << nTotalFlatDofsWithOverlap << std::endl;

        // Total 'real' number of DOFs (substract missing concentration DOFs on inner membrane elements)
        int nTotalDofsWithoutOverlap = (nNodesXCoarse*nNodesYCoarse * (NUMBER_OF_SPECIES+1))
            - ((params.nMembraneElements()-1) * nNodesXCoarse * NUMBER_OF_SPECIES);
        debug_info << "Total #DOFs of the sequential problem (i.e., without overlap): "
            << nTotalDofsWithoutOverlap << std::endl;
        dofs[i] = nTotalDofsWithoutOverlap;
        debug_info << "===============================" << std::endl;


        Ax1RefinementMapGridFunction<GV,PHYSICS> fineToCoarseElemenIndexGF(gv,physics,nNodesXCoarse-1,
            nNodesYCoarse-1);
        typename Ax1RefinementMapGridFunction<GV,PHYSICS>::Traits::DomainType xlocal(0.5);
        std::vector<int> fineToCoarseMap(gv.size(0));
        std::multimap<int,int> coarseToFineMap;
        for (typename PHYSICS::ElementIterator_All eit=gv.template begin<0>(); eit!=gv.template end<0>(); ++eit)
        {
          int elemIndex = physics.getElementIndex(*eit);
          typename Ax1RefinementMapGridFunction<GV,PHYSICS>::Traits::RangeType coarseIndex;
          fineToCoarseElemenIndexGF.evaluate(*eit,xlocal,coarseIndex);
          fineToCoarseMap[elemIndex] = coarseIndex[0];
          coarseToFineMap.insert(std::make_pair(coarseIndex[0],elemIndex));

          //debug_jochen << "fine elem #" << elemIndex << " <-> " << " coarse elem #" << coarseIndex << std::endl;
        }
        debug_jochen << "Done setting up fineToCoarseMap!" << std::endl;

        // Now loop over coarse grid and set up map with associations between fine element index and
        // coarse entity, which can then be used in a special gridfunction
        typename PHYSICS::ElementMapper ceMapper(gv_coarse);
        std::map<int,typename PHYSICS::ElementPointer> fineToCoarseElementMap;
        for(typename PHYSICS::ElementIterator_All ceit = gv_coarse.template begin<0>();
            ceit != gv_coarse.template end<0>(); ++ceit)
        {
          int ceIndex = ceMapper.map(*ceit);
          //debug_jochen << "Coarse entity #" << ceIndex << std::endl;

          int countFineElements = 0;

          std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator>
            fineElems = coarseToFineMap.equal_range(ceIndex);
          for (std::multimap<int, int>::iterator it = fineElems.first; it != fineElems.second; ++it)
          {
            typename PHYSICS::ElementPointer ep(ceit);
            std::pair<int,typename PHYSICS::ElementPointer> pair(it->second,ep);

            fineToCoarseElementMap.insert(pair);
            countFineElements++;
          }
          //debug_jochen << "#fine elements belonging to coarse entity #" << countFineElements << std::endl;
//          for(int i=0; i<fineToCoarseMap.size(); i++)
//          {
//            if(fineToCoarseMap[i] == ceIndex)
//            {
//              typename PHYSICS::ElementPointer ep(ceit);
//              std::pair<int,typename PHYSICS::ElementPointer> pair(i,ep);
//
//              fineToCoarseElementMap.insert(pair);
//              countFineElements++;
//            }
//          }

          //debug_jochen << "nFine = " << countFineElements << std::endl;
          if(countFineElements == 0)
            DUNE_THROW(Dune::Exception, "Could not assign coarse entity #" << ceIndex << " @"
                << ceit->geometry().center() << " to any fine entity!");
        }
        debug_jochen << " map size: " << fineToCoarseElementMap.size() << std::endl;
        debug_jochen << " nElements: " << physics.nElements() << std::endl;
        assert(fineToCoarseElementMap.size() == physics.nElements());
        debug_jochen << "Done setting up fine to coarse entity map!" << std::endl;

        // Test mapping
//        for (typename PHYSICS::ElementIterator_All eit=gv.template begin<0>(); eit!=gv.template end<0>(); ++eit)
//        {
//          int elemIndex = physics.getElementIndex(*eit);
//
//          typename PHYSICS::ElementPointer cep = fineToCoarseElementMap.at(elemIndex);
//
//          //debug_jochen << "Fine entity #" << elemIndex << " @ " << eit->geometry().center()
//          //    << " -> coarse entity @" << cep->geometry().center() << std::endl;
//        }

        Ax1RefinementMapGridFunction_Subdomain<SubGV,PHYSICS> fineToCoarseElemenIndexGF_elec(elecGV,physics,
            nNodesXCoarse-1,nNodesYCoarse-1-params.nMembraneElements());
        std::vector<int> fineToCoarseMap_elec(elecGV.size(0));
        std::multimap<int,int> coarseToFineMap_elec;
        const bool calculateConcentrationErrors = params.general.get("calculateConcentrationErrors",false);
        if(calculateConcentrationErrors)
        {
          for (typename PHYSICS::SubDomainElementIterator_All seit=elecGV.template begin<0>();
              seit!=elecGV.template end<0>(); ++seit)
          {
            int elemIndex = physics.getSubDomainElementIndex(*seit);
            typename Ax1RefinementMapGridFunction<GV,PHYSICS>::Traits::RangeType coarseIndex;
            fineToCoarseElemenIndexGF_elec.evaluate(*seit,xlocal,coarseIndex);
            fineToCoarseMap_elec[elemIndex] = coarseIndex[0];
            coarseToFineMap_elec.insert(std::make_pair(coarseIndex[0],elemIndex));

            //debug_jochen << "fine elem #" << elemIndex << " <-> " << " coarse elem #" << coarseIndex << std::endl;
          }
          debug_jochen << "Done setting up fineToCoarseMap_elec!" << std::endl;
          debug_jochen << "Size: " << fineToCoarseMap_elec.size() << std::endl;
        }

        int countGlobal = 0;
        int countLast = 0;
        typename PHYSICS::SubDomainElementMapper ceMapper_elec(elecGV_coarse);
        std::map<int,typename PHYSICS::SubDomainElementPointer> fineToCoarseElementMap_elec;
        if(calculateConcentrationErrors)
        {
          for(typename PHYSICS::SubDomainElementIterator_All cseit = elecGV_coarse.template begin<0>();
              cseit != elecGV_coarse.template end<0>(); ++cseit)
          {
            int ceIndex = ceMapper_elec.map(*cseit);
            //debug_jochen << "Coarse entity #" << ceIndex << std::endl;

            int countFineElements = 0;
            std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator>
              fineElems = coarseToFineMap_elec.equal_range(ceIndex);
            for (std::multimap<int, int>::iterator it = fineElems.first; it != fineElems.second; ++it)
            {
              typename PHYSICS::SubDomainElementPointer ep(cseit);
              std::pair<int,typename PHYSICS::SubDomainElementPointer> pair(it->second,ep);

              fineToCoarseElementMap_elec.insert(pair);
              countFineElements++;
            }
//            for(int i=0; i<fineToCoarseMap_elec.size(); i++)
//            {
//              if(fineToCoarseMap_elec[i] == ceIndex)
//              {
//                typename PHYSICS::SubDomainElementPointer ep(cseit);
//                std::pair<int,typename PHYSICS::SubDomainElementPointer> pair(i,ep);
//
//                fineToCoarseElementMap_elec.insert(pair);
//                countFineElements++;
//              }
//            }
            if(countFineElements != countLast)
            {
              debug_jochen << "Coarse element @" << cseit->geometry().center() << ", unusual # of fine "
                  << "elements: " << countFineElements << std::endl;
            }

            countLast = countFineElements;
            countGlobal += countFineElements;

            //debug_jochen << "nFine = " << countFineElements << std::endl;
            if(countFineElements == 0)
              DUNE_THROW(Dune::Exception, "Could not assign elec coarse entity #" << ceIndex << " @"
                  << cseit->geometry().center() << " to any elec fine entity!");
          }

          debug_jochen << " map size: " << fineToCoarseElementMap_elec.size() << std::endl;
          debug_jochen << " nElements: " << elecGV.size(0) << std::endl;
          debug_jochen << " total number of assigned fine elements: " << countGlobal << std::endl;
          assert(fineToCoarseElementMap_elec.size() == elecGV.size(0));
          debug_jochen << "Done setting up fine to elec coarse entity map!" << std::endl;
        }

        // Test mapping
        Dune::FieldVector<double,2> one(1.0);
        Dune::FieldVector<double,2> zero(0.0);
        int elemIndex = -1;
//        for (typename PHYSICS::SubDomainElementIterator_All seit=elecGV.template begin<0>();
//            seit!=elecGV.template end<0>(); ++seit)
//        {
//          elemIndex = physics.getSubDomainElementIndex(*seit);
//
//          typename PHYSICS::SubDomainElementPointer cep = fineToCoarseElementMap_elec.at(elemIndex);
//
//          Dune::FieldVector<double,2> local = cep->geometry().local(seit->geometry().center());
//
//          if(! (Tools::lessOrEqualThan(local,one) && Tools::greaterOrEqualThan(local,zero)))
//          {
//            DUNE_THROW(Dune::Exception, "Coarse element @" << cep->geometry().center()
//                << " does not contain fine element @" << seit->geometry().center() << "!");
//          }
//
//          //debug_jochen << "Elec fine entity #" << elemIndex << " @ " << seit->geometry().center()
//          //    << " -> elec coarse entity @" << cep->geometry().center() << std::endl;
//        }
//        debug_jochen << "Last used elemIndex: " << elemIndex << std::endl;


        std::vector<int> equiMembraneGroups(1, -12);
        bool doInterpolate = false;
        bool useGridConvergenceMode = true;

        // Ok, now put the map that we just set up into a new gridfunction and use it for
        // interpolate the coarse GF onto the fine grid
        typedef Ax1CoarseToFineGridMap<GV,typename OutputTraits::DGF_POT,PHYSICS> GF_POT_COARSE2FINE;
        GF_POT_COARSE2FINE dgfPot_fine(gv,dgfPot_coarse,physics,fineToCoarseElementMap,equiMembraneGroups);

        // All concentrations
        typedef Ax1CoarseToFineGridMap<SubGV,typename OutputTraits::DGF_CON,PHYSICS> GF_CON_COARSE2FINE;
        GF_CON_COARSE2FINE dgfCon_fine(elecGV,dgfCon_coarse,physics,fineToCoarseElementMap_elec,equiMembraneGroups);

        // Sodium only
        //typedef Ax1CoarseToFineGridMap<SubGV,DGF_NA,PHYSICS> GF_CON_COARSE2FINE;
        //GF_CON_COARSE2FINE dgfNa_fine(elecGV,dgfNa_coarse,physics,fineToCoarseElementMap_elec,equiMembraneGroups);
        // Potassium only
        //typedef Ax1CoarseToFineGridMap<SubGV,DGF_K,PHYSICS> GF_CON_COARSE2FINE;
        //GF_CON_COARSE2FINE dgfK_fine(elecGV,dgfK_coarse,physics,fineToCoarseElementMap_elec,equiMembraneGroups);
        // Chloride
        //typedef Ax1CoarseToFineGridMap<SubGV,DGF_CL,PHYSICS> GF_CON_COARSE2FINE;
        //GF_CON_COARSE2FINE dgfCl_fine(elecGV,dgfCl_coarse,physics,fineToCoarseElementMap_elec,equiMembraneGroups);

        // The commented interpolation GFs cannot be used, as they rely on a coarse grid
        // with only one element in x-direction
//        typedef Ax1CoarseToFineGridTransferGridFunction<GV,typename OutputTraits::DGF_POT,PHYSICS>
//          GF_POT_COARSE2FINE;
//        GF_POT_COARSE2FINE dgfPot_fine(gv,dgfPot_coarse,physics,nNodesXCoarse-1,nNodesYCoarse-1,
//           equiMembraneGroups,doInterpolate,useGridConvergenceMode);
//         typedef Ax1CoarseToFineGridTransferGridFunction<SubGV,typename OutputTraits::DGF_CON,PHYSICS>
//           GF_CON_COARSE2FINE;
//         GF_CON_COARSE2FINE dgfCon_fine(elecGV,dgfCon_coarse,physics,nNodesXCoarse-1,nNodesYCoarse-1,
//           equiMembraneGroups,doInterpolate,useGridConvergenceMode);

        // This interpolation GF is horribly slow, as it pretty much iterates over the whole
        // grid each time evaluate() is called, which is a lot!
        // TODO: Write an interpolation GF which does this grid iteration procedure ONCE and saves
        // a fine-to-coarse-entity mapping, probably as a map from cell index (int) to the coarse
        // entity (EntityPointer). This should be doable in reasonable time!
//        typedef Ax1CoarseGridInterpolationGridFunction<GV,typename OutputTraits::DGF_POT,PHYSICS>
//                  GF_POT_COARSE2FINE;
//        GF_POT_COARSE2FINE dgfPot_fine(gv,*dgfPot_coarse[0],physics);


         // Now that we have the interpolation function, create GF for L2 / max error!
        const int integrationOrder = params.general.get("integrationOrder",8);

         // GFS_POT or GFS_POT_SUB?
//         typedef DifferenceSquaredAdapter<GF_POT_COARSE2FINE,typename Traits::DGF_POT> GF_L22_POT;
//         GF_L22_POT l22Pot(dgfPot_fine,dgfPot);
//
//         typename GF_L22_POT::Traits::RangeType l2Error(0.0);
//         const int integrationOrder = 1;
//         Acme2CylGeometryTools::integrateGridFunctionOverCylinder(l22Pot, l2Error,integrationOrder);
//         l2Error = std::sqrt(l2Error);
//         debug_jochen << "L2 error for grid #" << i << ": " << l2Error << std::endl;


         debug_jochen << "--- Starting error calculations in cylinder coordinates..." << std::endl;
         double l2ErrPot = ErrorNorms::l2Norm(dgfPot_fine,dgfPot,integrationOrder);
         debug_jochen << "=== [pot] ErrorNorms::l2Error() for grid #" << i << ": " << l2ErrPot << std::endl;
         double l2ErrCon = -1;
         if(calculateConcentrationErrors)
         {
           // All concentrations
           l2ErrCon = ErrorNorms::l2Norm(dgfCon_fine,dgfConElec,integrationOrder);
           // Sodium only
           //double l2ErrCon = ErrorNorms::l2Norm(dgfNa_fine,dgfNa,integrationOrder);
           //double l2ErrCon = ErrorNorms::l2Norm(dgfK_fine,dgfK,integrationOrder);
           //double l2ErrCon = ErrorNorms::l2Norm(dgfCl_fine,dgfCl,integrationOrder);
           debug_jochen << "=== [con] ErrorNorms::l2Error() for grid #" << i << ": " << l2ErrCon << std::endl;
         }

         double maxErrPot = ErrorNorms::maxNorm(dgfPot_fine,dgfPot,integrationOrder);
         debug_jochen << "=== [pot] ErrorNorms::maxError() for grid #" << i << ": " << maxErrPot << std::endl;

         double maxErrCon = -1;
         if(calculateConcentrationErrors)
         {
           maxErrCon = ErrorNorms::maxNorm(dgfCon_fine,dgfConElec,integrationOrder);
           //double maxErrCon = ErrorNorms::maxNorm(dgfNa_fine,dgfNa,integrationOrder);
           //double maxErrCon = ErrorNorms::maxNorm(dgfK_fine,dgfK,integrationOrder);
           //double maxErrCon = ErrorNorms::maxNorm(dgfCl_fine,dgfCl,integrationOrder);

           debug_jochen << "=== [con] ErrorNorms::maxError() for grid #" << i << ": " << maxErrCon << std::endl;
         }

         if(params.general.get("calculateIn2D",false))
         {
           debug_jochen << std::endl;
           debug_jochen << "--- Starting error calculations in plain 2D coordinates..." << std::endl;
           double l2ErrPot2D = ErrorNorms::l2Norm<GF_POT_COARSE2FINE,typename Traits::DGF_POT,false>
             (dgfPot_fine,dgfPot,integrationOrder);
           debug_jochen << "=== [pot] ErrorNorms::l2Error() for grid #" << i << ": " << l2ErrPot2D << std::endl;
           if(calculateConcentrationErrors)
           {
             double l2ErrCon2D = ErrorNorms::l2Norm<GF_CON_COARSE2FINE,typename Traits::DGF_CON_SUB,false>
               (dgfCon_fine,dgfConElec,integrationOrder);
             //double l2ErrCon2D = ErrorNorms::l2Norm<GF_CON_COARSE2FINE,DGF_NA,false>
             //   (dgfNa_fine,dgfNa,integrationOrder);
             //double l2ErrCon2D = ErrorNorms::l2Norm<GF_CON_COARSE2FINE,DGF_K,false>
             //    (dgfK_fine,dgfK,integrationOrder);
             //double l2ErrCon2D = ErrorNorms::l2Norm<GF_CON_COARSE2FINE,DGF_CL,false>
             //     (dgfCl_fine,dgfCl,integrationOrder);
             debug_jochen << "=== [con] ErrorNorms::l2Error() for grid #" << i << ": " << l2ErrCon2D << std::endl;
             // (Max errors are of course the same as in cylinder coordinates, as no integration is done)
           }
         }

         l2ErrorsPot[i] = l2ErrPot;
         l2ErrorsCon[i] = l2ErrCon;
         maxErrorsPot[i] = maxErrPot;
         maxErrorsCon[i] = maxErrCon;
      }
      debug_info << std::endl << std::endl << std::endl;
      debug_info << " ============================================================= " << std::endl;
      debug_info << " SUMMARY" << std::endl;
      debug_info << " ============================================================= " << std::endl;

      for(int i=1; i<grids.size(); i++)
      {
        debug_jochen << "Grids " << (i-1) << " - " << i << ":" << std::endl;
        debug_jochen << "  - DOFs: " << dofs[i-1] << " - " << dofs[i] << " | factor "
            << (dofs[i]/dofs[i-1]) << std::endl;
        debug_jochen << "  - [pot] L2 errors: " << l2ErrorsPot[i-1] << " - " << l2ErrorsPot[i] << " | factor "
            << (l2ErrorsPot[i-1]/l2ErrorsPot[i]) << std::endl;
        debug_jochen << "  - [pot] max errors: " << maxErrorsPot[i-1] << " - " << maxErrorsPot[i] << " | factor "
            << (maxErrorsPot[i-1]/maxErrorsPot[i]) << std::endl;
        debug_jochen << "  - [con] L2 errors: " << l2ErrorsCon[i-1] << " - " << l2ErrorsCon[i] << " | factor "
                      << (l2ErrorsCon[i-1]/l2ErrorsCon[i]) << std::endl;
        debug_jochen << "  - [con] max errors: " << maxErrorsCon[i-1] << " - " << maxErrorsCon[i] << " | factor "
            << (maxErrorsCon[i-1]/maxErrorsCon[i]) << std::endl;
        double conv_order_l2_pot = (std::log2(l2ErrorsPot[i-1])-std::log2(l2ErrorsPot[i]));
        /// (std::log(dofs[i-1])-std::log(dofs[i]));
        double conv_order_l2_con = (std::log2(l2ErrorsCon[i-1])-std::log2(l2ErrorsCon[i]));
        double conv_order_max_pot = (std::log2(maxErrorsPot[i-1])-std::log2(maxErrorsPot[i]));
        /// (std::log(dofs[i-1])-std::log(dofs[i]));
        double conv_order_max_con = (std::log2(maxErrorsCon[i-1])-std::log2(maxErrorsCon[i]));

        debug_jochen << "  - [pot] convergence order (L2 norm): " << conv_order_l2_pot << std::endl;
        debug_jochen << "  - [pot] convergence order (max norm): " << conv_order_max_pot << std::endl;
        debug_jochen << "  - [con] convergence order (L2 norm): " << conv_order_l2_con << std::endl;
        debug_jochen << "  - [con] convergence order (max norm): " << conv_order_max_con << std::endl;
        debug_jochen << std::endl;
      }
      // =============================================================================================


      // ======= Save simulation state ================================================================
      if(physics.getParams().doSaveState())
      {
        std::string saveFilename = physics.getParams().getSaveFilename();
        acme2_cylOutput.saveState(time,dt,saveFilename);
        debug_info << "=============================================================" << std::endl;
        debug_info << "Saved simulation state to file " << saveFilename << std::endl;
        debug_info << "=============================================================" << std::endl;
      }
      // ==============================================================================================

      // Print #DOFs at simulation end
      debug_info << "U.N(): " << unew.N() << std::endl;
      debug_info << "U.flatsize(): " << unew.flatsize() << std::endl;
      debug_info << "con/pot GFS size: " << gfsConSub.size() << " / " << gfsPotSub.size() << std::endl;

      debug_info << "Total Acme2CylOutput time: " << acme2_cylOutput.getTotalOutputTime() << "s" << std::endl;
    }

  private:
    std::vector<Dune::shared_ptr<Grid> >& grids;

    GV& gv;
    PHYSICS& physics;

    SubGV& elecGV;
    SubGV& membGV;
};

#endif /* DUNE_AX1_ACME2CYL_SETUP_HH */
