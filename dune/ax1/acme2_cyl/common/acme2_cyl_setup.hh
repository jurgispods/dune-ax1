#ifndef DUNE_AX1_ACME2CYL_SETUP_HH
#define DUNE_AX1_ACME2CYL_SETUP_HH

#include <dune/common/array.hh>

#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/common/clock.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/conforming.hh>
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
#include <dune/pdelab/ordering/permutedordering.hh>
#include <dune/pdelab/ordering/permutationordering.hh>

#include <dune/ax1/common/ax1_discretegridfunction.hh>
#include <dune/ax1/common/ax1_linearproblem.hh>
#include <dune/ax1/common/ax1_newton.hh>
#include <dune/ax1/common/ax1_parallelconstraintshelper.hh>
#include <dune/ax1/common/ax1_simulationdata.hh>
#include <dune/ax1/common/ax1_simulationstate.hh>
#include <dune/ax1/common/ax1_solution_container.hh>
#include <dune/ax1/common/ax1_solverbackend.hh>
#include <dune/ax1/common/chargedensitygridfunction.hh>
#include <dune/ax1/common/charge_layer_mori.hh>
#include <dune/ax1/common/charge_layer_mori_implicit.hh>
#include <dune/ax1/common/ionicstrengthgridfunction.hh>
#include <dune/ax1/common/membranefluxgridfunction.hh>
#include <dune/ax1/common/membranefluxgridfunction_implicit.hh>
#include <dune/ax1/common/poisson_boltzmann_concentrationgridfunction.hh>
#include <dune/ax1/common/poisson_boltzmann_rhs_gridfunction.hh>
#include <dune/ax1/common/power_parameters.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/vectordiscretegridfunctiongradient.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_boundary.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_geometrytools.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_initial.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_output.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_solutionvectors.hh>
#include <dune/ax1/acme2_cyl/common/poisson_boltzmann_parameters.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_simulation.hh>
#include <dune/ax1/acme2_cyl/configurations/default/default_nernst_planck_parameters.hh>
#include <dune/ax1/acme2_cyl/configurations/default/default_poisson_parameters.hh>
#include <dune/ax1/acme2_cyl/operator/convectiondiffusionfem.hh>
#include <dune/ax1/acme2_cyl/operator/poisson_boltzmann_operator.hh>
#include <dune/ax1/acme2_cyl/operator/acme2_cyl_operator_fully_coupled.hh>
#include <dune/ax1/acme2_cyl/operator/acme2_cyl_toperator.hh>

template<class Grid, class GV, class PHYSICS, class SubGV>
class Acme2CylSetup
{
  public:

    template<typename DGF_CON, typename DGF_POT, typename SolutionContainer,bool implicit>
    struct MembraneFluxClassSelector
    {
      typedef MembraneFluxGridFunction<DGF_CON,DGF_POT,PHYSICS,SolutionContainer> GF_MEMB_FLUX;
      typedef ChargeLayerMori<DGF_CON,DGF_POT,PHYSICS,SolutionContainer> GF_MORI_FLUX;
    };

    template<typename DGF_CON, typename DGF_POT, typename SolutionContainer>
    struct MembraneFluxClassSelector<DGF_CON,DGF_POT,SolutionContainer,true>
    {
      typedef MembraneFluxGridFunctionImplicit<DGF_CON,DGF_POT,PHYSICS,SolutionContainer> GF_MEMB_FLUX;
      typedef ChargeLayerMoriImplicit<DGF_CON,DGF_POT,PHYSICS,SolutionContainer> GF_MORI_FLUX;
    };

    //! Traits class providing function spaces and related data types for the acme2_cyl use case
    template<int NUMBER_OF_SPECIES>
    struct Acme2CylTraits
    {
      typedef Acme2CylSetup<Grid,GV,PHYSICS,SubGV> SETUP;

      typedef PHYSICS Physics;

      typedef GV GridView;
      typedef Grid GridType;
      typedef typename GV::Grid::ctype Coord;
      typedef double Real;

      typedef SubGV SubGridView;
      typedef typename SubGV::Grid SubGridType;

#if USE_PARALLEL==1
#if USE_OVERLAP==1
      //typedef Dune::PDELab::OverlappingConformingDirichletConstraints CONSTRAINTS;
      typedef Ax1OverlappingDirichletContraints CONSTRAINTS;
#else
      typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints<GV> CONSTRAINTS;
#endif
#else
      typedef Dune::PDELab::ConformingDirichletConstraints CONSTRAINTS;
#endif
      typedef Dune::PDELab::ConformingDirichletConstraints SEQUENTIAL_CONSTRAINTS;
      typedef Dune::PDELab::ISTLVectorBackend<> VBE_SINGLE;
      typedef Dune::PDELab::ISTLVectorBackend<> VBE_CON;

#ifdef MULTIPLE_MEMBRANE_ELEMENTS
      typedef Dune::PDELab::ISTLVectorBackend<> VBE;
#else
      //typedef Dune::PDELab::ISTLVectorBackend<> VBE;
      typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::static_blocking,AX1_BLOCKSIZE> VBE;
#endif

      const static std::size_t feOrder = 1;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,Real,feOrder> FEM_CON;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,Real,feOrder> FEM_POT;

      // Nernst-Planck GFS (on electrolyte subdomain)
      typedef Dune::PDELab::GridFunctionSpace<SubGV,FEM_CON,CONSTRAINTS,VBE_SINGLE> GFS_SINGLE_CON;

      typedef Dune::PDELab::PowerGridFunctionSpace<
        GFS_SINGLE_CON, NUMBER_OF_SPECIES, VBE_CON, Dune::PDELab::EntityBlockedOrderingTag> GFS_CON;

      // Poisson GFS (on electrolyte or membrane subdomain)
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM_POT,CONSTRAINTS,VBE_SINGLE> GFS_POT;

      // Choose ordering based on chosen setup
#ifdef MULTIPLE_MEMBRANE_ELEMENTS
      //typedef Dune::PDELab::LexicographicOrderingTag OrderingTag
      typedef Dune::PDELab::ordering::Permuted<Dune::PDELab::LexicographicOrderingTag> OrderingTag;
#else
      typedef Dune::PDELab::InterleavedOrderingTag OrderingTag;
      // Fallback ordering for cases where InterleavedOrderingTag will not work
      //typedef Dune::PDELab::LexicographicOrderingTag OrderingTag;
#endif
      // Multidomain GFS
      typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<GridType,VBE,OrderingTag,
        GFS_CON,GFS_POT> MultiGFS;

      // Extract SubGFS for use in gridfunctions; use these after creating MultiGFS!
      typedef typename Dune::PDELab::GridFunctionSubSpace<MultiGFS,Dune::TypeTree::TreePath<0> > GFS_CON_SUB;
      typedef typename Dune::PDELab::GridFunctionSubSpace<MultiGFS,Dune::TypeTree::TreePath<1> > GFS_POT_SUB;

      // Extract solution vector type
      //typedef typename Dune::PDELab::BackendVectorSelector<GFS_CON,Real>::Type U_CON;
      //typedef typename Dune::PDELab::BackendVectorSelector<GFS_POT,Real>::Type U_POT;
      typedef typename Dune::PDELab::BackendVectorSelector<MultiGFS,Real>::Type U;

      // Helper type: Newton solution container
      typedef Ax1SolutionContainer<U, GFS_POT_SUB, GFS_CON_SUB> SolutionContainer;

      // Various gridfunctions
      typedef Dune::PDELab::VectorDiscreteGridFunction <GFS_CON_SUB,U> DGF_CON_SUB;
      typedef Dune::PDELab::Ax1VectorDiscreteGridFunctionGradient <GFS_CON_SUB,U> DGF_CON_GRAD_SUB;

      typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_CON_SUB, PHYSICS> DGF_CON;
      typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_CON_GRAD_SUB, PHYSICS> DGF_CON_GRAD;

      typedef Dune::PDELab::DiscreteGridFunction<GFS_POT_SUB,U> DGF_POT;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_POT_SUB,U> DGF_POT_GRAD;

      typedef ChargeDensityGridFunction<DGF_CON, PHYSICS> GF_CD;
      typedef IonicStrengthGridFunction<DGF_CON, PHYSICS> GF_IS;

      typedef typename MembraneFluxClassSelector<DGF_CON,DGF_POT,SolutionContainer,
          PHYSICS::Traits::useImplicitMembraneFlux>::GF_MEMB_FLUX GF_MEMB_FLUX;
      typedef typename MembraneFluxClassSelector<DGF_CON,DGF_POT,SolutionContainer,
          PHYSICS::Traits::useImplicitMembraneFlux>::GF_MORI_FLUX GF_MORI_FLUX;

      // Custom grid function for initial concentration values
      typedef typename PHYSICS::Traits::INITIAL_CON INITIAL_GF_CON; // not obsolete
      typedef typename PHYSICS::Traits::INITIAL_POT INITIAL_GF_POT; // obsolete for now

      // Parameter classes
      typedef typename PHYSICS::Traits::template NERNST_PLANCK_PARAMETERS<PHYSICS,GF_MEMB_FLUX,GF_MORI_FLUX> PARAMETERS_SINGLE_CON;
      //typedef std::vector<Dune::shared_ptr<PARAMETERS_SINGLE_CON> > PARAMETERS_CON;
      typedef PowerParameters<PARAMETERS_SINGLE_CON,NUMBER_OF_SPECIES> PARAMETERS_CON;
      typedef typename PHYSICS::Traits::template POISSON_PARAMETERS<PHYSICS> PARAMETERS_POT;

      // Generic grid functions for initial values (subdomains)
      //typedef InitialCon<GV,Real,NUMBER_OF_SPECIES,PHYSICS,INITIAL_GF_CON> INITIAL_CON;
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PARAMETERS_SINGLE_CON> INITIAL_SINGLE_CON;
      typedef Dune::PDELab::PowerGridFunction<INITIAL_SINGLE_CON, NUMBER_OF_SPECIES> INITIAL_CON;

      //typedef InitialPot<GV,Real,PHYSICS,INITIAL_GF_POT> INITIAL_POT;
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PARAMETERS_POT> INITIAL_POT;

      typedef Dune::PDELab::CompositeGridFunction<INITIAL_CON,INITIAL_POT> INITIAL_ELEC;

      // Helper grid function to calculate equilibirum concentrations from stationary Poisson-Boltzmann potential
      //typedef PoissonBoltzmannRHSGridFunction<DGF_POT, INITIAL_CON, PHYSICS> GF_PB_RHS;
      //typedef PoissonBoltzmannConcentrationGridFunction<INITIAL_CON, DGF_POT, SubGV, PHYSICS> GF_PB_CON;

      typedef typename PHYSICS::Traits::template ELEC_OPERATOR<PARAMETERS_CON,PARAMETERS_POT,FEM_CON,FEM_POT,
          PHYSICS::Traits::useMembraneContributions> ELEC_OPERATOR;

    };

    // Choose domain and range field type
    typedef GV GridView;
    typedef typename GV::Grid GridType;
    typedef typename GV::Grid::ctype Coord;
    typedef double Real;

    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

    // SubGrid stuff
    typedef typename SubGV::Grid SubGrid;
    typedef Acme2CylTraits<NUMBER_OF_SPECIES> Traits;

    static const bool writeIfNotConverged = true;

    Acme2CylSetup(Grid& grid_, GV& gv_, PHYSICS& physics_,
        SubGV& elecGV_, SubGV& membGV_)
    : grid(grid_),
      gv(gv_),
      physics(physics_),
      elecGV(elecGV_),
      membGV(membGV_)
    {
      if(physics.getParams().general.get("useMoriOperatorSplit",false))
        DUNE_THROW(Dune::Exception, "This application uses fully-coupled Newton scheme, but the "
            << "incompatible flag 'useMoriOperatorSplit' is set!");
    }

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
      typename Traits::FEM_CON femCon(gv);
      typename Traits::FEM_POT femPot(gv);
      const int lfsPotSize = femPot.maxLocalSize();
      debug_jochen << "lfsPotSize = " << lfsPotSize << std::endl;

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

      // Create MultiGFS with permutation ordering (identity for now)
      typename Traits::OrderingTag permutationOrderingTag;
      typename Traits::MultiGFS multigfs(grid,permutationOrderingTag,gfsCon,gfsPot);

      if(! params.doReorderMatrix())
      {
        DUNE_THROW(Dune::Exception, "It is necessary to use matrix reordering when using more than one membrane element layer!");
      }

      if(params.doReorderMatrix())
      {
        debug_info << " === Reordering matrix! === " << std::endl;
        multigfs.update();

        // Permutation for DOFs
        std::vector<std::size_t> permutation, inv_permutation;

        // Set initial permutation to identity
        permutation.resize(gfsCon.size()+gfsPot.size());
        inv_permutation.resize(gfsCon.size()+gfsPot.size());

        // Reorder matrix by grouping together DOFs at vertices
        // Old method (error-prone!)
  //      if(params.doReorderMatrix())
  //      {
  //        Tools::setupDOFPermutation(physics, gfsCon.size(), gfsPot.size(), permutation, inv_permutation);
  //      }

        debug_info << " == Setting up DOF permutation vector... " << std::endl;
        // Use multigfs object to set up the real permutation and update in ordering
        physics.setupDOFPermutation(multigfs, permutation, inv_permutation);

        // Old method when using PermutationOrderingTag
        //multigfs.orderingTag().updatePermutation(inv_permutation);

        debug_info << " == Handing over permutation vector to MultiGFS ordering... " << std::endl;
        // Following: new method using Steffen's PermutedOrderingTag
        std::vector<std::size_t>& tagPerm = multigfs.orderingTag().permutation();
        if(tagPerm.size() != inv_permutation.size())
        {
          debug_verb << "[MME permutation] Permuted ordering vector has to be resized from size " << tagPerm.size()
              << " to size " << inv_permutation.size() << std::endl;
          tagPerm.resize(inv_permutation.size());
        }
        debug_verb << "[MME permutation] Using the following permutation:" << std::endl;
        for(int i=0; i<tagPerm.size(); ++i)
        {
          tagPerm[i] = inv_permutation[i];
          debug_verb << i << " -> " << inv_permutation[i] << std::endl;
        }
        debug_verb << "[MME permutation] Updated ordering permutation!" << std::endl;
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
      typename Traits::MultiGFS multigfs(grid,interleavedOrderingTag,gfsCon,gfsPot);

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


      // ========== A bunch of gridfunctions ==========================================================
      typename Traits::DGF_CON_SUB dgfConElec(gfsConSub,unew);
      typename Traits::DGF_CON_SUB dgfOldConElec(gfsConSub,uold);
      typename Traits::DGF_CON_GRAD_SUB dgfGradConElec(gfsConSub, unew);

      typename Traits::DGF_CON dgfCon(gv, dgfConElec, physics);
      typename Traits::DGF_CON dgfOldCon(gv, dgfOldConElec, physics);
      typename Traits::DGF_CON_GRAD dgfGradCon(gv, dgfGradConElec, physics);

      typename Traits::GF_CD gfChargeDensity(dgfCon, dgfOldCon, physics);
      typename Traits::GF_IS gfIonicStrength(dgfCon, physics);

      typename Traits::DGF_POT dgfPot(gfsPotSub, unew);
      typename Traits::DGF_POT dgfOldPot(gfsPotSub, uold);
      typename Traits::DGF_POT_GRAD dgfGradPot(gfsPotSub, unew);
      // ==============================================================================================


      // ============== Future feature: fully fully-implicit method ===================================
      /* For al truly fully-implicit method, the membrane boundary conditions need to be updated in each
       * Newton iteration. For this, we need the currect Newton solution vector, as the evolution of
       * the membrane flux depends on the membrane potential difference, which is a global informatio
       * we cannot access from the local operator, where each element only sees its own local solution
       * vector. Therefore, as a hack, create a container which contains a pointer to the Newton solution
       * vector. The Newton will update this pointer at the beginning of each iteration; the local operator
       * can use this vector to calculate the global measures like membrane potential or membrane interface
       * concentrations.
       * More specifically, we put a reference to this object into the membrane flux gridfunction class
       * instead of the local operator for now, which should be sufficient.
       */
      typename Traits::SolutionContainer newtonSolution(gfsPotSub, gfsConSub);
      // Set initial pointer
      newtonSolution.setSolutionConNew(&unew);
      newtonSolution.setSolutionPotNew(&unew);
      newtonSolution.setSolutionConOld(&uold);
      newtonSolution.setSolutionPotOld(&uold);
      // ==============================================================================================


      // ========== Define parameter classes containing the model problem =============================
      // Custom initial gridfunction fpor concentrations
      typename Traits::INITIAL_GF_CON initialGFCon(gv,physics.getParams());

      // Central class: Membrane flux
      typename Traits::GF_MEMB_FLUX gfMembFlux(dgfOldCon, dgfOldPot, physics, newtonSolution, tEquilibrium);

      std::string loadBoundaryLocation = params.boundary.get("loadBoundary","bottom");
      bool loadMembraneFlux = !physics.getParams().isBoundaryDirichlet_Concentration(loadBoundaryLocation)
              && (loadBoundaryLocation == "membrane")
              && physics.getParams().boundary.get("useTimeDependentBoundaryValuesCon",false);

      // Deactivate membrane flux calculation completely when it is loaded from an external file
      if(loadMembraneFlux)
      {
        gfMembFlux.setActive(false);
      }

      // New Mori flux class
      typename Traits::GF_MORI_FLUX gfMoriFlux(dgfOldCon, dgfOldPot, physics, newtonSolution, tEquilibrium);

      // Hack: Calculate Mori flux without actually using it. Useful for generating flux data that shall be used
      // as a boundary condition in a later simulation
      if(params.general.get("forceMoriFluxCalculation", false))
      {
        gfMoriFlux.setActive(true);
      }
      // Another hack: Do not calculate Mori flux even when solving Mori equations.
      // Useful when loading complete (ionic+Mori) fluxed from previous simulation
      if(params.general.get("disableMoriFluxCalculation", false))
      {
        gfMoriFlux.setActive(false);
      }

      debug_verb << "Initializing parameter classes..." << std::endl;

      // New handling: Vector of parameter classes for better compliance with convection diffusion
      // parameter classes interface
      typename Traits::PARAMETERS_CON parametersCon;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        Dune::shared_ptr<typename Traits::PARAMETERS_SINGLE_CON> ptr(
          new typename Traits::PARAMETERS_SINGLE_CON(j,gv,physics,initialGFCon,gfMembFlux,gfMoriFlux,tEquilibrium));
        parametersCon.push_back(ptr);
      }

      for(int jj=0; jj<parametersCon.size(); jj++)
      {
        debug_verb << "Count " << jj << ": " << parametersCon[jj].use_count() << std::endl;
      }
      debug_verb << "- Initialized concentration parameters." << std::endl;

      // default: PC = e * e * N_A * LS * LS / ( eps0 * k * T );
      Real poissonConstant = physics.getPoissonConstant();
      const Real lengthScale =  physics.getLengthScale();
      const Real timeScale = physics.getTimeScale();
      if(physics.getParams().useMori())
      {
        // Mori: PC = e * e * N_A * LS * LS * LS / (k * T * TS);
        //poissonConstant *= (con_eps0 * LS / TS);

        // Update: Actually the factor e/(kT) should be thrown out in order to match the left hand side units!
        // => PC = e * N_A * LS * LS * LS / TS;
        poissonConstant = con_e * con_mol * lengthScale * lengthScale * lengthScale / timeScale;
      }
      typename Traits::PARAMETERS_POT parametersPot(gv,physics,poissonConstant);

      debug_verb << "- Initialized potential parameters." << std::endl;
      debug_verb << "DONE initializing parameter classes!" << std::endl;

      // In case of loading boundary values from external file: Load initial values once at the beginning
      debug_info << "# Calling initial prepareNextTimeStep(" << (time) << ") on parameter classes" << std::endl;
      parametersPot.prepareNextTimeStep(time);
      for(int j=0; j<parametersCon.size(); ++j)
      {
        parametersCon[j]->prepareNextTimeStep(time);
      }
      // ==============================================================================================


      // ====== Create grid functions and interpolate initial values onto coefficient vectors =========
      // Generic wrapper class for initial grid function

      // Initial concentrations
      std::vector<Dune::shared_ptr<typename Traits::INITIAL_SINGLE_CON> > initialSingleCon;
      typename Traits::INITIAL_CON initialCon;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        Dune::shared_ptr<typename Traits::INITIAL_SINGLE_CON> ptr(new typename Traits::INITIAL_SINGLE_CON
            (gv,*parametersCon[j]));
        initialSingleCon.push_back(ptr);

        initialCon.setChild(j, *initialSingleCon[j]);
      }
      initialCon.setTime(time);

      // Initial potential
      typename Traits::INITIAL_POT initialPot(gv,parametersPot);
      initialPot.setTime(time);

      // Composite initial grid function for the electrolyte subdomain
      typename Traits::INITIAL_ELEC initialElec(initialCon,initialPot);
      // ==============================================================================================

      // ==============================================================================================
      // Define and instantiate output class
      typedef Acme2CylOutput<Traits> ACME2CYL_OUTPUT;
      ACME2CYL_OUTPUT acme2_cylOutput(gv,elecGV,membGV,multigfs,gfsConSub,gfsPotSub,solutionVectors,
          gfMembFlux,gfMoriFlux,physics,parametersPot,parametersCon,2*degreePot,2*degreeCon);
      // ==============================================================================================

      // ============== BOUNDARY CONDITION TYPES ======================================================
      //typedef BCTypeSingleCon<typename Traits::PARAMETERS_SINGLE_CON,PHYSICS> BCType_SINGLE_CON; // old
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<typename Traits::PARAMETERS_SINGLE_CON>
        BCType_SINGLE_CON;

      typedef Dune::PDELab::PowerConstraintsParameters<BCType_SINGLE_CON,NUMBER_OF_SPECIES> BCType_CON;
      typedef BCTypePot<typename Traits::PARAMETERS_POT,PHYSICS> BCType_POT;

      typedef Dune::PDELab::CompositeConstraintsParameters<BCType_CON, BCType_POT> BCType_ELEC;

      // Create NUMBER_OF_SPECIES boundary condition type classes
      BCType_CON bctypeCon;
      for(int i=0; i<NUMBER_OF_SPECIES; ++i)
      {
        BCType_SINGLE_CON* bctypeSingleCon = new BCType_SINGLE_CON(*parametersCon[i]);
        bctypeCon.setChild(i,*bctypeSingleCon);
      }
      BCType_POT bctypePot(gv,parametersPot,physics);

      BCType_ELEC bctypeElec(bctypeCon, bctypePot);
      // ==============================================================================================


      // ============== Define local operators ========================================================
      // ##### Spatial part #####
      // ### Standard FEM (CG) fully-coupled Poisson-Nernst-Planck local operator on electrolyte subdomain
      typedef typename Traits::ELEC_OPERATOR LOP_ELEC;
      //const int intorderadd = (USE_CYLINDER_COORDINATES ? 2 : 0);
      const int intorderadd = params.general.get("quadratureOrderAdd",0);
      LOP_ELEC lopElec(parametersCon, parametersPot, lfsPotSize, intorderadd);

      // ### Standard FEM (CG) Poisson local operator on membrane subdomain
      typedef Dune::PDELab::ConvectionDiffusionFEM<typename Traits::PARAMETERS_POT,typename Traits::FEM_POT,
          PHYSICS::Traits::useMembraneContributions> LOP_MEMB;
      LOP_MEMB lopMemb(parametersPot, lfsPotSize, intorderadd);

      // ##### Temporal part #####
      typedef NernstPlanckTimeLocalOperator<PHYSICS,typename Traits::FEM_CON> TLOP_ELEC;
      TLOP_ELEC tlop(physics, intorderadd);
      // ==============================================================================================


      // =========== dune-multidomain setup stuff =====================================================
      typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> EC;
      typedef Dune::PDELab::MultiDomain::SubDomainSubsetCondition<Grid> SC;
      EC conditionElec(0); //only on subdomain 0
      EC conditionMemb(1); //only on subdomain 1
      SC conditionAll;    //on all subdomains

      // Empty constraints for subproblems; what is this good for anyway!?
#if USE_PARALLEL == 1 && USE_OVERLAP == 0
      typename Traits::CONSTRAINTS empty_constraints(gv);
#else
      typename Traits::CONSTRAINTS empty_constraints;
#endif

      typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<typename Traits::MultiGFS,typename Traits::MultiGFS,
                LOP_ELEC,EC,typename Traits::GFS_CON,typename Traits::GFS_POT> ElecSubProblem;
      ElecSubProblem elecSubProblem(lopElec,conditionElec);

      typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<typename Traits::MultiGFS,typename Traits::MultiGFS,
          LOP_MEMB,EC,typename Traits::GFS_POT> MembSubProblem;
      MembSubProblem membSubProblem(lopMemb,conditionMemb);

      typedef Dune::PDELab::MultiDomain::SubProblem<typename Traits::MultiGFS,typename Traits::MultiGFS,
          TLOP_ELEC,EC,0> TimeSubProblem;
      TimeSubProblem timeSubProblem(tlop,conditionElec);
      // ==============================================================================================


      // =================== BOUNDARY CONDITIONS ======================================================
      typedef typename Traits::MultiGFS::template ConstraintsContainer<Real>::Type CC;
      CC cc;

      auto md_constraints = Dune::PDELab::MultiDomain::template constraints<Real>(multigfs,
            Dune::PDELab::MultiDomain::constrainSubProblem(elecSubProblem, bctypeElec),
            Dune::PDELab::MultiDomain::constrainSubProblem(membSubProblem, bctypePot));

      md_constraints.assemble(cc);
      debug_info << multigfs.size() << " DOF, " << cc.size() << " restricted" << std::endl;

      // Create one additional contraints container, which is only used for the very special case of
      // interpolating time-dependet Dirichlet boundary values in the overlapping parallel case.
      CC ccWithoutOverlap;
      // For the assembly of this container, set the flag in the helper class to false; this is possible
      // because the contraints object is handed to the GFS by reference
      constraints.setOverlapIsDirichlet(false);
      md_constraints.assemble(ccWithoutOverlap);
      // Restore default value of flag in the helper class (don't know if this is necessary)
      constraints.setOverlapIsDirichlet(true);
      debug_info << "Without overlap: " << ccWithoutOverlap.size() << " restricted" << std::endl;

      // The following boundary value classes are not used, as they are assumed to be fulfilled by the
      // initial conditions are not allowed to change over time!
      // ###### Dirichlet values ########
      typedef DirichletValuesSingleCon<typename Traits::PARAMETERS_CON> DirichletValuesSingleCon;

      // Create NUMBER_OF_SPECIES Dirichlet value classes
      typedef Dune::PDELab::PowerGridFunction<DirichletValuesSingleCon,NUMBER_OF_SPECIES> DirichletValuesCon;
      DirichletValuesCon dirichletValuesCon;

      for(int i=0; i<NUMBER_OF_SPECIES; ++i)
      {
        DirichletValuesSingleCon* dirichletValuesSingleCon
          = new DirichletValuesSingleCon(gv,parametersCon,i);
        dirichletValuesCon.setChild(i,*dirichletValuesSingleCon);
      }
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter
          <typename Traits::PARAMETERS_POT> DirichletValuesPot;
      DirichletValuesPot dirichletValuesPot(gv,parametersPot);

      // Is this necessary? Initial conditions should fulfill the boundary conditions anyways,
      // otherwise we have an ill-posed problem!
      //Dune::PDELab::interpolate(dirichletValuesCon_Inside,subGfsCon_Inside,unewCon_Inside);
      //Dune::PDELab::copy_nonconstrained_dofs(ccCon_Inside,uoldCon_Inside,unewCon_Inside);
      //Dune::PDELab::interpolate(dirichletValuesPot,gfsPot,uPot);
      //Output::printSingleCoefficientVector(uPot, "uPot t=0");
      // ==============================================================================================

      // ============== Make grid operator ============================================================
      //typedef Dune::PDELab::ISTLMatrixBackend MBE;
      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
      // Use 9 entries per line for standard Q1 finite element stencil
      MBE mbe(9);

      typedef Dune::PDELab::MultiDomain::GridOperator<typename Traits::MultiGFS,typename Traits::MultiGFS,
        MBE,Real,Real,Real,CC,CC,ElecSubProblem,MembSubProblem> GO0;
      typedef Dune::PDELab::MultiDomain::GridOperator<typename Traits::MultiGFS,typename Traits::MultiGFS,
        MBE,Real,Real,Real,CC,CC,TimeSubProblem> GO1;
      typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;

      GO0 go0(multigfs,multigfs,cc,cc,mbe,elecSubProblem,membSubProblem);
      GO1 go1(multigfs,multigfs,cc,cc,mbe,timeSubProblem);
      IGO igo(go0,go1);

      /* ATTENTION! Do _not_ try to instantiate Jacobian matrix here! The one step gridoperator needs the
       * time stepping method (use igo.setMethod(myTimeStepper) before any pattern generation or assembly can be done!
       */

      // Here, initial values coming from the both subdomain initial GFs are written to uold
      Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,uold,initialElec,elecSubProblem);
      // Only interpolate initial values on membrane subdomain when there exist pure interior membrane vertices
      if(params.nMembraneElements() > 1)
      {
        Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,uold,initialPot,membSubProblem);
      }
      unew = uold;

      //Output::printRawCoefficientVector(unew, "unew");
      // ==============================================================================================

      // ========== Select a linear solver backend ====================================================
      const int maxLinIt = params.general.get("linearSolverMaxIt", 5000);
#if USE_PARALLEL==1
#if USE_OVERLAP==1
      assert(gv.grid().overlapSize(0) > 0);

#ifdef MULTIPLE_MEMBRANE_ELEMENTS
//      // Old method when doing manually reordering of DOFs
//      typedef Dune::PDELab::Ax1ReorderedISTLBackend_OVLP_BCGS_ILU0<typename Traits::MultiGFS,CC> LS;
//      // Hand over DOF index permutations to this backend!
//      LS ls(multigfs,cc,perm,inv_perm,maxLinIt,1);


      // AMG does not make sense with a non-blocked matrix, so it is omitted here!
      //typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILU0<typename Traits::MultiGFS,CC> LS;
      //LS ls(multigfs,cc,maxLinIt,5);
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILUn<typename Traits::MultiGFS,CC> LS;
      LS ls(multigfs,cc,1,maxLinIt,5);
      //typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<typename Traits::MultiGFS,CC> LS;
      //LS ls(multigfs,cc,maxLinIt,5);
      //typedef Dune::PDELab::ISTLBackend_OVLP_GMRES_ILU0<typename Traits::MultiGFS,CC> LS;
      //LS ls(multigfs,cc,maxLinIt,5,40,true);
#else
      //typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILU0<typename Traits::MultiGFS,CC> LS;
      //LS ls(multigfs,cc,maxLinIt,5);
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILUn<typename Traits::MultiGFS,CC> LS;
      LS ls(multigfs,cc,1,maxLinIt,5);
      //typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<typename Traits::MultiGFS,CC> LS;
      //LS ls(multigfs,cc,maxLinIt,5);
      //typedef Dune::PDELab::ISTLBackend_OVLP_GMRES_ILU0<typename Traits::MultiGFS,CC> LS;
      //LS ls(multigfs,cc,maxLinIt,5,40,true);
      //typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<IGO> LS;
      //LS ls(multigfs,maxLinIt,5);
      // This is the preferred preconditioner/linear solver combination for 'hard' parallel problems
      //typedef Dune::PDELab::ISTLBackend_GMRES_AMG_ILU0<IGO> LS;
      //LS ls(multigfs,maxLinIt,5);

//      typename LS::Parameters amgParams(ls.parameters());
//      amgParams.setAccumulate(Dune::Amg::atOnceAccu);
//      amgParams.setDefaultValuesAnisotropic(Grid::dimension);
//      amgParams.setSkipIsolated(true);
//      ls.setParameters(amgParams);
//      Tools::printAmgParams(ls.parameters());
#endif
#else
      assert(gv.grid().overlapSize(0) == 0);
      //typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<IGO> LS;
      //LS ls(multigfs);
      typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_ILUn<IGO> LS;
      LS ls(igo,maxLinIt,1,1);
#endif
#else
      // funktionieren bei: PowerGFS lexicographic / PowerGFS blockwise<1,1,1> ?
      //typedef Dune::PDELab::ISTLBackend_SEQ_LOOP_Jac LS; // nö / ?

      //typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS; // nö / nö
      //typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILUn LS; // nö / ?
      //typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS; // nö / ?
      //typedef Dune::PDELab::ISTLBackend_SEQ_CG_Jac LS; // nö / ?
      //typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<IGO> LS; // nö / ?
      //typedef Dune::PDELab::ISTLBackend_SEQ_MINRES_SSOR LS; // nö / nö
      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_Jac LS; // nö / nö
      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS; // nö / nö

      // !!! Watch them constructor parameters !!!
      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0 LS; //nö / ja
      //LS ls(5000,5);
      typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILUn LS; // hardly / ja
      LS ls(1,1.0,maxLinIt,5);

      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<IGO,true> LS; // nö / ja
      //LS ls(maxLinIt);
      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SOR<IGO,true> LS; // ? / ja
      //LS ls(maxLinIt);
      //typedef Dune::PDELab::ISTLBackend_SEQ_LS_AMG_SSOR<IGO> LS; // nö / ja
      //LS ls(maxLinIt);
      //typedef Dune::PDELab::ISTLBackend_SEQ_LS_AMG_SOR<IGO> LS; //nö / ja
      //LS ls(maxLinIt);

      //typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
      //LS ls(5); // verbose = 1
#endif
      // This is a huge, fully demangled type name and therefore looks annoying in the output;
      // but it is very important to know which linear solver was used!
      debug_info << "Using linear solver " << Tools::getTypeName(ls) << std::endl;
#if 1
      // ==============================================================================================

      // ========= Solver for nonlinear problem per stage ==========================================
      typedef Ax1Newton<ACME2CYL_OUTPUT,LOP_ELEC,IGO,LS,typename Traits::U> PDESOLVER;
      std::string newtonPrefix(params.getOutputPrefix());
      newtonPrefix += "octave/acme2_cyl";
      PDESOLVER pdesolver(acme2_cylOutput,lopElec,igo,ls,newtonPrefix.c_str(),newtonSolution);
      pdesolver.setReassembleThreshold(0.0);
      pdesolver.setVerbosityLevel(5); // 2
      pdesolver.setReduction(params.getReduction()); // 1e-10
      pdesolver.setAbsoluteLimit(params.getAbsLimit()); // 1e-10
      pdesolver.setMinLinearReduction(params.general.get("newtonMinLinReduction",1e-5)); // has no effect when using direct solver like SuperLU
      pdesolver.setFixedLinearReduction(params.general.get("newtonFixedLinReduction",false));
      pdesolver.setMaxIterations(params.general.get("newtonMaxIt",50));
      pdesolver.setForceIteration(params.general.get("newtonForceIteration", false));
      // This version doesn't throw an exception when line search does not converge. One such case where this situation
      // can occur is when the very first Newton iteration results in an increase of the defect, while the following
      // iterations have the normal behaviour of reducing the defect.
      pdesolver.setLineSearchStrategy(PDESOLVER::hackbuschReuskenAcceptBestNoThrow);
      //pdesolver.setLineSearchStrategy(PDESOLVER::hackbuschReuskenAcceptBest);
      if(params.general.get("disableLineSearch",false)) // Optionally disable line search
        pdesolver.setLineSearchStrategy(PDESOLVER::noLineSearch);
      pdesolver.setLineSearchMaxIterations(10); // 10 (50)
      pdesolver.setPrintMatrix(params.doPrintMatrix());
      pdesolver.setPrintRhs(params.doPrintRhs());
      pdesolver.setRowPreconditioner(params.useRowNormPreconditioner());
      // Reordering is now done by a custom PDELab ordering, manual reordering is disabled!
      pdesolver.setReorderMatrix(false);
      pdesolver.setFullyImplicit(params.boundary.get("fullyImplicitMembraneFlux",false));
//      // Setup DOF permutations (necessary when reordering matrix in Ax1Newton)
//      if(params.doReorderMatrix())
//      {
//        pdesolver.setPermutation(perm);
//        pdesolver.setInversePermutation(inv_perm);
//      }
      pdesolver.setPrintResidual(params.doPrintResidual());
      // ==============================================================================================


      // ========== time-stepper ======================================================================
      //Dune::PDELab::Alexander2Parameter<Real> timeStepper;
      //Dune::PDELab::Alexander3Parameter<Real> timeStepper;
      Dune::PDELab::ImplicitEulerParameter<Real> timeStepper;
      // theta=0.5 => Crank-Nicholson
      //Dune::PDELab::OneStepThetaParameter<Real> timeStepper(0.5);
      //Dune::PDELab::ExplicitEulerParameter<Real> timeStepper;
      //Dune::PDELab::RK4Parameter<Real> timeStepper;
      //Dune::PDELab::HeunParameter<Real> timeStepper;
      typedef Dune::PDELab::OneStepMethod<Real,IGO,PDESOLVER,typename Traits::U,typename Traits::U> SOLVER;
      SOLVER solver(timeStepper,igo,pdesolver);
      solver.setVerbosityLevel(0); // 2

      // In case we need access to the pdesolver later on
      const PDESOLVER& pdesolverRef = solver.getPDESolver();
      // ==============================================================================================

      {
        igo.setMethod(timeStepper);
        typename IGO::Traits::Jacobian jac(igo);
        debug_info << " +++++++++++++++++++++++ Pattern statistics ++++++++++++++++++++++++++ " << std::endl;
        debug_info << jac.patternStatistics() << std::endl;
        debug_info << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
      }

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

        // Enforcing boundary conditions is removed; use the loaded values and leave them untouched!
        /*
        // Set Dirichlet boundary values from the saved data just loaded
        Real xReference = params.xMin() + 0.5 * (params.xMax() - params.xMin());
        Real yReference = params.yMin() + 0.5 * (params.yMax() - params.yMin());
        std::vector<Real> boundaryValues = physics.getBoundaryValues(dgfPot, xReference, yReference);
        initialGFPot.setPotValues(boundaryValues);
        debug_jochen << "Boundary values: ";
        Output::printVector(boundaryValues);

        // Enforce boundary conditions
        typename Traits::U uold_save = uold;
        Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,uold,initialElec,elecSubProblem);
        Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,uold,initialPot,membSubProblem);
        Dune::PDELab::copy_nonconstrained_dofs(cc,uold_save,uold);
        typename Traits::U unew_save = unew;
        Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,unew,initialElec,elecSubProblem);
        Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,unew,initialPot,membSubProblem);
        Dune::PDELab::copy_nonconstrained_dofs(cc,unew_save,unew);
        */

        // Keep loaded time and dt when continuing a loaded simulation
        if(not params.doContinueSimulation())
        {
          time = 0.0;
          dt = dtstart;
        } else {
          if(tend > 0 && tend < time)
            DUNE_THROW(Dune::Exception, "Start time (" << time << ") > end time (" << tend << ")!");
        }
        if(params.doEquilibration() && tEquilibrium > time)
        {
          dt = physics.getParams().dtEquilibrium();
          debug_info << "== Starting with initial dt = " << dt << " until tEquilibrium " << tEquilibrium
              << " ==" << std::endl;
        }
      }
      // ==============================================================================================

      // ========== Initial output ====================================================================
      std::vector<Real> debyeLength(2,0.0);
      physics.getDebyeLength(debyeLength);

      debug_info << std::endl;
      debug_info << "@@@@@@@ [Intracellular] Debye length / length scale " << debyeLength[0] / physics.getLengthScale()
                 << std::endl;
      debug_info << "@@@@@@@ [Extracellular] Debye length / length scale " << debyeLength[1] / physics.getLengthScale()
                 << std::endl;
      if(physics.getParams().dYMin() >= (debyeLength[0]/physics.getLengthScale())
          || physics.getParams().dYMin() >= (debyeLength[1]/physics.getLengthScale()))
      {
        debug_warn << "WARNING: Grid does not resolve Debye length [dy_min = " << physics.getParams().dYMin() << "]!" << std::endl;
      }
      debug_verb << "" << std::endl;

      debug_jochen << "Diffusion length for this initial dt: " << 2*std::sqrt(physics.getDiffCoeff(Na, 0) * dt) << std::endl;

      // Test cylinder coordinates stuff
      /*
      Acme2CylGeometryCheck<Grid> acme2CylChecks(params);

      // Integrate over elements (--> volumes)
      typename Traits::DGF_CON::Traits::RangeType sum(0.0);
      Acme2CylGeometryTools::integrateGridFunctionOverCylinder(dgfOldCon,sum,2);
      debug_jochen << "Integral over concentrations: " << std::endl;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        debug_jochen << "  " << ION_NAMES[j] << ": " << sum[j] << std::endl;
      }

      debug_jochen << std::endl;
      debug_jochen << std::endl;

      // Integrate over intersections (--> areas)
      typename Traits::GF_MEMB_FLUX::Traits::RangeType sum2(0.0);

      Acme2CylGeometryTools::integrateBoundaryGridFunctionOverCylinder(gfMembFlux,sum2,2);
      debug_jochen << "Integral over memb fluxes: " << std::endl;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        debug_jochen << "  " << ION_NAMES[j] << ": " << sum2[j] << std::endl;
      }
      */
      // ==============================================================================================

      //Output::printRawCoefficientVector(unew, "unew after load", 30);

      // ========== Run simulation ====================================================================
      Acme2CylSimulation<Traits,ACME2CYL_OUTPUT>::run(time, dt, dtstart, tend, tEquilibrium,
          physics, gv, membGV,
          solver,
          parametersPot,
          parametersCon,
          multigfs,
          uold, unew,
          cc,ccWithoutOverlap,
          gfMembFlux, gfMoriFlux, dgfCon, dgfPot, dgfGradPot,
          acme2_cylOutput, solutionVectors,
          initialElec, elecSubProblem);
      // ==============================================================================================

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
#endif
    }

  private:


    Grid& grid;

    GV& gv;
    PHYSICS& physics;

    SubGV& elecGV;
    SubGV& membGV;
};

#endif /* DUNE_AX1_ACME2CYL_SETUP_HH */
