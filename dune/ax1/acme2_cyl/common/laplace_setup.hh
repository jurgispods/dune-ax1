#ifndef DUNE_AX1_LAPLACE_SETUP_HH
#define DUNE_AX1_LAPLACE_SETUP_HH

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
#include <dune/ax1/common/ionicstrengthgridfunction.hh>
#include <dune/ax1/common/membranefluxgridfunction.hh>
#include <dune/ax1/common/poisson_boltzmann_concentrationgridfunction.hh>
#include <dune/ax1/common/poisson_boltzmann_rhs_gridfunction.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/vectordiscretegridfunctiongradient.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_boundary.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_geometrytools.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_initial.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_solutionvectors.hh>
#include <dune/ax1/acme2_cyl/common/laplace_simulation.hh>
#include <dune/ax1/acme2_cyl/common/laplace_output.hh>
#include <dune/ax1/acme2_cyl/configurations/laplace/laplace_parameters.hh>
#include <dune/ax1/acme2_cyl/operator/convectiondiffusionfem.hh>

#define USE_NEWTON 0

template<class Grid, class GV, class PHYSICS, class SubGV>
class LaplaceSetup
{
  public:
    
    //! Traits class providing function spaces and related data types for the acme2_cyl use case
    template<int NUMBER_OF_SPECIES>
    struct LaplaceTraits
    {
      typedef LaplaceSetup<Grid,GV,PHYSICS,SubGV> SETUP;

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
      typedef Dune::PDELab::ISTLVectorBackend<> VBE;

      typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,GV::dimension> FEM;

      // Laplace/Poisson GFS
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CONSTRAINTS,VBE> GFS;

      // Multidomain GFS
    typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<GridType,VBE,
        Dune::PDELab::LexicographicOrderingTag,GFS> MultiGFS;

      // Extract SubGFS for use in gridfunctions; use these after creating MultiGFS!
      typedef typename Dune::PDELab::GridFunctionSubSpace<MultiGFS,Dune::PDELab::TypeTree::TreePath<0> > GFS_POT_SUB;

      // Extract solution vector type
      typedef typename Dune::PDELab::BackendVectorSelector<MultiGFS,Real>::Type U;

      // Various gridfunctions
      typedef Dune::PDELab::DiscreteGridFunction<GFS_POT_SUB,U> DGF_POT;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_POT_SUB,U> DGF_POT_GRAD;

      // Poisson parameter class with template parameter 'isLaplace" set to true
      typedef LaplaceParameters<GV,Real,PHYSICS> PARAMETERS_POT;

      // Grid functions for initial values (subdomains)
      typedef typename PHYSICS::Traits::INITIAL_POT INITIAL_GF_POT; //obsolete
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PARAMETERS_POT> INITIAL_POT;
      typedef INITIAL_POT INITIAL_ELEC;

      // Boundary values (obsolete)
      typedef typename PHYSICS::Traits::POISSON_BOUNDARY BOUNDARY_POT;
    };

    // Choose domain and range field type
    typedef GV GridView;
    typedef typename GV::Grid GridType;
    typedef typename GV::Grid::ctype Coord;
    typedef double Real;

    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

    typedef LaplaceTraits<NUMBER_OF_SPECIES> Traits;

    static const bool writeIfNotConverged = true;

    LaplaceSetup(Grid& grid_, GV& gv_, PHYSICS& physics_,
        SubGV& elecGV_, SubGV& membGV_)
    : grid(grid_),
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
      const int degreePot = 1;

			// Finite element map for CG
      typename Traits::FEM fem;

      // =================== GFS SETUP ================================================================
#if USE_PARALLEL == 1 && USE_OVERLAP == 0
      typename Traits::CONSTRAINTS constraints(gv);
#else
      typename Traits::CONSTRAINTS constraints;
#endif

      // GFS (only one on the whole grid)
      typename Traits::GFS gfs(gv,fem,constraints);
#if USE_PARALLEL == 1 && USE_OVERLAP == 0
      if(gv.grid().ghostSize(0) > 0)
      {
        constraints.compute_ghosts(gfs);
      }
#endif
      // multigfs creation for LexicographicOrderingTag
      typename Traits::MultiGFS multigfs(grid,gfs);
      debug_verb << "Instantiated MultiGFS." << std::endl;

      debug_verb << "Setting up concentration Sub GFS..." << std::endl;
      // SubGridFunctionSpaces for use in gridfunctions
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

      // Total 'real' number of DOFs
      int nTotalDofsWithoutOverlap = params.nNodes() * 1;
      debug_info << "Total #DOFs of the sequential problem (i.e., without overlap): " << nTotalDofsWithoutOverlap << std::endl;
      debug_info << "===============================" << std::endl;


      // Now let's put 'em all into one structure!
      Acme2CylSolutionVectors<Traits> solutionVectors(uold, unew);
      // ==============================================================================================

      // ====== Create grid functions =================================================================
      typename Traits::DGF_POT dgfPot(gfsPotSub, unew);
      typename Traits::DGF_POT dgfOldPot(gfsPotSub, uold);
      typename Traits::DGF_POT_GRAD dgfGradPot(gfsPotSub, unew);
      // ==============================================================================================


      // ========== Define Parameter classes containing the model problem =============================
      // In the default acme2_cyl(_par) setup, the Poisson constant is tailored to using a charge density
      // rhs in the Poisson equation; in volume conductor theory, however, the rhs is a _current_ density
      // and we need to correct for this, as well as for the units of the loaded boundary (membrane) flux
      Real poissonConstant = physics.getPoissonConstant();
      poissonConstant *= (con_eps0 * physics.getLengthScale() / physics.getTimeScale());
      typename Traits::PARAMETERS_POT parametersPot(physics,poissonConstant);
      parametersPot.setTime(time); // Necessary for loading initial boundary values from file
      // ==============================================================================================


      // ==============================================================================================
      // Define and instantiate output class
      typedef LaplaceOutput<Traits> LAPLACE_OUTPUT;
      LAPLACE_OUTPUT laplaceOutput(gv,elecGV,membGV,multigfs,gfsPotSub,solutionVectors,
          physics,2*degreePot);
      // ==============================================================================================

      // ============== BOUNDARY CONDITION TYPES ======================================================
      typedef BCTypePot<typename Traits::PARAMETERS_POT,PHYSICS> BCType_POT;
      BCType_POT bctypePot(gv,parametersPot,physics);

      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter
          <typename Traits::PARAMETERS_POT> BCTypePot;
      BCTypePot bctype(gv, parametersPot);
      // ==============================================================================================

      // ============== BOUNDARY / INITIAL VALUES======================================================
      typename Traits::INITIAL_POT initialPot(gv,parametersPot);
      // ==============================================================================================

      // ============== Define local operators ========================================================
      // ### Standard FEM (CG) Laplace local operator
      const int intorderadd = 0;
      // Template parameter 'useMembraneContributions' set to false for this setup!
      typedef Dune::PDELab::ConvectionDiffusionFEM<typename Traits::PARAMETERS_POT,typename Traits::FEM,
          true> LOP;
      LOP lop(parametersPot, intorderadd);
      // ==============================================================================================


      // =========== dune-multidomain setup stuff =====================================================
      typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> EC;
      EC conditionElec(0); //only on subdomain 0
      EC conditionMemb(1); //only on subdomain 1

      typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<typename Traits::MultiGFS,typename Traits::MultiGFS,
          LOP,EC,typename Traits::GFS> ElecSubProblem;
      ElecSubProblem elecSubProblem(lop,conditionElec);

      typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<typename Traits::MultiGFS,typename Traits::MultiGFS,
          LOP,EC,typename Traits::GFS> MembSubProblem;
      MembSubProblem membSubProblem(lop,conditionMemb);
      // ==============================================================================================

      // =================== BOUNDARY CONDITIONS ======================================================
      typedef typename Traits::MultiGFS::template ConstraintsContainer<Real>::Type CC;
      CC cc;

      auto md_constraints = Dune::PDELab::MultiDomain::template constraints<Real>(multigfs,
            Dune::PDELab::MultiDomain::constrainSubProblem(elecSubProblem, bctypePot),
            Dune::PDELab::MultiDomain::constrainSubProblem(membSubProblem, bctypePot));

      md_constraints.assemble(cc);
      debug_info << multigfs.size() << " DOF, " << cc.size() << " restricted" << std::endl;

      // Create one additional contraints container, which is only used for the very special case of
      // interpolating time-dependent Dirichlet boundary values in the overlapping parallel case.
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
      typedef Dune::PDELab::ISTLMatrixBackend MBE;
      //typedef Dune::PDELab::GridOperator<GFS_POT,GFS_POT,STATIONARY_LOP_POT,MBE,Real,Real,Real,CC_POT,CC_POT>
      //  STATIONARY_GO_POT;
      //STATIONARY_GO_POT stationaryGoPot(gfsPot,ccPot,gfsPot,ccPot,stationaryLopPot);

      typedef Dune::PDELab::MultiDomain::GridOperator<typename Traits::MultiGFS,typename Traits::MultiGFS,
        MBE,Real,Real,Real,CC,CC,ElecSubProblem,MembSubProblem> GO;
      GO go(multigfs,multigfs,cc,cc,elecSubProblem,membSubProblem);

      // Here, initial values coming from the both subdomain initial GFs are written to uold
      Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,uold,initialPot,elecSubProblem);
      Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,uold,initialPot,membSubProblem);
      unew = uold;
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
      // TODO Test CG solvers for Laplace problem
      //typedef Dune::PDELab::ISTLBackend_OVLP_CG_SSORk<typename Traits::MultiGFS,CC> LS;
      //LS ls(multigfs,cc,maxLinIt,5,5);
      //typedef Dune::PDELab::ISTLBackend_OVLP_CG_ILU0<typename Traits::MultiGFS,CC> LS;
      //LS ls(multigfs,cc,maxLinIt,5);
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILU0<typename Traits::MultiGFS,CC> LS;
      LS ls(multigfs,cc,maxLinIt,5);
      //typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILUn<typename Traits::MultiGFS,CC> LS;
      //LS ls(multigfs,cc,1,maxLinIt,5);
      //typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<typename Traits::MultiGFS,CC> LS;
      //LS ls(multigfs,cc,maxLinIt,5);
      //typedef Dune::PDELab::ISTLBackend_OVLP_GMRES_ILU0<typename Traits::MultiGFS,CC> LS;
      //LS ls(multigfs,cc,maxLinIt,5,40,true);
      //typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<IGO> LS;
      //LS ls(multigfs,maxLinIt,5);
      // This is the preferred preconditioner/linear solver combination for 'hard' parallel problems
//      typedef Dune::PDELab::ISTLBackend_GMRES_AMG_ILU0<GO> LS;
//      LS ls(multigfs,maxLinIt,5);
//
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
      typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_ILUn<GO> LS;
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

      // ========= Solver for linear problem ==========================================================
#if USE_NEWTON==1
      typedef Dune::PDELab::Newton<GO,LS,typename Traits::U> PDESOLVER;
      PDESOLVER pdesolver(go,ls);
      pdesolver.setReassembleThreshold(0.0);
      pdesolver.setVerbosityLevel(5); // 2
      pdesolver.setReduction(params.getReduction()); // 1e-10
      pdesolver.setAbsoluteLimit(params.getAbsLimit()); // 1e-10
      pdesolver.setMinLinearReduction(params.general.get("newtonMinLinReduction",1e-5)); // has no effect when using direct solver like SuperLU
      pdesolver.setFixedLinearReduction(params.general.get("newtonFixedLinReduction",false));
      pdesolver.setMaxIterations(params.general.get("newtonMaxIt",50));
      pdesolver.setForceIteration(params.general.get("newtonForceIteration", false));
      pdesolver.setLineSearchStrategy(PDESOLVER::noLineSearch);
      pdesolver.setLineSearchMaxIterations(10); // 10 (50)
#else
      typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,typename Traits::U> PDESOLVER;
      PDESOLVER pdesolver(go,ls,params.getReduction());
      pdesolver.setHangingNodeModifications(false);
#endif
      // ==============================================================================================

      // ========== Load saved state ======================================================================
      if(physics.getParams().doLoadState())
      {
        //solutionVectors.printHostGridVectors();
        std::map<std::string,std::string> loadFileNames = laplaceOutput.loadState(time,dt);
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

      //Output::printRawCoefficientVector(unew, "unew");

      // ========== Run simulation ====================================================================
      LaplaceSimulation<Traits,LAPLACE_OUTPUT>::run(time, dt, dtstart, tend, tEquilibrium,
          physics, gv, membGV,
          pdesolver,
          multigfs,
          uold, unew,
          cc,ccWithoutOverlap,
          dgfPot, dgfGradPot,
          laplaceOutput, solutionVectors,
          initialPot, elecSubProblem);
      // ==============================================================================================

      // ======= Save simulation state ================================================================
      if(physics.getParams().doSaveState())
      {
        std::string saveFilename = physics.getParams().getSaveFilename();
        laplaceOutput.saveState(time,dt,saveFilename);
        debug_info << "=============================================================" << std::endl;
        debug_info << "Saved simulation state to file " << saveFilename << std::endl;
        debug_info << "=============================================================" << std::endl;
      }
      // ==============================================================================================

      // Print #DOFs at simulation end
      debug_info << "U.N(): " << unew.N() << std::endl;
      debug_info << "U.flatsize(): " << unew.flatsize() << std::endl;
      debug_info << "GFS size: " << gfsPotSub.size() << std::endl;

      debug_info << "Total Acme2CylOutput time: " << laplaceOutput.getTotalOutputTime() << "s" << std::endl;
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
