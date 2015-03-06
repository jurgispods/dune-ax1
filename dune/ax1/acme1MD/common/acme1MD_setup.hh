#ifndef DUNE_AX1_ACME1MD_SETUP_HH
#define DUNE_AX1_ACME1MD_SETUP_HH

#include <dune/common/array.hh>

#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/finiteelementmap/p1fem.hh>
#include <dune/pdelab/finiteelementmap/pk1dbasis.hh>
#include <dune/pdelab/finiteelementmap/opbfem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/function/selectcomponent.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/backend/seqistlsolverbackend.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

#include <dune/ax1/common/ax1_newton.hh>
#include <dune/ax1/common/ax1_linearproblem.hh>
#include <dune/ax1/common/ax1_discretegridfunction.hh>
#include <dune/ax1/common/ax1_simulationdata.hh>
#include <dune/ax1/common/ax1_simulationstate.hh>
#include <dune/ax1/common/chargedensitygridfunction.hh>
#include <dune/ax1/common/ionicstrengthgridfunction.hh>
#include <dune/ax1/common/membranefluxgridfunction.hh>
#include <dune/ax1/common/poisson_boltzmann_concentrationgridfunction.hh>
#include <dune/ax1/common/poisson_boltzmann_rhs_gridfunction.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/vectordiscretegridfunctiongradient.hh>
#include <dune/ax1/acme1MD/common/acme1MD_boundary.hh>
#include <dune/ax1/acme1MD/common/acme1MD_initial.hh>
#include <dune/ax1/acme1MD/common/acme1MD_output.hh>
#include <dune/ax1/acme1MD/common/acme1MD_solutionvectors.hh>
#include <dune/ax1/acme1MD/common/convectiondiffusionfem.hh>
#include <dune/ax1/acme1MD/common/convectiondiffusiondg.hh>
#include <dune/ax1/acme1MD/common/nernst_planck_parameters.hh>
#include <dune/ax1/acme1MD/common/poisson_parameters.hh>
#include <dune/ax1/acme1MD/common/poisson_boltzmann_parameters.hh>
#include <dune/ax1/acme1MD/common/poisson_boltzmann_operator.hh>
#include <dune/ax1/acme1MD/common/proportionalflowcoupling.hh>
#include <dune/ax1/acme1MD/fully-implicit/acme1MD_experimental_operator_fully_implicit.hh>
#include <dune/ax1/acme1MD/fully-implicit/acme1MD_fully_coupled.hh>
#include <dune/ax1/acme1MD/operator-split/acme1MD_operator_split.hh>
#include <dune/ax1/acme1MD/operator-split/acme1MD_iteration.hh>
#include <dune/ax1/acme1MD/operator-split/acme1MD_toperator.hh>
#include <dune/ax1/acme1MD/operator-split/nernst_planck_power_operator.hh>

template<class Grid, class GV, class PHYSICS, class SubGV>
class Acme1MDSetup
{
  public:
    
    //! Traits class providing function spaces and related data types for the acme1MD use case
    template<int NUMBER_OF_SPECIES>
    struct Acme1MDTraits
    {
      typedef Acme1MDSetup<Grid,GV,PHYSICS,SubGV> SETUP;

      typedef PHYSICS Physics;

      typedef GV GridView;
      typedef Grid GridType;
      typedef typename GV::Grid::ctype Coord;
      typedef double Real;

      typedef SubGV SubGridView;
      typedef typename SubGV::Grid SubGridType;

      typedef Dune::PDELab::ConformingDirichletConstraints CONSTRAINTS;
      typedef Dune::PDELab::ISTLVectorBackend<1> VBE;

      typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM_CON;
      typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM_POT;

      // Nernst-Planck GFS (on electrolyte subdomain)
      typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper PowerGFSOrdering;
      typedef Dune::PDELab::GridFunctionSpace<SubGV,FEM_CON,CONSTRAINTS,VBE> GFS_SINGLE_CON;
      typedef Dune::PDELab::PowerGridFunctionSpace<
        GFS_SINGLE_CON, NUMBER_OF_SPECIES, PowerGFSOrdering> GFS_CON;

      // Poisson GFS (on electrolyte or membrane subdomain)
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM_POT,CONSTRAINTS,VBE> GFS_POT;

      // Composite GFS for electrolyte
      typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper CompositeGFSOrdering;

      // Multidomain GFS
      typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<GridType,VBE,GFS_CON,GFS_POT> MultiGFS;

      // Extract solution vector types
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CON,Real>::Type U_CON;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_POT,Real>::Type U_POT;
      typedef typename Dune::PDELab::BackendVectorSelector<MultiGFS,Real>::Type U;

      typedef Acme1MDOutput<GV,SubGV,GFS_CON,GFS_POT,U_CON,U_POT,PHYSICS> ACME1MD_OUTPUT;
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
    typedef Acme1MDTraits<NUMBER_OF_SPECIES> Traits;

    static const bool useRowPreconditioner = true;
    static const bool writeIfNotConverged = true;
    static const bool resetTime = true;

    Acme1MDSetup(Grid& grid_, GV& gv_, PHYSICS& physics_,
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

      Real tEquilibrium = physics.getParams().tEquilibrium();
      if(tEquilibrium > time)
      {
        dt = physics.getParams().dtEquilibrium();
        debug_info << "== Starting with initial dt = " << dt << " until tEquilibirum " << tEquilibrium
            << " ==" << std::endl;
      }

      const int dim = GV::dimension;
      const int degreeCon=1;
      const int degreePot=1;

			// Finite element map for CG
      typename Traits::FEM_CON femCon(degreeCon);
      typename Traits::FEM_POT femPot(degreePot);

      // =================== GFS SETUP ================================================================
      // Single-component GFS for one ion species
      typename Traits::GFS_SINGLE_CON gfsSingleCon(elecGV,femCon);
      // Power grid function space for all ion species
      typename Traits::GFS_CON gfsCon(gfsSingleCon);

      // Full GFS for potential on electrolyte and membrane subdomain
      typename Traits::GFS_POT gfsPot(gv,femPot);

      // Multidomain GFS
      typename Traits::MultiGFS multigfs(grid,gfsCon,gfsPot);

      // Types for different grid functions
      typedef Dune::PDELab::VectorDiscreteGridFunction <typename Traits::GFS_CON,typename Traits::U_CON> DGF_CON_SUB;
      typedef Dune::PDELab::VectorDiscreteGridFunctionGradient <typename Traits::GFS_CON,typename Traits::U_CON> DGF_CON_GRAD_SUB;

      typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_CON_SUB, PHYSICS> DGF_CON;
      typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_CON_GRAD_SUB, PHYSICS> DGF_CON_GRAD;

      typedef Dune::PDELab::DiscreteGridFunction<typename Traits::GFS_POT,typename Traits::U_POT> DGF_POT;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<typename Traits::GFS_POT,typename Traits::U_POT> DGF_POT_GRAD;

      typedef ChargeDensityGridFunction<DGF_CON, PHYSICS> GF_CD;
      typedef IonicStrengthGridFunction<DGF_CON, PHYSICS> GF_IS;

      typedef MembraneFluxGridFunction<DGF_CON,DGF_POT,PHYSICS> GF_MEMB_FLUX;

      // Grid functions for initial values (subdomains)
      typedef typename PHYSICS::Traits::INITIAL_CON INITIAL_GF_CON;
      typedef InitialCon<GV,Real,NUMBER_OF_SPECIES,PHYSICS,INITIAL_GF_CON> INITIAL_CON;
      typedef typename PHYSICS::Traits::INITIAL_POT INITIAL_GF_POT;
      typedef InitialPot<GV,Real,PHYSICS,GF_CD,INITIAL_GF_POT> INITIAL_POT;

      // Helper grid function to calculate equilibirum concentrations from stationary Poisson-Boltzmann potential
      //typedef PoissonBoltzmannRHSGridFunction<DGF_POT, INITIAL_CON, PHYSICS> GF_PB_RHS;
      //typedef PoissonBoltzmannConcentrationGridFunction<INITIAL_CON, DGF_POT, SubGV, PHYSICS> GF_PB_CON;
      // ==============================================================================================


      // =========== Define solution vectors ==========================================================
      // Coefficient vectors on subdomains
      typename Traits::U_CON uoldCon(gfsCon,0.0);
      typename Traits::U_CON unewCon(gfsCon,0.0);
      debug_jochen << "U_CON.N(): " << unewCon.N() << std::endl;
      debug_jochen << "U_CON.flatsize(): " << unewCon.flatsize() << std::endl;

      typename Traits::U_POT uoldPot(gfsPot,0.0);
      typename Traits::U_POT unewPot(gfsPot,0.0);
      debug_jochen << "U_POT.N(): " << unewPot.N() << std::endl;
      debug_jochen << "U_POT.flatsize(): " << unewPot.flatsize() << std::endl;

      typename Traits::U uold(multigfs,0.0);
      typename Traits::U unew = uold;
      debug_jochen << "U.N(): " << unew.N() << std::endl;
      debug_jochen << "U.flatsize(): " << unew.flatsize() << std::endl;

      // Now let's put 'em all into one structure! //TODO alle rin!
      Acme1MDSolutionVectors<Traits> solutionVectors(uold, unew, uoldCon, unewCon, uoldPot, unewPot);
      // ==============================================================================================

      // ====== Create grid functions and interpolate initial values onto coefficient vectors =========
      // Initial concentrations
      INITIAL_GF_CON initialGFCon(gv,physics.getParams());

      // Generic wrapper class for initial grid function
      INITIAL_CON initialCon(gv,physics,initialGFCon);
      initialCon.setTime(time);

      DGF_CON_SUB dgfConElec(gfsCon,unewCon);
      DGF_CON_SUB dgfOldConElec(gfsCon,uoldCon);
      DGF_CON_GRAD_SUB dgfGradConElec(gfsCon, unewCon);

      DGF_CON dgfCon(gv, dgfConElec, physics);
      DGF_CON dgfOldCon(gv, dgfOldConElec, physics);
      DGF_CON_GRAD dgfGradCon(gv, dgfGradConElec, physics);

      GF_CD gfChargeDensity(dgfCon, dgfOldCon, physics);
      GF_IS gfIonicStrength(dgfCon, physics);

      // Initial potential
      INITIAL_GF_POT initialGFPot(gv,physics.getParams());
      INITIAL_POT initialPot(gv,physics,gfChargeDensity,initialGFPot,2*degreePot);
      initialPot.setTime(time);

      DGF_POT dgfPot(gfsPot, unewPot);
      DGF_POT dgfOldPot(gfsPot, uoldPot);
      DGF_POT_GRAD dgfGradPot(gfsPot, unewPot);

      // Flag 'true' for updating channel states in each time step
      GF_MEMB_FLUX gfMembFlux(dgfOldCon, dgfOldPot, physics, true, tEquilibrium);

      // Composite initial grid function for the electrolyte subdomain
      typedef Dune::PDELab::CompositeGridFunction<INITIAL_CON,INITIAL_POT> INITIAL_ELEC;
      INITIAL_ELEC initialElec(initialCon,initialPot);

      // Define and instantiate output class
      typename Traits::ACME1MD_OUTPUT acme1MDOutput(gv,elecGV,membGV,gfsCon,gfsPot,unewCon,unewPot,
          gfMembFlux,physics,2*degreePot,2*degreeCon);
      // ==============================================================================================


      // ========== Define Parameter classes containing the model problem =============================
      typedef typename PHYSICS::Traits::NERNST_PLANCK_BOUNDARY BOUNDARY_CON;
      BOUNDARY_CON boundaryCon(initialGFCon);
      typedef NernstPlanckParameters<GV,Real,PHYSICS,BOUNDARY_CON,GF_MEMB_FLUX> PARAMETERS_CON;
      PARAMETERS_CON parametersCon(gv,physics,boundaryCon,gfMembFlux,tEquilibrium);

      typedef typename PHYSICS::Traits::POISSON_BOUNDARY BOUNDARY_POT;
      BOUNDARY_POT boundaryPot;
      typedef PoissonParameters<GV,Real,PHYSICS,BOUNDARY_POT> PARAMETERS_POT;
      PARAMETERS_POT parametersPot(physics,boundaryPot);
      // ==============================================================================================

      // ============== BOUNDARY CONDITION TYPES ======================================================
      typedef BCTypeSingleCon<PARAMETERS_CON,PHYSICS> BCType_SINGLE_CON;
      typedef Dune::PDELab::PowerConstraintsParameters<BCType_SINGLE_CON,NUMBER_OF_SPECIES> BCType_CON;
      typedef BCTypePot<PARAMETERS_POT,PHYSICS> BCType_POT;

      // Create NUMBER_OF_SPECIES boundary condition type classes
      BCType_CON bctypeCon;
      for(int i=0; i<NUMBER_OF_SPECIES; ++i)
      {
        BCType_SINGLE_CON* bctypeSingleCon = new BCType_SINGLE_CON(gv,parametersCon,physics,i);
        bctypeCon.setChild(i,*bctypeSingleCon);
      }
      BCType_POT bctypePot(gv,parametersPot,physics);

      typedef Dune::PDELab::CompositeConstraintsParameters<BCType_CON, BCType_POT> BCType_ELEC;
      BCType_ELEC bctypeElec(bctypeCon, bctypePot);
      // ==============================================================================================


      // ============== Define local operators ========================================================
      // ##### Stationary operator #####
      //typedef PoissonBoltzmannLocalOperator<PHYSICS,PARAMETERS_POT,INITIAL_CON> STATIONARY_LOP_POT;
      //STATIONARY_LOP_POT stationaryLopPot(physics, parametersPot, initialCon);

//      GF_PB_RHS gfPoissonBoltzmannRHS(dgfPot, initialCon, physics);
//      typedef PoissonBoltzmannParameters<GV,Real,PHYSICS,GF_PB_RHS,BOUNDARY_POT> PARAMETERS_PB;
//      PARAMETERS_PB parametersPB(physics,gfPoissonBoltzmannRHS,boundaryPot);
//
//      typedef Dune::PDELab::ConvectionDiffusionFEM<PARAMETERS_PB,FEM_POT> STATIONARY_LOP_POT;
//      STATIONARY_LOP_POT stationaryLopPot(parametersPB);

      // ##### Spatial part #####
      // ### Standard FEM (CG) fully-coupled Poisson-Nernst-Planck local operator on electrolyte subdomain
      typedef Acme1MDExperimentalOperatorFullyImplicit<PARAMETERS_CON,PARAMETERS_POT,typename Traits::FEM_CON,
          typename Traits::FEM_POT> LOP_ELEC;
      LOP_ELEC lopElec(parametersCon, parametersPot);

      // ### Standard FEM (CG) Poisson local operator on membrane subdomain
      typedef Dune::PDELab::ConvectionDiffusionFEM<PARAMETERS_POT,typename Traits::FEM_POT> LOP_MEMB;
      LOP_MEMB lopMemb(parametersPot);

      // ##### Temporal part #####
      typedef NernstPlanckTimeLocalOperator<PHYSICS> TLOP_ELEC;
      TLOP_ELEC tlop(physics);
      // ==============================================================================================


      // =========== dune-multidomain setup stuff =====================================================
      typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> EC;
      typedef Dune::PDELab::MultiDomain::SubDomainSupersetCondition<Grid> SC;
      EC conditionElec(0); //only on subdomain 0
      EC conditionMemb(1); //only on subdomain 1
      SC conditionAll;    //on all subdomains

      // Empty constraints for subproblems
      typename Traits::CONSTRAINTS constraints;

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
      //typedef Traits::MultiGFS::ConstraintsContainer<Real>::Type CC;
      typedef typename Traits::MultiGFS::template ConstraintsContainer<Real>::Type CC;
      CC cc;

      auto md_constraints = Dune::PDELab::MultiDomain::template constraints<Real>(multigfs,
            Dune::PDELab::MultiDomain::constrainSubProblem(elecSubProblem, bctypeElec),
            Dune::PDELab::MultiDomain::constrainSubProblem(membSubProblem, bctypePot));

      md_constraints.assemble(cc);
      std::cout << multigfs.size() << " DOF, " << cc.size() << " restricted" << std::endl;

      // ###### Dirichlet values ########
      typedef DirichletValuesSingleCon<PARAMETERS_CON> DirichletValuesSingleCon;

      // Create NUMBER_OF_SPECIES Dirichlet value classes
      typedef Dune::PDELab::PowerGridFunction<DirichletValuesSingleCon,NUMBER_OF_SPECIES> DirichletValuesCon;
      DirichletValuesCon dirichletValuesCon;

      for(int i=0; i<NUMBER_OF_SPECIES; ++i)
      {
        DirichletValuesSingleCon* dirichletValuesSingleCon
          = new DirichletValuesSingleCon(gv,parametersCon,i);
        dirichletValuesCon.setChild(i,*dirichletValuesSingleCon);
      }
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PARAMETERS_POT> DirichletValuesPot;
      DirichletValuesPot dirichletValuesPot(gv,parametersPot);

      // Is this necessary? Initial conditions should fulfill the boundary conditions anyways,
      // otherwise we have an ill-posed problem!
      //Dune::PDELab::interpolate(dirichletValuesCon_Inside,subGfsCon_Inside,unewCon_Inside);
      //Dune::PDELab::copy_nonconstrained_dofs(ccCon_Inside,uoldCon_Inside,unewCon_Inside);
      //Dune::PDELab::interpolate(dirichletValuesPot,gfsPot,uPot);
      //Output::printSingleCoefficientVector(uPot, "uPot t=0");
      // ==============================================================================================

      // ============== Make grid operator ============================================================
      typedef typename Traits::VBE::MatrixBackend MBE;
      //typedef Dune::PDELab::GridOperator<GFS_POT,GFS_POT,STATIONARY_LOP_POT,MBE,Real,Real,Real,CC_POT,CC_POT>
      //  STATIONARY_GO_POT;
      //STATIONARY_GO_POT stationaryGoPot(gfsPot,ccPot,gfsPot,ccPot,stationaryLopPot);

      typedef Dune::PDELab::MultiDomain::GridOperator<typename Traits::MultiGFS,typename Traits::MultiGFS,
        MBE,Real,Real,Real,CC,CC,ElecSubProblem,MembSubProblem> GO0;
      typedef Dune::PDELab::MultiDomain::GridOperator<typename Traits::MultiGFS,typename Traits::MultiGFS,
        MBE,Real,Real,Real,CC,CC,TimeSubProblem> GO1;
      typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
      GO0 go0(multigfs,multigfs,cc,cc,elecSubProblem,membSubProblem);
      GO1 go1(multigfs,multigfs,cc,cc,timeSubProblem);
      IGO igo(go0,go1);

      Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,uold,initialElec,elecSubProblem);
      Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,uold,initialPot,membSubProblem);
      unew = uold;
      // ==============================================================================================

      // ========== Select a linear solver backend ====================================================
      //typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS;
      //LS ls(5000,false);
      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
      //LS ls(5000,false);
      typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
      LS ls(0); // verbose = 1
      // ==============================================================================================

      // ========= Solver for (non-)linear problem per stage ==========================================
      // Stationary Poisson-Boltzmann solver
//      typedef Ax1Newton<Acme1MDOutput,STATIONARY_GO_POT,LS_POT,U_POT> STATIONARY_SOLVER_POT;
//      std::string stationaryOutputPrefix = "octave/";
//      stationaryOutputPrefix += stationaryLopPot.getName();
//      STATIONARY_SOLVER_POT stationarySolverPot(acme1MDOutput,stationaryGoPot,lsPot,
//          stationaryOutputPrefix.c_str());
//      stationarySolverPot.setReassembleThreshold(0.0);
//      stationarySolverPot.setVerbosityLevel(0); // (2)
//      stationarySolverPot.setAbsoluteLimit(1e-8);
//      stationarySolverPot.setReduction(1e-6);
//      stationarySolverPot.setMinLinearReduction(1e-4);
//      stationarySolverPot.setMaxIterations(25);
//      stationarySolverPot.setLineSearchMaxIterations(10);
//      stationarySolverPot.setPrintMatrix(false);
//      stationarySolverPot.setPrintRhs(false);

      typedef Ax1Newton<typename Traits::ACME1MD_OUTPUT,IGO,LS,typename Traits::U> PDESOLVER;
      PDESOLVER pdesolver(acme1MDOutput,igo,ls,"octave/acme1MD_mat");
      pdesolver.setReassembleThreshold(0.0);
      pdesolver.setVerbosityLevel(2); // 2
      pdesolver.setReduction(physics.getParams().getReduction()); // 1e-10
      pdesolver.setAbsoluteLimit(physics.getParams().getAbsLimit()); // 1e-10
      pdesolver.setMinLinearReduction(1e-3);
      pdesolver.setMaxIterations(30);
      //pdesolver.setLineSearchStrategy(PDESOLVER::noLineSearch);
      pdesolver.setLineSearchMaxIterations(50); // 10
      pdesolver.setPrintMatrix(false);
      pdesolver.setPrintRhs(false);
      pdesolver.setRowPreconditioner(useRowPreconditioner);
      // ==============================================================================================


      // ========== time-stepper ======================================================================
      //Dune::PDELab::Alexander3Parameter<Real> timeStepper;
      //Dune::PDELab::Alexander3Parameter<Real> timeStepper;
      Dune::PDELab::ImplicitEulerParameter<Real> timeStepper;
      //Dune::PDELab::ExplicitEulerParameter<Real> timeStepper;
      //Dune::PDELab::RK4Parameter<Real> timeStepper;
      //Dune::PDELab::HeunParameter<Real> timeStepper;
      typedef Dune::PDELab::OneStepMethod<Real,IGO,PDESOLVER,typename Traits::U,typename Traits::U> SOLVER;
      SOLVER solver(timeStepper,igo,pdesolver);
      solver.setVerbosityLevel(0); // 2
      
      // In case we need access to the pdesolver later on
      const PDESOLVER& pdesolverRef = solver.getPDESolver();
      // ==============================================================================================
      

      // ========== Load saved state ======================================================================
      if(physics.getParams().doLoadState())
      {
        //solutionVectors.printHostGridVectors();
        std::string loadFilename = physics.getParams().getLoadFilename();
        loadState(time,dt,solutionVectors,loadFilename);
        debug_info << "======================================================================================"
            << std::endl;
        debug_info << "Loaded simulation state from file " << loadFilename << std::endl;
        debug_info << "======================================================================================"
            << std::endl;

        if(resetTime)
        {
          time = 0.0;
        }
        dt = dtstart;
        if(tEquilibrium > time)
        {
          dt = physics.getParams().dtEquilibrium();
          debug_info << "== Starting with initial dt = " << dt << " until tEquilibirum " << tEquilibrium
              << " ==" << std::endl;
        }
        //solutionVectors.printHostGridVectors();
      }
      // ==============================================================================================

      // ========== Initial output ====================================================================
      Tools::compositeToChildCoefficientVector(multigfs, unew, unewCon, 0);
      uoldCon = unewCon;
      Tools::compositeToChildCoefficientVector(multigfs, unew, unewPot, 1);
      uoldPot = unewPot;

      //Output::printMultipleComponentCoefficientVector(unewCon, NUMBER_OF_SPECIES);
      //Output::printSingleCoefficientVector(uPot, "uPot");

      double debyeLength = physics.getDebyeLength(gfIonicStrength);
      // TODO Change this when adding a second extracellular domain
      double h_elec = (2 * physics.getParams().xMax() - physics.getParams().dMemb()) / (elecGV.size(0));

      debug_verb << std::endl;
      debug_verb << "@@@@@@@ minimum Debye length / length scale " << debyeLength / physics.getLengthScale()
          << std::endl;
      if(h_elec >= (debyeLength/physics.getLengthScale()))
      {
        debug_warn << "WARNING: Grid does not resolve Debye length [h_elec = " << h_elec << "]!" << std::endl;
      }
      debug_verb << "" << std::endl;

      typename Traits::ACME1MD_OUTPUT::DiagnosticInfo& diagInfo = acme1MDOutput.getDiagInfo();
      typename Traits::ACME1MD_OUTPUT::DiagnosticInfo lastDiagInfo(diagInfo); // copy of diagInfo


      Real dpot(0.0);
      Real old_dt = dt;
      Real last_dt = dt;
      Real last_dpot(0.0);

      diagInfo.tEquilibrium = tEquilibrium;
      diagInfo.addDebugData(std::string("abs max pot change"), 0.0);
      diagInfo.addDebugData(std::string("abs max con change"), 0.0);
      diagInfo.addDebugData(std::string("rel max pot change"), 0.0);
      diagInfo.addDebugData(std::string("rel max con change"), 0.0);

      diagInfo.addDebugData(std::string("old_dt"), old_dt);
      diagInfo.addDebugData(std::string("last_dt"), last_dt);
      diagInfo.addDebugData(std::string("dpot"), dpot);
      diagInfo.addDebugData(std::string("last_dpot"), last_dpot);
      diagInfo.addDebugData(std::string("dpot_dt"), 0.0);
      diagInfo.addDebugData(std::string("last_dpot_dt"), 0.0);

      diagInfo.addDebugData(std::string("rel_dpot_dt"), 0.0);


      acme1MDOutput.writeStep(time);
      Tools::pecletNumber(gv, parametersCon, dgfGradPot);
      debug_info  << std::endl << "########## initial output done" << std::endl << std::endl;
      // ==============================================================================================


      // ========= time loop ==========================================================================
      Real outputTimeInterval = physics.getParams().getOutputTimeInterval();
      int outputCounter = 1;
      bool printEveryTimeStep = (dt > outputTimeInterval && !physics.getParams().useAdaptiveTimeStep());

      double tInj_start = physics.getParams().tInj_start();
      double tInj_end   = physics.getParams().tInj_end();

      debug_jochen << "tInj_start = " << tInj_start << ", tInj_end = " << tInj_end << std::endl;

      if(physics.getParams().doStimulation())
      {
        assert(tInj_start > tEquilibrium);
      }

      typename Traits::U uChange = unew;
      typename Traits::U_CON uConChange(unewCon);
      typename Traits::U_POT uPotChange(unewPot);
      typename Traits::U_CON uConChange_Rel(unewCon);
      typename Traits::U_POT uPotChange_Rel(unewPot);

      // Total number of time steps
      int totalTimeSteps = 0;

      // Number of iterations/time steps since last output
      int iterations = 0;
      int timeSteps = 0;

      while (time<tend-1e-8)
      {
        debug_verb << "TIME = " << time << std::endl;
        debug_verb << "Last dt = " << dt << std::endl;

        last_dt = old_dt;
        old_dt = dt;

        // Use default time step when equilibration phase is over
        if(std::abs(time-tEquilibrium) < 1e-6)
        {
          dt = dtstart;
          bool printEveryTimeStep = (dt > outputTimeInterval && !physics.getParams().useAdaptiveTimeStep());
        }

        // Use default time step when stimulation starts
        if(time+dt > tInj_start && time < tInj_start)
        {
          debug_jochen << "== Injection starts, resetting dt to start value " << dtstart << std::endl;
          dt = dtstart;
        }


        bool acceptTimeStep = not physics.getParams().useAdaptiveTimeStep()
            || (time < tEquilibrium || std::abs(time - tEquilibrium) < 1e-6);

        typename GF_MEMB_FLUX::Traits::RangeType oldMaxFlux = gfMembFlux.getMaxFlux();
        // Update channels
        gfMembFlux.updateChannels(membGV, time, dt);
        gfMembFlux.updateFlux(membGV, time, dt);
        while(not acceptTimeStep)
        {
          Real dpot_dt = uPotChange.base().infinity_norm() / dt;
          Real last_dpot_dt = last_dpot / last_dt;

          typename GF_MEMB_FLUX::Traits::RangeType maxFlux = gfMembFlux.getMaxFlux();

          debug_jochen << "==========================================" << std::endl;
          debug_jochen << "dt = " << dt << std::endl;
          typename GF_MEMB_FLUX::Traits::RangeFieldType oldMaxTotalFlux(0.0);
          typename GF_MEMB_FLUX::Traits::RangeFieldType maxTotalFlux(0.0);
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            debug_jochen << "Last time step / this time step (flux*dt) = "
                << (oldMaxFlux[j] * dt) << " / " << (maxFlux[j] * dt) << std::endl;
            oldMaxTotalFlux += oldMaxFlux[j];
            maxTotalFlux += maxFlux[j];
          }
          // Factor representing the ratio of the current flux with respect to the one from the last time step
          typename GF_MEMB_FLUX::Traits::RangeFieldType fluxFactor = maxFlux.one_norm() / oldMaxFlux.one_norm();

          debug_jochen << std::endl << "Old / this (max flux * dt) = "
              << (oldMaxTotalFlux * dt) << " / " << (maxTotalFlux * dt)
              << " [factor " << fluxFactor << "]" << std::endl;

          // (1) 'Soft' criterion: Try to adjust time step according to number Newton iterations
          if(diagInfo.iterations < 3 && diagInfo.iterations <= lastDiagInfo.iterations)
          {
            dt *= 1.2; // Carefully increase time step
          }
          if(diagInfo.iterations >= 5)
          {
            dt /= 2;  // Use half time step when we need too many Newton iterations
          }

          // (2) 'Harder' criterion: Bound time step according to change in potential
          Real rel_dpot_dt = std::abs((dpot_dt-last_dpot_dt)/(std::max(dpot_dt,last_dpot_dt)));
          diagInfo.addDebugData(std::string("rel_dpot_dt"), rel_dpot_dt);
          /*
          if(time+dt > tInj_start && time+dt < tInj_end)
          {
            if(rel_dpot_dt > 0.05)
            {
              dt /= 2;
            }
            if(rel_dpot_dt < 0.001)
            {
              dt *= 1.2;
            }
          }
          */

          // (3) 'Hardest' criterion: Bound maximum/minimum time step
          dt = std::min(dt, 0.1e-3 / physics.getTimeScale()); // maximum 0.1 ms
          dt = std::max(dt, 1.0e-6 / physics.getTimeScale()); // minimum 1 µs

          typename DGF_POT::Traits::RangeType minPotJump = 1e100;
          for(typename PHYSICS::SubDomainElementIterator sdeit = membGV.template begin<0>();
              sdeit != membGV.template end<0>(); ++sdeit)
          {
            typename DGF_POT::Traits::RangeType potJump;
            physics.getMembranePotentialJump(membGV.grid().multiDomainEntity(*sdeit), dgfPot, potJump);
            minPotJump = std::min(minPotJump, potJump);
          }
          // Check if membrane potential is above -50mV
          if((time+dt > tInj_start && time+dt < tInj_end)
              //|| (time+dt > 15e3+tInj_start && time+dt < 15e3+tInj_end)
              || physics.convertTo_mV(minPotJump) < 50.)
          {
            dt = std::min(dt, 0.01e-3 / physics.getTimeScale()); // maximum 0.01 ms during AP
          }

          debug_jochen << "Last time step #iterations = " << diagInfo.iterations
              << " => new time step dt = " << dt << std::endl;

          debug_jochen << "==========================================" << std::endl;
          // Update channels
          gfMembFlux.updateChannels(membGV, time, dt);
          gfMembFlux.updateFlux(membGV, time, dt);
          acceptTimeStep = true;
        }
        gfMembFlux.acceptTimeStep(true);

        // Update values representing last iteration
        last_dpot = dpot;
        lastDiagInfo = diagInfo; // save diagnostic information for usage in next time step

        diagInfo.clear();
        // Hack to let PDE solver be able to get the current time
        diagInfo.time = (time+dt);

        // !! Time-dependent constraints are currently not implemented in PDELab !!
        // evaluate constraints for current time step
//        bctypePot.setTime(time+dt);
//        bctypeCon.setTime(time+dt);
//        ccCon.clear();
//        Dune::PDELab::constraints( bctypeCon, gfsCon, ccCon );
//        ccPot.clear();
//        Dune::PDELab::constraints( bctypePot, gfsPot, ccPot );

        physics.setTimeStep(dt);
        diagInfo.dt = dt;
        
        try{

          // ========= DO ONE TIME STEP ===========
          Acme1MDFullyCoupled::timeStep(time,dt,uold,unew,solver,multigfs,unewCon,unewPot,acme1MDOutput.getDiagInfo());
          uChange = unew;
          uChange -= uold;
          debug_jochen << "L2  change in solution in this time step: " << uChange.base().two_norm() << std::endl;
          debug_jochen << "MAX change in solution in this time step: " << uChange.base().infinity_norm() << std::endl;

          uold = unew;

        } catch(Dune::Exception& e)
        {

//          // Write debug output to file for later inspection
//          if(writeIfNotConverged)
//          {
//            debug_verb << "write debug output" << std::endl;
//            for(int i=0; i<diagInfo.maxDiffCon.size(); i++)
//            {
//              std::stringstream infoStream;
//              Real iteration = i;
//              infoStream << "iteration: " << std::fixed << std::setw(3) << i;
//              Output::gnuplotDoubleAppend("acme1MD_debug.dat", iteration, diagInfo.maxDiffCon[i], diagInfo.maxDiffPot[i],
//                  physics.getParams().getOutputPrefix(), infoStream.str());
//
//            }
//          }

          time += dt;
          // Write out non-converged solution
          acme1MDOutput.writeStep(time);
          Tools::pecletNumber(gv, parametersCon, dgfGradPot);
          //Output::printMultipleComponentCoefficientVector(uCon, NUMBER_OF_SPECIES);
          throw e;
        }

        uConChange = unewCon;
        uConChange -= uoldCon;
        uPotChange = unewPot;
        uPotChange -= uoldPot;

        dpot = uPotChange.base().infinity_norm();

        uPotChange_Rel = unewPot;
        int n = unewPot.flatsize();
        for(int i=0; i<n; ++i)
        {
          Traits::VBE::access(uPotChange_Rel, i) /= Traits::VBE::access(uoldPot, i);
          Traits::VBE::access(uPotChange_Rel, i) -= 1;

          //debug_jochen << uoldPot[i] << " --> " << unewPot[i] << " (" << uPotChange_Rel[i] << std::endl;
        }
        uConChange_Rel = unewCon;
        for(int i=0; i<uConChange_Rel.N(); i++)
        {
          uConChange_Rel[i] /= uoldCon[i];
          uConChange_Rel[i] -= 1;
        }
        {
          Dune::ios_base_all_saver hackbraten(std::cout);
          debug_info << std::scientific << std::setprecision(16);
          debug_info << "L2  con change  = " << uConChange.base().two_norm() << std::endl;
          debug_info << "MAX con change  = " << uConChange.base().infinity_norm() << std::endl;
          debug_info << "L2  pot change  = " << uPotChange.base().two_norm() << std::endl;
          debug_info << "MAX pot change  = " << uPotChange.base().infinity_norm() << std::endl;

          diagInfo.addDebugData(std::string("abs max pot change"), uPotChange.base().infinity_norm());
          diagInfo.addDebugData(std::string("rel max pot change"), (uPotChange_Rel.base().infinity_norm()));
          diagInfo.addDebugData(std::string("abs max con change"), uConChange.base().infinity_norm());
          diagInfo.addDebugData(std::string("old_dt"), old_dt);
          diagInfo.addDebugData(std::string("last_dt"), last_dt);
          diagInfo.addDebugData(std::string("dpot"), dpot);
          diagInfo.addDebugData(std::string("last_dpot"), last_dpot);
          diagInfo.addDebugData(std::string("rel max con change"), (uConChange_Rel.base().infinity_norm()));
          diagInfo.addDebugData(std::string("dpot_dt"), uPotChange.base().infinity_norm()/dt);
          diagInfo.addDebugData(std::string("last_dpot_dt"), last_dpot/last_dt);

        }

        uoldCon = unewCon;
        uoldPot = unewPot;

        time += dt;
        iterations += diagInfo.iterations;
        timeSteps++;
        totalTimeSteps++;

        double diffToNextTimeStep = time - outputCounter * outputTimeInterval;
        if (printEveryTimeStep || std::abs(diffToNextTimeStep) < 1e-8 || diffToNextTimeStep > 0)
        {
        	++outputCounter;

        	//Output::printSingleCoefficientVectorDG(uPot, "pot");
          //Output::printMultipleComponentCoefficientVectorDG(unewCon, NUMBER_OF_SPECIES);

        	acme1MDOutput.writeStep(time);
        	Tools::pecletNumber(gv, parametersCon, dgfGradPot);

        	debug_info << std::endl << "########## output done, time: " << time
        	    << " [" << timeSteps << " steps, average #iterations: " << (iterations/timeSteps) << "] ###########"
        	    << std::endl << std::endl;

        	iterations = 0;
        	timeSteps = 0;
        }


      }
      // ==============================================================================================

      debug_info << "Total number of time steps: " << totalTimeSteps << std::endl;

      // ======= Save simulation state ================================================================
      if(physics.getParams().doSaveState())
      {
        std::string saveFilename = physics.getParams().getSaveFilename();
        saveState(time,dt,solutionVectors,saveFilename);
        debug_info << "=============================================================" << std::endl;
        debug_info << "Saved simulation state to file " << saveFilename << std::endl;
        debug_info << "=============================================================" << std::endl;
      }
      // ==============================================================================================

    }

  private:
    void saveState(double time, double dt, Acme1MDSolutionVectors<Traits>& solutionVectors, std::string filename)
    {
      Ax1SimulationData<Real> simulationData(gv.size(0), time, dt, filename);

      std::vector<Real> stdVec;

      // MD grid
      solutionVectors.uold.std_copy_to(stdVec);
      simulationData.addVector("uold", stdVec);

      solutionVectors.unew.std_copy_to(stdVec);
      simulationData.addVector("unew", stdVec);

      // Electrolyte subdomain
      solutionVectors.uoldCon.std_copy_to(stdVec);
      simulationData.addVector("uoldCon", stdVec);

      solutionVectors.unewCon.std_copy_to(stdVec);
      simulationData.addVector("unewCon", stdVec);

      solutionVectors.uoldPot.std_copy_to(stdVec);
      simulationData.addVector("uoldPot", stdVec);

      solutionVectors.unewPot.std_copy_to(stdVec);
      simulationData.addVector("unewPot", stdVec);

      Ax1SimulationState<Real> state;
      state.setData(simulationData);
      state.saveState();
    }

    void loadState(double& time, double& dt, Acme1MDSolutionVectors<Traits>& solutionVectors, std::string filename)
    {
      Ax1SimulationState<Real> state;
      state.loadState(filename);

      Ax1SimulationData<Real>& simulationData = state.getData();

      if(gv.size(0) != simulationData.getNElements())
      {
        DUNE_THROW(Dune::Exception,
            "Error loading saved simulation state from file '" << filename
            << "', was the saved state generated with the same grid and finite element order?");
      }

      time = simulationData.getTime();
      dt = simulationData.getDt();

      std::vector<Real> stdVec;

      // Multidomain grid
      stdVec = simulationData.getVector("uold").data;
      solutionVectors.uold.std_copy_from(stdVec);

      stdVec = simulationData.getVector("unew").data;
      solutionVectors.unew.std_copy_from(stdVec);

      // Electrolyte subdomain
      stdVec = simulationData.getVector("uoldCon").data;
      solutionVectors.uoldCon.std_copy_from(stdVec);

      stdVec = simulationData.getVector("unewCon").data;
      solutionVectors.unewCon.std_copy_from(stdVec);

      stdVec = simulationData.getVector("uoldPot").data;
      solutionVectors.uoldPot.std_copy_from(stdVec);

      stdVec = simulationData.getVector("unewPot").data;
      solutionVectors.unewPot.std_copy_from(stdVec);

    }

    Grid& grid;

    GV& gv;
    PHYSICS& physics;

    SubGV& elecGV;
    SubGV& membGV;
};

#endif /* DUNE_AX1_ACME1MD_SETUP_HH */
