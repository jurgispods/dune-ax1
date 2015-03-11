#ifndef DUNE_AX1_ACME0_SETUP_HH
#define DUNE_AX1_ACME0_SETUP_HH

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
#include <dune/ax1/common/chargedensitygridfunction.hh>
#include <dune/ax1/common/ionicstrengthgridfunction.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/vectordiscretegridfunctiongradient.hh>
#include <dune/ax1/acme0/common/acme0_boundary.hh>
#include <dune/ax1/acme0/common/acme0_initial.hh>
#include <dune/ax1/acme0/common/acme0_output.hh>
#include <dune/ax1/acme0/common/convectiondiffusionfem.hh>
#include <dune/ax1/acme0/common/convectiondiffusiondg.hh>
#include <dune/ax1/acme0/common/nernst_planck_parameters.hh>
#include <dune/ax1/acme0/common/poisson_parameters.hh>
#include <dune/ax1/acme0/operator-split/acme0_operator_split.hh>
#include <dune/ax1/acme0/operator-split/acme0_toperator.hh>
#include <dune/ax1/acme0/operator-split/nernst_planck_power_operator.hh>
#include <dune/ax1/acme0/fully-implicit/acme0_fully_coupled.hh>
#include <dune/ax1/acme0/fully-implicit/acme0_operator_fully_implicit.hh>
#include <dune/ax1/acme0/fully-implicit/acme0_experimental_operator_fully_implicit.hh> // fully fully-coupled!
#include <dune/ax1/acme0/fully-implicit/acme0_toperator_fully_implicit.hh>


#define USE_CG 1

template<class GV, class PHYSICS>
class Acme0Setup
{
  public:
    
    // Choose domain and range field type
    typedef GV GridView;
    typedef typename GV::Grid::ctype Coord;
    typedef double Real;
    static const bool useRowPreconditioner = false;
    
    const bool USE_OPERATOR_SPLIT;

    
    Acme0Setup(const GV& gv_, PHYSICS& physics_, bool opSplit_ = true) :
      gv(gv_),
      physics(physics_),
			USE_OPERATOR_SPLIT(opSplit_)
    {
#if USE_CG == 1
      if(physics_.getParams().useMembrane())
      {
        DUNE_THROW(Dune::Exception, "acme0 cannot be used with a membrane when using continuous finite elements!");
      }
#endif
    }


    void setup (double dtstart, double tend)
    {
      const int dim = GV::dimension;
      Real time = 0.0;
      const int degree=1;

      // We could as well choose NoConstraints for DG, as these are ignored anyway
      // => Use ConvectionDiffusion parameter class for setting boundary values!
      typedef Dune::PDELab::ConformingDirichletConstraints CONSTRAINTS; // constraints class (Dirichlet)
      typedef Dune::PDELab::ISTLVectorBackend<1> VBE;                   //vector backend

#if USE_CG==1
			// Finite element map for CG
      typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM_CON;
      FEM_CON femCon(degree);
      typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM_POT;
      FEM_POT femPot(degree);
#else
      // Finite element map for DG
      typedef Dune::PDELab::OPBLocalFiniteElementMap<Coord,Real,degree,dim,Dune::GeometryType::simplex> FEM_CON;
      FEM_CON femCon;
      typedef Dune::PDELab::OPBLocalFiniteElementMap<Coord,Real,degree,dim,Dune::GeometryType::simplex> FEM_POT;
      FEM_POT femPot;
#endif

      // =================== GFS SETUP ================================================================
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM_CON,CONSTRAINTS,VBE> GFS_SINGLE_CON;
      GFS_SINGLE_CON gfsSingleCon(gv,femCon);


      // Power grid function space for all ion species
      typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper PowerGFSOrdering;
      typedef Dune::PDELab::PowerGridFunctionSpace<
        GFS_SINGLE_CON, NUMBER_OF_SPECIES, PowerGFSOrdering> GFS_CON;
      GFS_CON gfsCon(gfsSingleCon);

      // Create grid function space for potential
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM_POT,CONSTRAINTS,VBE> GFS_POT;
      GFS_POT gfsPot(gv,femPot);

      // Now wrap grid function spaces into composite GFS
      //typedef Dune::PDELab::GridFunctionSpaceComponentBlockwiseMapper<3,1> CompositeGFSOrdering;
      //typedef Dune::PDELab::GridFunctionSpaceBlockwiseMapper CompositeGFSOrdering;
      typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper CompositeGFSOrdering;
      typedef Dune::PDELab::CompositeGridFunctionSpace<
          CompositeGFSOrdering,GFS_CON,GFS_POT> GFS;
      GFS gfs(gfsCon,gfsPot);

      typedef Dune::PDELab::GridFunctionSubSpace<GFS,0> GFS_SUB_CON;
      GFS_SUB_CON gfsSubCon(gfs);
      typedef Dune::PDELab::GridFunctionSubSpace<GFS,1> GFS_SUB_POT;
      GFS_SUB_POT gfsSubPot(gfs);

      // Extract solution vector types
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,Real>::Type U;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CON,Real>::Type U_CON;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_POT,Real>::Type U_POT;
      
      // Types for different grid functions
#ifdef ACME0_FULLY_IMPLICIT
      typedef Dune::PDELab::VectorDiscreteGridFunction  <GFS_SUB_CON,U> DGF_CON;
      typedef Dune::PDELab::DiscreteGridFunction<GFS_SUB_POT,U> DGF_POT;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_SUB_POT,U> DGF_POT_GRAD;
#else
      typedef Dune::PDELab::VectorDiscreteGridFunction  <GFS_CON,U_CON> DGF_CON;
      typedef Dune::PDELab::DiscreteGridFunction<GFS_POT,U_POT> DGF_POT;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_POT,U_POT> DGF_POT_GRAD;
#endif
      typedef ChargeDensityGridFunction<DGF_CON, PHYSICS> GF_CD;
      typedef IonicStrengthGridFunction<DGF_CON, PHYSICS> GF_IS;
      // ==============================================================================================

      // =========== Define solution vectors ==========================================================
      U_CON uoldCon(gfsCon,0.0);
      U_CON unewCon(gfsCon,0.0);
      U_POT uPot(gfsPot,0.0);
      U uold(gfs,0.0);
      U unew(gfs,0.0);
      // ==============================================================================================

      // ====== Create grid functions and interpolate initial values onto coefficient vectors =========
      // Inital concentrations
      typedef typename PHYSICS::Traits::INITIAL_CON INITIAL_CON;
      INITIAL_CON initialCon(gv,physics.getParams());
      Dune::PDELab::interpolate(initialCon,gfsCon,uoldCon);
      unewCon = uoldCon;

#ifdef ACME0_FULLY_IMPLICIT
      debug_verb << "@@@  fully-coupled case! -> use special DGF for concentrations" << std::endl;
      // coupled
      DGF_CON dgfCon(gfsSubCon, unew);

      // decoupled
      //DGF_CON dgfCon(gfsSubCon, uold);
#else
      DGF_CON dgfCon(gfsCon, unewCon);
#endif
      GF_CD gfChargeDensity(dgfCon, dgfCon, physics);
      GF_IS gfIonicStrength(dgfCon, physics);

      /*
      Output::printMultipleComponentCoefficientVectorDG(unewCon, NUMBER_OF_SPECIES);
      std::valarray<Real> conValArray;
      Tools::getSolutionVectorDG(dgfCon, 0, physics.getPosition(), conValArray);
      for(int i=0; i<conValArray.size(); ++i)
      {
        debug_verb << conValArray[i] << std::endl;
      }
      */

      // Initial potential
      typedef InitialPot<GV,Real,PHYSICS,GF_CD> INITIAL_POT;
      INITIAL_POT initialPot(gv,physics,gfChargeDensity);
      initialPot.setTime(time);
      Dune::PDELab::interpolate(initialPot,gfsPot,uPot);


#ifdef ACME0_FULLY_IMPLICIT
      debug_verb << "@@@  fully-coupled case! -> use special DGF for potential" << std::endl;
      // coupled
      DGF_POT_GRAD dgfGradPot(gfsSubPot, unew);
      DGF_POT      dgfPot(gfsSubPot, unew);

      // decoupled
      //DGF_POT_GRAD dgfGradPot(gfsSubPot, uold);
      //DGF_POT      dgfPot(gfsSubPot, uold);
#else
      DGF_POT      dgfPot(gfsPot, uPot);
      DGF_POT_GRAD dgfGradPot(gfsPot, uPot);
#endif
      
      typedef Dune::PDELab::CompositeGridFunction<INITIAL_CON,INITIAL_POT> INITIAL_U;
      INITIAL_U initialU(initialCon,initialPot);
      Dune::PDELab::interpolate(initialU,gfs,uold);
      unew = uold;

      // Define and instantiate output class
      typedef Acme0Output<GFS_CON,GFS_POT,U_CON,U_POT,PHYSICS> Acme0Output;
      Acme0Output acme0Output(gfsCon, gfsPot, unewCon, uPot, physics);
      // ==============================================================================================

      // ========== Define Parameter classes containing the model problem =============================
      typedef typename PHYSICS::Traits::NERNST_PLANCK_BOUNDARY BOUNDARY_CON;
      BOUNDARY_CON boundaryCon(initialCon);
      typedef NernstPlanckParameters<GV,Real,PHYSICS,DGF_POT_GRAD,BOUNDARY_CON> PARAMETERS_CON;
      PARAMETERS_CON parametersCon(physics,dgfGradPot,boundaryCon);

      typedef typename PHYSICS::Traits::POISSON_BOUNDARY BOUNDARY_POT;
      BOUNDARY_POT boundaryPot;
      typedef PoissonParameters<GV,Real,PHYSICS,GF_CD,INITIAL_POT,BOUNDARY_POT> PARAMETERS_POT;
      PARAMETERS_POT parametersPot(physics,gfChargeDensity,initialPot,boundaryPot);
      // ==============================================================================================

      // =================== BOUNDARY CONDITIONS ======================================================
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<PARAMETERS_CON> BCType_SINGLE_CON;
      typedef Dune::PDELab::PowerConstraintsParameters<BCType_SINGLE_CON,NUMBER_OF_SPECIES> BCType_CON;
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<PARAMETERS_POT> BCType_POT;
      typedef Dune::PDELab::CompositeConstraintsParameters<BCType_CON,BCType_POT> BCType;

      // Attention! With this construction, all species have the same BCType.
      // We cannot distinguish between species when using ConvectionDiffusionBoundaryConditionAdapter!
      BCType_SINGLE_CON bctypeSingleCon(gv,parametersCon);
      BCType_CON bctypeCon(bctypeSingleCon);
      BCType_POT bctypePot(gv,parametersPot);
      BCType bctype(bctypeCon,bctypePot);

      typedef typename GFS_CON::template ConstraintsContainer<Real>::Type CC_CON;
      CC_CON ccCon;
      typedef typename GFS_POT::template ConstraintsContainer<Real>::Type CC_POT;
      CC_POT ccPot;
      typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
      CC cc;

      Dune::PDELab::constraints( bctypePot, gfsPot, ccPot); // assemble potential constraints
      Dune::PDELab::constraints( bctypeCon, gfsCon, ccCon); // assemble concentration constraints
      Dune::PDELab::constraints( bctype, gfs, cc); // assemble all constraints

      debug_info << "constrained potential dofs="     << ccPot.size() << " of " << gfsPot.globalSize() << std::endl;
      debug_info << "constrained concentration dofs=" << ccCon.size() << " of " << gfsCon.globalSize() << std::endl;
      debug_info << "constrained dofs=" << cc.size() << " of " << gfs.globalSize() << std::endl;

      // Attention! With this construction, all species have the same Dirichlet values.
      // We cannot distinguish between species when using ConvectionDiffusionDirichletExtensionAdapter!
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PARAMETERS_CON> DirichletValuesSingleCon;
      DirichletValuesSingleCon dirichletValuesSingleCon(gv,parametersCon);
      typedef Dune::PDELab::PowerGridFunction<DirichletValuesSingleCon,NUMBER_OF_SPECIES> DirichletValuesCon;
      DirichletValuesCon dirichletValuesCon(dirichletValuesSingleCon);
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PARAMETERS_POT> DirichletValuesPot;
      DirichletValuesPot dirichletValuesPot(gv,parametersPot);

      Dune::PDELab::interpolate(dirichletValuesPot,gfsPot,uPot);
      // ==============================================================================================

      // ============== Make grid operator space =========================================
      // ##### Spatial part #####
#if USE_CG==1
      // ### Standard FEM (CG) Nernst-Planck local operator
      typedef Dune::PDELab::ConvectionDiffusionFEM<PARAMETERS_CON,FEM_CON> LOP_SINGLE_CON;
      LOP_SINGLE_CON lopSingleCon(parametersCon);
#else
      // ### DG Nernst-Planck local operator
      Dune::PDELab::ConvectionDiffusionDGMethod::Type mCon;
      mCon = Dune::PDELab::ConvectionDiffusionDGMethod::SIPG;
      Dune::PDELab::ConvectionDiffusionDGWeights::Type wCon;
      wCon = Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn;
      Real alphaCon = 2.0; // 2.0
      typedef Dune::PDELab::ConvectionDiffusionDG<PARAMETERS_CON,FEM_CON> LOP_SINGLE_CON;
      LOP_SINGLE_CON lopSingleCon(parametersCon,mCon,wCon,alphaCon);
#endif
      typedef NernstPlanckDGPowerOperator<PARAMETERS_CON,FEM_CON,LOP_SINGLE_CON> LOP_CON;
      LOP_CON lopCon(parametersCon,lopSingleCon);
	    
#if USE_CG==1
      // ### Standard FEM (CG) Poisson local operator
      typedef Dune::PDELab::ConvectionDiffusionFEM<PARAMETERS_POT,FEM_POT> LOP_POT;
      LOP_POT lopPot(parametersPot);
#else
	    // ### DG Poisson local operator
      Dune::PDELab::ConvectionDiffusionDGMethod::Type mPot;
      mPot = Dune::PDELab::ConvectionDiffusionDGMethod::SIPG;
      Dune::PDELab::ConvectionDiffusionDGWeights::Type wPot;
      wPot = Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn;
      Real alphaPot = 2.0; // 2.0
      typedef Dune::PDELab::ConvectionDiffusionDG<PARAMETERS_POT,FEM_POT> LOP_POT;
    	LOP_POT lopPot(parametersPot,mPot,wPot,alphaPot);
#endif
      
      // ### Combined local operator for fully coupled system
      //typedef Acme0OperatorFullyImplicit<PARAMETERS_CON,PARAMETERS_POT,FEM_CON,FEM_POT,LOP_SINGLE_CON,LOP_POT> LOP;
      //LOP lop(parametersCon,parametersPot,lopSingleCon,lopPot);
    	typedef Acme0ExperimentalOperatorFullyImplicit<PARAMETERS_CON,PARAMETERS_POT,FEM_CON,FEM_POT> LOP;
      LOP lop(parametersCon,parametersPot);
      
      // ##### Temporal part #####
      typedef NernstPlanckTimeLocalOperator<PHYSICS> TLOP_CON;
      TLOP_CON tlopCon(physics, 2);
      typedef Acme0TimeLocalOperator TLOP;
    	TLOP tlop(2);

      typedef VBE::MatrixBackend MBE;
      
      typedef Dune::PDELab::GridOperator<GFS_POT,GFS_POT,LOP_POT,MBE,Real,Real,Real,CC_POT,CC_POT> GO_POT;
      GO_POT goPot(gfsPot,ccPot,gfsPot,ccPot,lopPot);
      typedef Dune::PDELab::GridOperator<GFS_CON,GFS_CON,LOP_CON,MBE,Real,Real,Real,CC_CON,CC_CON> GO_CON0;
      GO_CON0 goCon0(gfsCon,ccCon,gfsCon,ccCon,lopCon);
      typedef Dune::PDELab::GridOperator<GFS_CON,GFS_CON,TLOP_CON,MBE,Real,Real,Real,CC_CON,CC_CON> GO_CON1;
      GO_CON1 goCon1(gfsCon,ccCon,gfsCon,ccCon,tlopCon);
      typedef Dune::PDELab::OneStepGridOperator<GO_CON0,GO_CON1> GO_CON;
      GO_CON goCon(goCon0,goCon1);

      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO0;
      GO0 go0(gfs,cc,gfs,cc,lop);
      typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,Real,Real,Real,CC,CC> GO1;
      GO1 go1(gfs,cc,gfs,cc,tlop);
      typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
      IGO igo(go0,go1);

      // ### CG Poisson local operator
      /*
      typedef PoissonOperator<BCT_POT,PHYSICS,DGF_CON> LOP_POT;
      LOP_POT lopPot(bctPot, physics, dgfCon, 2);    // spatial part
      //typedef Dune::PDELab::GridOperator<GFS_POT,GFS_POT,LOP_POT,MBE,Real,Real,Real,CC_POT,CC_POT> GO_POT;
      //GO_POT goPot(gfsPot,ccPot,gfsPot,ccPot,lopPot);
       */
      // ==============================================================================================


      // ========== Select a linear solver backend ====================================================
      //typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS_CON;
      //typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS_POT;
      //LS_CON lsCon(5000,false);
      //LS_POT lsPot(5000,false);
      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
      //LS ls(5000,false);
      typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS_CON;
      LS_CON lsCon(0); // verbose = 1
      typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS_POT;
      LS_POT lsPot(0); // verbose = 1
      typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
      LS ls(0); // verbose = 1
      // ==============================================================================================

      // ========= Solver for (non-)linear problem per stage ============================================
      // Poisson solver
      typedef Ax1StationaryLinearProblemSolver<GO_POT,LS_POT,U_POT> SOLVER_POT;
      SOLVER_POT solverPot(goPot,lsPot,1e-10,"octave/poisson");
      solverPot.setPrintMatrix(false);
      solverPot.setPrintRhs(false);
      solverPot.setRowPreconditioner(useRowPreconditioner);
      solverPot.setVerbosity(0); //2
      debug_verb << "Poisson PDE solver set up" << std::endl;
      //debug_verb << "============= Solving initial Poisson problem" << std::endl;
      //solverPot.apply(uPot);
      //debug_verb << "============= SOLVED initial Poisson problem" << std::endl;

      /*
      // ============ Compare with analytical/exact solution ===========================================
      //std::valarray<Real> x = physics.getCellCenterPositions();
      std::valarray<Real> pos;

      // Numerical solution
      std::valarray<Real> pot;
      //Tools::getSolutionVectorDG(dgfPot, 0, pos, pot);
      Tools::getSolutionVectorDGmean(dgfPot, pos, pot);

      // Analytical solution
      typedef typename PHYSICS::Traits::SOLUTION_CON ANALYTICAL_SOLUTION_CON;
      typedef typename PHYSICS::Traits::SOLUTION_POT ANALYTICAL_SOLUTION_POT;
      ANALYTICAL_SOLUTION_CON gfAnalyticalSolutionCon(gfsCon.gridview(),physics.getParams());
      ANALYTICAL_SOLUTION_POT gfAnalyticalSolutionPot(gfsPot.gridview(),physics.getParams());
      std::valarray<Real> potAnal;
      //Tools::getSolutionVectorDG(gfAnalyticalSolutionPot, 0, pos, potAnal);
      Tools::getSolutionVectorDGmean(gfAnalyticalSolutionPot, pos, potAnal);
      std::valarray<Real> naAnal;
      //Tools::getSolutionVectorDG(gfAnalyticalSolutionCon, 0, pos, naAnal);
      //Tools::getSolutionVectorDGmean(gfAnalyticalSolutionCon, pos, naAnal);


      // Exact solution
      std::valarray<Real> potExact;
      initialPot.setCalculateAll(true);
      //Tools::getSolutionVectorDG(initialPot, 0, pos, potExact);
      Tools::getSolutionVectorDGmean(initialPot, pos, potExact);
      std::valarray<Real> chargeDensity;
      Output::gnuplotInitialize("pot_exact.dat", physics.getParams().getOutputPrefix(),"# pot exact initial");
      Output::gnuplotAppendArray("pot_exact.dat", pos, physics.convertTo_mV(potExact),
      std::string("bla"), physics.getParams().getOutputPrefix());

      //Tools::getSolutionVectorDG(gfChargeDensity, 0, pos, chargeDensity);
      Tools::getSolutionVectorDGmean(gfChargeDensity, pos, chargeDensity);
      for(int i=0; i<potExact.size(); ++i)
      {
        //debug_verb << "x = " << pos[i] << ",";
        if(physics.getParams().hasAnalyticalSolution())
        {
          //debug_verb << "naAnal = " << naAnal[i] << ", ";
          //debug_verb << "potAnal = " << potAnal[i] << ", ";
        }
        //debug_verb << "potExact = " << potExact[i] << ", pot = " << pot[i] << ", diff = "
        //    << (potExact[i] - pot[i]) << " @ chargeDensity = " << chargeDensity[i] << std::endl;
      }
      initialPot.setCalculateAll(false);
      // =============================================================================================
      */

      typedef Ax1StationaryLinearProblemSolver<GO_CON,LS_CON,U_CON> PDESOLVER_CON;
      PDESOLVER_CON pdesolverCon(goCon,lsCon,1e-10,"octave/nernst_planck");
      pdesolverCon.setPrintMatrix(false);
      pdesolverCon.setPrintRhs(false);
      pdesolverCon.setRowPreconditioner(useRowPreconditioner);
      pdesolverCon.setVerbosity(0); //2
      debug_verb << "Nernst-Planck PDE solver set up" << std::endl;


      typedef Ax1Newton<Acme0Output,IGO,LS,U> PDESOLVER;
      PDESOLVER pdesolver(acme0Output,igo,ls,"octave/acme0-fully-implicit_mat");
      pdesolver.setReassembleThreshold(0.0);
      pdesolver.setVerbosityLevel(2); // 2
      pdesolver.setReduction(physics.getParams().getReduction()); // 1e-10
      pdesolver.setAbsoluteLimit(physics.getParams().getAbsLimit()); // 1e-10
      pdesolver.setMinLinearReduction(1e-4);
      pdesolver.setMaxIterations(100);
      //pdesolver.setLineSearchStrategy(PDESOLVER::noLineSearch);
      pdesolver.setLineSearchMaxIterations(10); // 10
      pdesolver.setPrintMatrix(false);
      pdesolver.setPrintRhs(false);
      pdesolver.setRowPreconditioner(useRowPreconditioner);

      /*typedef Dune::PDELab::CompositeGridFunction<INITIAL_CON,DGF_POT> INITIAL_U;
      INITIAL_U initialU(initialCon,dgfPot);
      Dune::PDELab::interpolate(initialU,gfs,uold);*/

      // !! Copy values of child coefficient vectors to parent coeffcient vector !!
      //Tools::childToCompositeCoefficientVector(gfs, unewCon, 0, uold);
      Tools::childToCompositeCoefficientVector(gfs, uPot, 1, uold);
      unew = uold;
      //Output::printCoefficientVectorDG(unew, NUMBER_OF_SPECIES+1);
      // ==============================================================================================
    
    
      // ========== time-stepper ======================================================================
      //Dune::PDELab::Alexander2Parameter<Real> timeStepperCon;
      Dune::PDELab::ImplicitEulerParameter<Real> timeStepperCon;
      //Dune::PDELab::ExplicitEulerParameter<Real> timeStepperCon;
      //Dune::PDELab::RK4Parameter<Real> timeStepperCon;
      //Dune::PDELab::HeunParameter<Real> timeStepperCon;
      typedef Dune::PDELab::OneStepMethod<Real,GO_CON,PDESOLVER_CON,U_CON,U_CON> SOLVER_CON;
      SOLVER_CON solverCon(timeStepperCon,goCon,pdesolverCon);
      solverCon.setVerbosityLevel(0); // 2
      
      //Dune::PDELab::Alexander2Parameter<Real> method;
      Dune::PDELab::ImplicitEulerParameter<Real> timeStepper;
      //Dune::PDELab::RK4Parameter<Real> timeStepper;
      typedef Dune::PDELab::OneStepMethod<Real,IGO,PDESOLVER,U,U> SOLVER;
      SOLVER solver(timeStepper,igo,pdesolver);
      solver.setVerbosityLevel(0); // 2

      const PDESOLVER_CON& pdesolverConRef = solverCon.getPDESolver();
      const PDESOLVER&     pdesolverRef    = solver.getPDESolver();
      // ==============================================================================================
      
      // ========== Initial output ====================================================================
      double debyeLength = physics.getDebyeLength(gfIonicStrength);
      debug_verb << "" << std::endl;
      debug_verb << "@@@@@@@ minimum Debye length / length scale " << debyeLength / physics.getLengthScale()<< std::endl;
      debug_verb << "" << std::endl;

      //Tools::compositeToChildCoefficientVector(gfs, unew, uCon, 0);
      //Tools::compositeToChildCoefficientVector(gfs, unew, uPot, 1);
      //Output::printSingleCoefficientVectorDG(uPot, "pot");
      //Output::printMultipleComponentCoefficientVectorDG(unewCon, NUMBER_OF_SPECIES);
      acme0Output.writeStep(time);
      Tools::pecletNumber(gv, parametersCon, dgfGradPot);
      debug_info  << std::endl << "########## initial output done" << std::endl << std::endl;
      // ==============================================================================================

      // ========= time loop ==========================================================================
      // Setup solution method
      //typedef Acme0OperatorSplit<U_CON,U_POT,DirichletValuesCon,DirichletValuesPot,SOLVER_CON,
      typedef Acme0OperatorSplit<U_CON,U_POT,INITIAL_CON,INITIAL_POT,SOLVER_CON,
                          SOLVER_POT,GFS_CON,GFS_POT,PHYSICS,GV,GO_CON,GO_POT,Real> Acme0OperatorSplit;
      //Acme0OperatorSplit acme0OperatorSplit(uoldCon, unewCon, uPot, dirichletValuesCon, dirichletValuesPot,
      Acme0OperatorSplit acme0OperatorSplit(uoldCon, unewCon, uPot, initialCon, initialPot,
                         solverCon, solverPot, gfsCon, gfsPot, physics, gv, goCon, goPot);
      acme0OperatorSplit.setReduction(physics.getParams().getReduction());
      acme0OperatorSplit.setAbsLimit(physics.getParams().getAbsLimit());
      acme0OperatorSplit.setMaxIterations(100);
      acme0OperatorSplit.setUseDefect(false);

      Real dt = dtstart;
      Real outputTimeInterval = physics.getParams().getOutputTimeInterval();
      int outputCounter = 1;
      bool printEveryTimeStep = (dt > outputTimeInterval && !physics.getParams().useAdaptiveTimeStep());

      while (time<tend-1e-8)
      {
        //debug_verb << "TIME= " << time << std::endl;
        //debug_verb << "outputCounter= " << outputCounter << std::endl;

        //


        // Hack to let PDE solver be able to get the current time
        acme0Output.getDiagInfo().time = (time+dt);

        // !! Time-dependent constraints are currently not implemented for ConvectionDiffusionBoundaryConditionAdapter !!
        // evaluate constraints for current time step
        // (This doesn't do anything for DG, as there are no 'real' constraints)
//        if(USE_OPERATOR_SPLIT)
//        {
//          bctypePot.setTime(time+dt);
//          bctypeCon.setTime(time+dt);
//	        ccCon.clear();
//  	      Dune::PDELab::constraints( bctypeCon, gfsCon, ccCon );
//    	    ccPot.clear();
//      	  Dune::PDELab::constraints( bctypePot, gfsPot, ccPot );
//    	  } else {
//    	    bctype.setTime(time+dt);
//          cc.clear();
//          Dune::PDELab::constraints(bctype,gfs,cc);
//        }

        if(physics.getParams().useAdaptiveTimeStep())
        {
          Real safety_factor = 0.1;
          dt = safety_factor * Tools::getTimeStep(gv, physics, dgfGradPot);
          debug_verb << "Calculated time step: " << dt << std::endl;
        }
        physics.setTimeStep(dt);
        acme0Output.getDiagInfo().dt = dt;
        
        if(USE_OPERATOR_SPLIT)
        {
        	acme0OperatorSplit.timeStep(time,dt,acme0Output.getDiagInfo());
        	uoldCon = unewCon;
        } else {
          Acme0FullyCoupled::timeStep(time, dt, uold, initialU, unew, solver,
                                      gfs, unewCon, uPot, acme0Output.getDiagInfo());
          U uChange = unew;
          uChange -= uold;
          debug_jochen << "L2  change in solution in this time step: " << uChange.base().two_norm() << std::endl;
          debug_jochen << "MAX change in solution in this time step: " << uChange.base().infinity_norm() << std::endl;

	        uold = unew;
	      }
        
        // Permanently use new vector of gating particles
        physics.getMembrane().updateState();
        time += dt;

        //debug_verb << "time = " << time << std::endl;
        double diffToNextTimeStep = time - outputCounter * outputTimeInterval;
        if (printEveryTimeStep || std::abs(diffToNextTimeStep) < 1e-8 || diffToNextTimeStep > 0)
        {
        	++outputCounter;

        	//Output::printSingleCoefficientVectorDG(uPot, "pot");
          //Output::printMultipleComponentCoefficientVectorDG(unewCon, NUMBER_OF_SPECIES);

        	acme0Output.writeStep(time);
        	Tools::pecletNumber(gv, parametersCon, dgfGradPot);
        	acme0Output.getDiagInfo().clear();

        	debug_info
			      << std::endl
        	  << "########## output done, time: " << time << " ###########"
        	  << std::endl << std::endl;
        }
      }
      // ==============================================================================================
    }



  private:
    const GV& gv;
    PHYSICS& physics;

};

#endif /* DUNE_AX1_ACME0_SETUP_HH */
