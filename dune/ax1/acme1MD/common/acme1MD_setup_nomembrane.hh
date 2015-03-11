#ifndef DUNE_AX1_ACME1MD_SETUP_NOMEMBRANE_HH
#define DUNE_AX1_ACME1MD_SETUP_NOMEMBRANE_HH

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
#include <dune/ax1/common/chargedensitygridfunction.hh>
#include <dune/ax1/common/ionicstrengthgridfunction.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/vectordiscretegridfunctiongradient.hh>
#include <dune/ax1/acme1MD/common/acme1MD_boundary.hh>
#include <dune/ax1/acme1MD/common/acme1MD_initial.hh>
#include <dune/ax1/acme1MD/common/acme1MD_output.hh>
#include <dune/ax1/acme1MD/common/convectiondiffusionfem.hh>
#include <dune/ax1/acme1MD/common/convectiondiffusiondg.hh>
#include <dune/ax1/acme1MD/common/nernst_planck_parameters.hh>
#include <dune/ax1/acme1MD/common/poisson_parameters.hh>
#include <dune/ax1/acme1MD/operator-split/acme1MD_operator_split_nomembrane.hh>
#include <dune/ax1/acme1MD/operator-split/acme1MD_toperator.hh>
#include <dune/ax1/acme1MD/operator-split/nernst_planck_power_operator.hh>

#define USE_CG 1

template<class GV, class PHYSICS>
class Acme1MDSetup
{
  public:
    
    // Choose domain and range field type
    typedef GV GridView;
    typedef typename GV::Grid::ctype Coord;
    typedef double Real;
    static const bool useRowPreconditioner = false;

    Acme1MDSetup(const GV& gv_, PHYSICS& physics_) :
      gv(gv_),
      physics(physics_)
    {}


    void setup (double dtstart, double tend)
    {
      const int dim = GV::dimension;
      Real time = 0.0;
      Real tEquilibrium = 1e3;
      const int degreeCon=1;
      const int degreePot=2;

      // We could as well choose NoConstraints for DG, as these are ignored anyway
      // => Use ConvectionDiffusion parameter class for setting boundary values!
      typedef Dune::PDELab::ConformingDirichletConstraints CONSTRAINTS; // constraints class (Dirichlet)
      typedef Dune::PDELab::ISTLVectorBackend<1> VBE;                   //vector backend

			// Finite element map for CG
      typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM_CON;
      FEM_CON femCon(degreeCon);
      typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM_POT;
      FEM_POT femPot(degreePot);

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

      // Extract solution vector types
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CON,Real>::Type U_CON;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_POT,Real>::Type U_POT;
      
      // Types for different grid functions
      typedef Ax1VectorDiscreteGridFunction  <PHYSICS,GFS_CON,U_CON> DGF_CON;
      typedef Dune::PDELab::DiscreteGridFunction<GFS_POT,U_POT> DGF_POT;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_POT,U_POT> DGF_POT_GRAD;
      typedef ChargeDensityGridFunction<DGF_CON, PHYSICS> GF_CD;
      typedef IonicStrengthGridFunction<DGF_CON, PHYSICS> GF_IS;
      // ==============================================================================================

      // =========== Define solution vectors ==========================================================
      U_CON uoldCon(gfsCon,0.0);
      U_CON unewCon(gfsCon,0.0);
      U_POT uPot(gfsPot,0.0);
      // ==============================================================================================

      // ====== Create grid functions and interpolate initial values onto coefficient vectors =========
      // Inital concentrations
      typedef typename PHYSICS::Traits::INITIAL_CON INITIAL_CON;
      INITIAL_CON initialCon(gv,physics.getParams());
      Dune::PDELab::interpolate(initialCon,gfsCon,uoldCon);
      unewCon = uoldCon;

      DGF_CON dgfCon(physics, gfsCon, unewCon);
      GF_CD gfChargeDensity(dgfCon, dgfCon, physics);
      GF_IS gfIonicStrength(dgfCon, physics);

      // Initial potential
      typedef InitialPot<GV,Real,PHYSICS,GF_CD> INITIAL_POT;
      INITIAL_POT initialPot(gv,physics,gfChargeDensity,2*degreePot);
      initialPot.setTime(time);
      Dune::PDELab::interpolate(initialPot,gfsPot,uPot);

      DGF_POT      dgfPot(gfsPot, uPot);
      DGF_POT_GRAD dgfGradPot(gfsPot, uPot);
      
      // Define and instantiate output class
      typedef Acme1MDOutput<GFS_CON,GFS_POT,U_CON,U_POT,PHYSICS> Acme1MDOutput;
      Acme1MDOutput acme1MDOutput(gfsCon, gfsPot, unewCon, uPot, physics, 2*degreePot, 2*degreeCon);
      // ==============================================================================================

      // ========== Define Parameter classes containing the model problem =============================
      typedef typename PHYSICS::Traits::NERNST_PLANCK_BOUNDARY BOUNDARY_CON;
      BOUNDARY_CON boundaryCon(initialCon);
      typedef NernstPlanckParameters<GV,Real,PHYSICS,DGF_POT_GRAD,BOUNDARY_CON> PARAMETERS_CON;
      PARAMETERS_CON parametersCon(gv,physics,dgfGradPot,boundaryCon,tEquilibrium);

      typedef typename PHYSICS::Traits::POISSON_BOUNDARY BOUNDARY_POT;
      BOUNDARY_POT boundaryPot;
      typedef PoissonParameters<GV,Real,PHYSICS,GF_CD,BOUNDARY_POT> PARAMETERS_POT;
      PARAMETERS_POT parametersPot(physics,gfChargeDensity,boundaryPot);
      // ==============================================================================================

      // =================== BOUNDARY CONDITIONS ======================================================
      // ###### BCType ########
      typedef BCTypeSingleCon<PARAMETERS_CON> BCType_SINGLE_CON;
      typedef Dune::PDELab::PowerConstraintsParameters<BCType_SINGLE_CON,NUMBER_OF_SPECIES> BCType_CON;
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<PARAMETERS_POT> BCType_POT;

      // Create NUMBER_OF_SPECIES boundary condition type classes
      BCType_CON bctypeCon;
      for(int i=0; i<NUMBER_OF_SPECIES; ++i)
      {
        BCType_SINGLE_CON* bctypeSingleCon = new BCType_SINGLE_CON(gv,parametersCon,i);
        bctypeCon.setChild(i,*bctypeSingleCon);
      }
      BCType_POT bctypePot(gv,parametersPot);

      typedef typename GFS_CON::template ConstraintsContainer<Real>::Type CC_CON;
      CC_CON ccCon;
      typedef typename GFS_POT::template ConstraintsContainer<Real>::Type CC_POT;
      CC_POT ccPot;

      Dune::PDELab::constraints( bctypePot, gfsPot, ccPot); // assemble potential constraints
      Dune::PDELab::constraints( bctypeCon, gfsCon, ccCon); // assemble concentration constraints
      debug_info << "constrained potential dofs="     << ccPot.size() << " of " << gfsPot.globalSize() << std::endl;
      debug_info << "constrained concentration dofs=" << ccCon.size() << " of " << gfsCon.globalSize() << std::endl;

      // ###### Dirichlet values ########
      typedef DirichletValuesSingleCon<PARAMETERS_CON> DirichletValuesSingleCon;

      // Create NUMBER_OF_SPECIES Dirichlet value classes
      typedef Dune::PDELab::PowerGridFunction<DirichletValuesSingleCon,NUMBER_OF_SPECIES> DirichletValuesCon;
      DirichletValuesCon dirichletValuesCon;
      for(int i=0; i<NUMBER_OF_SPECIES; ++i)
      {
        DirichletValuesSingleCon* dirichletValuesSingleCon = new DirichletValuesSingleCon(gv,parametersCon,i);
        dirichletValuesCon.setChild(i,*dirichletValuesSingleCon);
      }
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PARAMETERS_POT> DirichletValuesPot;
      DirichletValuesPot dirichletValuesPot(gv,parametersPot);

      // Is this necessary? Initial conditions should fulfill the boundary conditions anyways,
      // otherwise we have an ill-posed problem!
      //Dune::PDELab::interpolate(dirichletValuesCon,gfsCon,unewCon);
      //Dune::PDELab::copy_nonconstrained_dofs(ccCon,uoldCon,unewCon);
      Dune::PDELab::interpolate(dirichletValuesPot,gfsPot,uPot);
      //Output::printSingleCoefficientVector(uPot, "initial pot");
      // ==============================================================================================

      // ============== Make grid operator space =========================================
      // ##### Spatial part #####
      // ### Standard FEM (CG) Nernst-Planck local operator
      typedef Dune::PDELab::ConvectionDiffusionFEM<PARAMETERS_CON,FEM_CON> LOP_SINGLE_CON;
      LOP_SINGLE_CON lopSingleCon(parametersCon);
      typedef NernstPlanckDGPowerOperator<PARAMETERS_CON,FEM_CON,LOP_SINGLE_CON> LOP_CON;
      LOP_CON lopCon(parametersCon,lopSingleCon);
	    
      // ### Standard FEM (CG) Poisson local operator
      typedef Dune::PDELab::ConvectionDiffusionFEM<PARAMETERS_POT,FEM_POT> LOP_POT;
      LOP_POT lopPot(parametersPot);
      
      // ##### Temporal part #####
      typedef NernstPlanckTimeLocalOperator<PHYSICS> TLOP_CON;
      TLOP_CON tlopCon(physics);

      typedef VBE::MatrixBackend MBE;
      
      typedef Dune::PDELab::GridOperator<GFS_POT,GFS_POT,LOP_POT,MBE,Real,Real,Real,CC_POT,CC_POT> GO_POT;
      GO_POT goPot(gfsPot,ccPot,gfsPot,ccPot,lopPot);
      typedef Dune::PDELab::GridOperator<GFS_CON,GFS_CON,LOP_CON,MBE,Real,Real,Real,CC_CON,CC_CON> GO_CON0;
      GO_CON0 goCon0(gfsCon,ccCon,gfsCon,ccCon,lopCon);
      typedef Dune::PDELab::GridOperator<GFS_CON,GFS_CON,TLOP_CON,MBE,Real,Real,Real,CC_CON,CC_CON> GO_CON1;
      GO_CON1 goCon1(gfsCon,ccCon,gfsCon,ccCon,tlopCon);
      typedef Dune::PDELab::OneStepGridOperator<GO_CON0,GO_CON1> GO_CON;
      GO_CON goCon(goCon0,goCon1);
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
      debug_verb << "============= Solving initial Poisson problem" << std::endl;
      solverPot.apply(uPot);
      debug_verb << "============= SOLVED initial Poisson problem" << std::endl;

      // =============================================================================================

      typedef Ax1StationaryLinearProblemSolver<GO_CON,LS_CON,U_CON> PDESOLVER_CON;
      PDESOLVER_CON pdesolverCon(goCon,lsCon,1e-10,"octave/nernst_planck");
      pdesolverCon.setPrintMatrix(false);
      pdesolverCon.setPrintRhs(false);
      pdesolverCon.setRowPreconditioner(useRowPreconditioner);
      pdesolverCon.setVerbosity(0); //2
      debug_verb << "Nernst-Planck PDE solver set up" << std::endl;
      // ==============================================================================================
    
    
      // ========== time-stepper ======================================================================
      Dune::PDELab::Alexander2Parameter<Real> timeStepperCon;
      //Dune::PDELab::ImplicitEulerParameter<Real> timeStepperCon;
      //Dune::PDELab::ExplicitEulerParameter<Real> timeStepperCon;
      //Dune::PDELab::RK4Parameter<Real> timeStepperCon;
      //Dune::PDELab::HeunParameter<Real> timeStepperCon;
      typedef Dune::PDELab::OneStepMethod<Real,GO_CON,PDESOLVER_CON,U_CON,U_CON> SOLVER_CON;
      SOLVER_CON solverCon(timeStepperCon,goCon,pdesolverCon);
      solverCon.setVerbosityLevel(0); // 2
      
      const PDESOLVER_CON& pdesolverConRef = solverCon.getPDESolver();
      // ==============================================================================================
      
      // ========== Initial output ====================================================================
      double debyeLength = physics.getDebyeLength(gfIonicStrength);
      debug_verb << "" << std::endl;
      debug_verb << "@@@@@@@ minimum Debye length / length scale " << debyeLength / physics.getLengthScale()<< std::endl;
      debug_verb << "" << std::endl;

      acme1MDOutput.writeStep(time);
      Tools::pecletNumber(gv, parametersCon);
      debug_info  << std::endl << "########## initial output done" << std::endl << std::endl;
      // ==============================================================================================

      // ========= time loop ==========================================================================
      // Setup solution method
      typedef Acme1MDOperatorSplit<U_CON,U_POT,DirichletValuesCon,DirichletValuesPot,SOLVER_CON,
      //typedef Acme1MDOperatorSplit<U_CON,U_POT,INITIAL_CON,INITIAL_POT,SOLVER_CON,
                          SOLVER_POT,GFS_CON,GFS_POT,PHYSICS,GV,GO_CON,GO_POT,Real> Acme1MDOperatorSplit;
      Acme1MDOperatorSplit acme1MDOperatorSplit(uoldCon, unewCon, uPot, dirichletValuesCon, dirichletValuesPot,
      //Acme1MDOperatorSplit acme1MDOperatorSplit(uoldCon, unewCon, uPot, initialCon, initialPot,
                         solverCon, solverPot, gfsCon, gfsPot, physics, gv, goCon, goPot);
      acme1MDOperatorSplit.setReduction(physics.getParams().getReduction());
      acme1MDOperatorSplit.setAbsLimit(physics.getParams().getAbsLimit());
      acme1MDOperatorSplit.setMaxIterations(100);
      acme1MDOperatorSplit.setUseDefect(false);

      Real dt = dtstart;
      Real outputTimeInterval = physics.getParams().getOutputTimeInterval();
      int outputCounter = 1;
      bool printEveryTimeStep = (dt > outputTimeInterval && !physics.getParams().useAdaptiveTimeStep());

      while (time<tend-1e-8)
      {
        //debug_verb << "TIME= " << time << std::endl;
        //debug_verb << "outputCounter= " << outputCounter << std::endl;

        // Hack to let PDE solver be able to get the current time
        acme1MDOutput.getDiagInfo().time = (time+dt);

        // !! Time-dependent constraints are currently not implemented in PDELab !!
        // evaluate constraints for current time step
//        bctypePot.setTime(time+dt);
//        bctypeCon.setTime(time+dt);
//        ccCon.clear();
//        Dune::PDELab::constraints( bctypeCon, gfsCon, ccCon );
//        ccPot.clear();
//        Dune::PDELab::constraints( bctypePot, gfsPot, ccPot );

        if(physics.getParams().useAdaptiveTimeStep())
        {
          Real safety_factor = 0.1;
          dt = safety_factor * Tools::getTimeStep(gv, physics, dgfGradPot);
          debug_verb << "Calculated time step: " << dt << std::endl;
        }
        physics.setTimeStep(dt);
        acme1MDOutput.getDiagInfo().dt = dt;
        
        acme1MDOperatorSplit.timeStep(time,dt,acme1MDOutput.getDiagInfo());


        U_CON uChange(unewCon);
        uChange -= uoldCon;

        //Output::printSingleCoefficientVector(uChange, "[Na] update");
        debug_info << "con change = " << uChange.base().two_norm() << std::endl;

        uoldCon = unewCon;
        
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

        	acme1MDOutput.writeStep(time);
        	Tools::pecletNumber(gv, parametersCon);
        	acme1MDOutput.getDiagInfo().clear();

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

#endif /* DUNE_AX1_ACME1MD_SETUP_NOMEMBRANE_HH */
