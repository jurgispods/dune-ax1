#ifndef DUNE_AX1_ACME0_PK_FULLY_IMPLICIT_HH
#define DUNE_AX1_ACME0_PK_FULLY_IMPLICIT_HH

#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/finiteelementmap/p1fem.hh>
#include <dune/pdelab/finiteelementmap/pk1dbasis.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/finiteelementmap/opbfem.hh>
#include <dune/pdelab/function/selectcomponent.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/backend/seqistlsolverbackend.hh>
#include <dune/pdelab/stationary/linearproblem.hh>


#include <dune/ax1/common/vectordiscretegridfunctiongradient.hh>
#include <dune/ax1/acme0/fully-implicit/acme0_operator_fully_implicit.hh>
#include <dune/ax1/acme0/fully-implicit/acme0_toperator_fully_implicit.hh>
#include <dune/ax1/acme0/common/acme0_output.hh>
#include <dune/ax1/acme0/common/acme0_initial.hh>
#include <dune/ax1/acme0/common/acme0_boundary.hh>
#include <dune/ax1/acme0/common/nernst_planck_parameters.hh>
#include <dune/ax1/acme0/common/poisson_parameters.hh>
#include <dune/ax1/acme0/common/convectiondiffusiondg.hh>
#include <dune/ax1/common/ax1_linearproblem.hh>
#include <dune/ax1/common/ax1_newton.hh>
#include <dune/ax1/common/output.hh>


template<class GV, class PHYSICS>
class Acme0PkFullyImplicit
{
  public:
    
    // Choose domain and range field type
    typedef GV GridView;
    typedef typename GV::Grid::ctype Coord;
    typedef double Real;
    
    
    Acme0PkFullyImplicit(const GV& gv_, PHYSICS& physics_) :
      gv(gv_),
      physics(physics_),
      tolCon(physics.getParams().getToleranceCon()),
      tolPot(physics.getParams().getTolerancePot())
    {}
    
    void acme0_Pk (double dtstart, double tend)
    {
      const int dim = GV::dimension;
      Real time = 0.0;
      const int degree=1;

      typedef Dune::PDELab::ConformingDirichletConstraints CONSTRAINTS; // constraints class (Dirichlet)
      typedef Dune::PDELab::ISTLVectorBackend<1> VBE;                   // block size = numOfEquations

      // =================== GFS SETUP ================================================================
      // Create grid function space for concentration of one single ion species
      //typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM_CON;
      //FEM_CON femCon(1);
      typedef Dune::PDELab::OPBLocalFiniteElementMap<Coord,Real,degree,dim,Dune::GeometryType::simplex> FEM_CON_DG;
      FEM_CON_DG femConDG;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM_CON_DG,CONSTRAINTS,VBE> GFS_SINGLE_CON;
      GFS_SINGLE_CON gfsSingleCon(gv,femConDG);

      // Power grid function space for all ion species
      typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper PowerGFSOrdering;
      //typedef Dune::PDELab::GridFunctionSpaceBlockwiseMapper PowerGFSOrdering;
      //typedef Dune::PDELab::GridFunctionSpaceComponentBlockwiseMapper<1> PowerGFSOrdering;

      typedef Dune::PDELab::PowerGridFunctionSpace<
      GFS_SINGLE_CON, NUMBER_OF_SPECIES, PowerGFSOrdering> GFS_CON;
      GFS_CON gfsCon(gfsSingleCon);

      // Create grid function space for potential
      //typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM_POT;
      //FEM_POT femPot(1);
      typedef Dune::PDELab::OPBLocalFiniteElementMap<Coord,Real,degree,dim,Dune::GeometryType::simplex> FEM_POT_DG;
      FEM_POT_DG femPotDG;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM_POT_DG,CONSTRAINTS,VBE> GFS_POT;
      GFS_POT gfsPot(gv,femPotDG);

      // Now wrap grid function spaces into composite GFS
      //typedef Dune::PDELab::GridFunctionSpaceComponentBlockwiseMapper<3,1> CompositeGFSOrdering;
      //typedef Dune::PDELab::GridFunctionSpaceBlockwiseMapper CompositeGFSOrdering;
      typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper CompositeGFSOrdering;

      typedef Dune::PDELab::CompositeGridFunctionSpace<
          CompositeGFSOrdering,GFS_CON,GFS_POT> GFS;
      GFS gfs(gfsCon,gfsPot);

      // Extract solution vector types
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,Real>::Type U;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CON,Real>::Type U_CON;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_POT,Real>::Type U_POT;

      // Types for different grid functions
      typedef Dune::PDELab::VectorDiscreteGridFunction<GFS_CON,U_CON> DGF_CON;
      typedef Dune::PDELab::DiscreteGridFunction<GFS_POT,U_POT> DGF_POT;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_POT,U_POT> DGF_POT_GRAD;
      typedef ChargeDensityGridFunction<DGF_CON,PHYSICS> GF_CD;
      // ==============================================================================================

      // =================== BOUNDARY CONDITIONS=======================================================
      typedef BCTypeSingleCon BCT_SINGLE_CON;
      typedef Dune::PDELab::PowerConstraintsParameters<BCTypeSingleCon,NUMBER_OF_SPECIES> BCT_CON;
      typedef BCTypePot BCT_POT;
      typedef Dune::PDELab::CompositeConstraintsParameters<BCT_CON,BCT_POT> BCT;

      BCT_SINGLE_CON bctSingleCon;
      BCT_CON bctCon(bctSingleCon);
      BCT_POT bctPot;
      BCT bct(bctCon,bctPot);

      typedef typename GFS_POT::template ConstraintsContainer<Real>::Type CC_POT;
      CC_POT ccPot;
      Dune::PDELab::constraints( bctPot, gfsPot, ccPot ); // assemble constraints

      typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
      CC cc;
      Dune::PDELab::constraints( bct, gfs, cc ); // assemble constraints
      debug_info << "constrained dofs=" << cc.size() << " of " << gfs.globalSize() << std::endl;
      // ==============================================================================================

      // =========== Define solution vector types =====================================================
      U uold(gfs,0.0);
      U unew(gfs,0.0);
      U_CON uCon(gfsCon,0.0);
      U_POT uPot(gfsPot,0.0);
      // ==============================================================================================

      // =========== Make FE function with initial value ==============================================
      //typedef Dune::PDELab::CompositeGridFunction<INITIAL_CON,INITIAL_POT> INITIAL_U;

      // Inital concentrations
      //typedef InitialCon<GV,Real,NUMBER_OF_SPECIES> INITIAL_CON;
      typedef typename PHYSICS::Traits::INITIAL_CON INITIAL_CON;
      INITIAL_CON initialCon(gv,physics.getParams());
      Dune::PDELab::interpolate(initialCon,gfsCon,uCon);

      DGF_CON dgfCon(gfsCon, uCon);
      GF_CD   gfChargeDensity(dgfCon, physics);

      // Initial potential
      typedef InitialPot<GV,Real,PHYSICS,GF_CD> INITIAL_POT;
      INITIAL_POT initialPot(gv,physics,gfChargeDensity);
      initialPot.setTime(time);
      Dune::PDELab::interpolate(initialPot,gfsPot,uPot);

      DGF_POT dgfPot(gfsPot, uPot);
      DGF_POT_GRAD dgfGradPot(gfsPot, uPot);

      typedef Dune::PDELab::CompositeGridFunction<INITIAL_CON,INITIAL_POT> INITIAL_U;
      INITIAL_U initialU(initialCon,initialPot);
      Dune::PDELab::interpolate(initialU,gfs,uold);
      // ==============================================================================================


      // ============== Make grid operator  ===========================================================
      typedef NernstPlanckParameters<GV,Real,PHYSICS,DGF_POT_GRAD> PARAMETERS_CON;
      PARAMETERS_CON parametersCon(physics,dgfGradPot);
      Dune::PDELab::ConvectionDiffusionDGMethod::Type mCon;
      mCon = Dune::PDELab::ConvectionDiffusionDGMethod::SIPG;
      Dune::PDELab::ConvectionDiffusionDGWeights::Type wCon;
      wCon = Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn;
      Real alphaCon = 2.0;

      typedef PoissonParameters<GV,Real,PHYSICS,GF_CD,INITIAL_POT> PARAMETERS_POT;
      PARAMETERS_POT parametersPot(physics,gfChargeDensity,initialPot);
      Dune::PDELab::ConvectionDiffusionDGMethod::Type mPot;
      mPot = Dune::PDELab::ConvectionDiffusionDGMethod::SIPG;
      Dune::PDELab::ConvectionDiffusionDGWeights::Type wPot;
      wPot = Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn;
      Real alphaPot = 2.0;

      typedef Acme0OperatorFullyImplicit<PARAMETERS_CON,PARAMETERS_POT,FEM_CON_DG,FEM_POT_DG> LOP;
      LOP lop(parametersCon, parametersPot, mCon, mPot, wCon, wPot, alphaCon, alphaPot);    // spatial part
      typedef Acme0TimeLocalOperator TLOP;
      TLOP tlop(2);                          // temporal part
      typedef VBE::MatrixBackend MBE;
      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO0;
      GO0 go0(gfs,cc,gfs,cc,lop);
      typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,Real,Real,Real,CC,CC> GO1;
      GO1 go1(gfs,cc,gfs,cc,tlop);
      typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
      IGO igo(go0,go1);

      // ### Poisson operator for initial Poisson solve
      typedef Dune::PDELab::ConvectionDiffusionDG<PARAMETERS_POT,FEM_POT_DG> LOP_POT_DG;
      LOP_POT_DG lopPotDG(parametersPot,mPot,wPot,alphaPot); // alpha=2.0 [wtf is this?]
      typedef Dune::PDELab::GridOperator<GFS_POT,GFS_POT,LOP_POT_DG,MBE,Real,Real,Real,CC_POT,CC_POT> GO_POT;
      GO_POT goPot(gfsPot,ccPot,gfsPot,ccPot,lopPotDG);
      debug_verb << "Poisson grid operator set up" << std::endl;
      // ==============================================================================================


      // ========== Select a linear solver backend ====================================================
      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
      //LS ls(5000,false);
      typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
      LS ls(0); // verbose = 1
      // ==============================================================================================

      // ========= Solver for non-linear problem per stage ============================================
      // Poisson solver
      typedef Ax1StationaryLinearProblemSolver<GO_POT,LS,U_POT> SLP_POT;
      SLP_POT solverPot(goPot,ls,1e-10,"octave/fully-implicit-poisson_mat");
      solverPot.setPrintMatrix(false);
      solverPot.setPrintRhs(false);
      //solverPot.setRowPreconditioner(true);
      debug_verb << "============= Solving initial Poisson problem" << std::endl;
      solverPot.apply(uPot);
      debug_verb << "============= SOLVED initial Poisson problem" << std::endl;

      typedef Ax1Newton<IGO,LS,U> PDESOLVER;
      PDESOLVER pdesolver(igo,ls,"octave/fully-implicit_mat");
      pdesolver.setReassembleThreshold(0.0);
      pdesolver.setVerbosityLevel(0); // 2
      pdesolver.setReduction(std::min(tolCon,tolPot)); // 1e-10
      pdesolver.setMinLinearReduction(1e-4);
      pdesolver.setMaxIterations(25);
      pdesolver.setLineSearchMaxIterations(10); // 10
      pdesolver.setPrintMatrix(false);
      pdesolver.setPrintRhs(false);
      //pdesolver.setRowPreconditioner(true);

      /*typedef Dune::PDELab::CompositeGridFunction<INITIAL_CON,DGF_POT> INITIAL_U;
      INITIAL_U initialU(initialCon,dgfPot);
      Dune::PDELab::interpolate(initialU,gfs,uold);*/

      Tools::childToCompositeCoefficientVector(gfs, uCon, 0, uold);
      Tools::childToCompositeCoefficientVector(gfs, uPot, 1, uold);
      unew = uold;
      // ==============================================================================================


      // ========== time-stepper ======================================================================
      //Dune::PDELab::Alexander2Parameter<Real> method;
      Dune::PDELab::ImplicitEulerParameter<Real> timeStepper;
      //Dune::PDELab::RK4Parameter<Real> timeStepper;
      Dune::PDELab::OneStepMethod<Real,IGO,PDESOLVER,U,U> solver(timeStepper,igo,pdesolver);
      solver.setVerbosityLevel(0); // 2
      // ==============================================================================================

      // ========== graphics for initial output =======================================================
      typedef Acme0Output<GFS_CON,GFS_POT,U_CON,U_POT,PHYSICS> Acme0Output;
      Acme0Output acme0Output(gfsCon, gfsPot, uCon, uPot, physics);
      typename Acme0Output::DiagnosticInfo diagInfo;
      Tools::compositeToChildCoefficientVector(gfs, unew, uCon, 0);
      Tools::compositeToChildCoefficientVector(gfs, unew, uPot, 1);
      Output::printSingleCoefficientVectorDG(uPot, "pot");
      Output::printMultipleComponentCoefficientVectorDG(uCon, NUMBER_OF_SPECIES);
      acme0Output.writeStep(time, diagInfo);
      Tools::pecletNumber(gv, parametersCon);
      debug_info  << std::endl << "########## initial output done" << std::endl << std::endl;
      // ==============================================================================================


      // ========= time loop ==========================================================================
      Real dt = dtstart;
      Real outputTimeInterval = physics.getParams().getOutputTimeInterval();
      int outputCounter = 1;
      
      while (time<tend-1e-8)
      {
        // if constraints depend on time or solution
        bctPot.setTime(time+dt);
        // evaluate constraints for current time step
        cc.clear();
        Dune::PDELab::constraints(bct,gfs,cc);

        //const PDESOLVER& newton = osm.getPDESolver();
        if(physics.getParams().useAdaptiveTimeStep())
        {
          Real safety_factor = 0.1;
          dt = safety_factor * Tools::getTimeStep(gv, physics, dgfGradPot);
          debug_verb << "Calculated time step: " << dt << std::endl;
        }
        physics.setTimeStep(dt);

        //debug_verb << "================= Solving fully implicit system..." << std::endl;
        solver.apply(time,dt,uold,initialU,unew);
        // The following two lines are essential and must NOT be removed!
        // Grid functions (e.g. charge density) rely on the current value of child GFS coefficient vectors
        Tools::compositeToChildCoefficientVector(gfs, unew, uCon, 0);
        Tools::compositeToChildCoefficientVector(gfs, unew, uPot, 1);
        //debug_verb << "================= Solved fully implicit system" << std::endl;

        // Permanently use new vector of gating particles
        physics.getMembrane().updateState();

        uold = unew;
        time += dt;

        if (std::abs(time - outputCounter * outputTimeInterval) < 1e-8)
        {
        	++outputCounter;
        	
	        //Output::printSingleCoefficientVectorDG(uPot, "pot");
          //Output::printMultipleComponentCoefficientVectorDG(uCon, NUMBER_OF_SPECIES);
    	    
    	    // Fill diagnosticInfo
    	    diagInfo.iterations = solver.getPDESolver().result().iterations;
          diagInfo.dt = dt;
    	    
      	  acme0Output.writeStep(time, diagInfo);
        	Tools::pecletNumber(gv, parametersCon);

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
     const Real tolCon;
     const Real tolPot;
};

#endif /* DUNE_AX1_ACME0_PK_FULLY_IMPLICIT_HH */
