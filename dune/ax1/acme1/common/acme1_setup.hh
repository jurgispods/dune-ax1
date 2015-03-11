#ifndef DUNE_AX1_ACME1_SETUP_HH
#define DUNE_AX1_ACME1_SETUP_HH

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
#if USE_SUBGRID==2
#include <dune/ax1/common/ax1_subgrid_tools.hh>
#endif
#if USE_SUBGRID==1
#include <dune/ax1/common/ax1_griddata_transfer.hh>
#endif
#include <dune/ax1/common/ax1_discretegridfunction.hh>
#include <dune/ax1/common/ax1_simulationdata.hh>
#include <dune/ax1/common/ax1_simulationstate.hh>
#include <dune/ax1/common/chargedensitygridfunction.hh>
#include <dune/ax1/common/ionicstrengthgridfunction.hh>
#include <dune/ax1/common/poisson_boltzmann_concentrationgridfunction.hh>
#include <dune/ax1/common/poisson_boltzmann_rhs_gridfunction.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/vectordiscretegridfunctiongradient.hh>
#include <dune/ax1/acme1/common/acme1_boundary.hh>
#include <dune/ax1/acme1/common/acme1_initial.hh>
#include <dune/ax1/acme1/common/acme1_output.hh>
#include <dune/ax1/acme1/common/acme1_solutionvectors.hh>
#include <dune/ax1/acme1/common/convectiondiffusionfem.hh>
#include <dune/ax1/acme1/common/convectiondiffusiondg.hh>
#include <dune/ax1/acme1/common/nernst_planck_parameters.hh>
#include <dune/ax1/acme1/common/poisson_parameters.hh>
#include <dune/ax1/acme1/common/poisson_boltzmann_parameters.hh>
#include <dune/ax1/acme1/common/poisson_boltzmann_operator.hh>
#include <dune/ax1/acme1/operator-split/acme1_operator_split.hh>
#include <dune/ax1/acme1/operator-split/acme1_iteration.hh>
#include <dune/ax1/acme1/operator-split/acme1_toperator.hh>
#include <dune/ax1/acme1/operator-split/nernst_planck_power_operator.hh>

#define USE_CG 1

template<class GV, class PHYSICS, class SubGV>
class Acme1Setup
{
  public:
    
    //! Traits class providing function spaces and related data types for the acme1 use case
    template<int NUMBER_OF_SPECIES>
    struct Acme1Traits
    {
      typedef Acme1Setup<GV,PHYSICS,SubGV> SETUP;

      typedef PHYSICS Physics;

      typedef GV GridView;
      typedef typename GV::Grid GridType;
      typedef typename GV::Grid::ctype Coord;
      typedef double Real;

      typedef SubGV SubGridView;
      typedef typename SubGV::Grid SubGridType;

      typedef Dune::PDELab::ConformingDirichletConstraints CONSTRAINTS;
      typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
      typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM_CON;
      typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM_POT;

      typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper PowerGFSOrdering;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM_CON,CONSTRAINTS,VBE> GFS_SINGLE_CON;
      typedef Dune::PDELab::PowerGridFunctionSpace<
        GFS_SINGLE_CON, NUMBER_OF_SPECIES, PowerGFSOrdering> GFS_CON;

      typedef Dune::PDELab::GridFunctionSpace<SubGV,FEM_CON,CONSTRAINTS,VBE> SUB_GFS_SINGLE_CON;
      typedef Dune::PDELab::PowerGridFunctionSpace<
        SUB_GFS_SINGLE_CON, NUMBER_OF_SPECIES, PowerGFSOrdering> SUB_GFS_CON;

      typedef Dune::PDELab::GridFunctionSpace<GV,FEM_POT,CONSTRAINTS,VBE> GFS_POT;
      typedef Dune::PDELab::GridFunctionSpace<SubGV,FEM_POT,CONSTRAINTS,VBE> SUB_GFS_POT;

      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CON,Real>::Type U_CON;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_POT,Real>::Type U_POT;

      typedef typename Dune::PDELab::BackendVectorSelector<SUB_GFS_CON,Real>::Type USUB_CON;
      typedef typename Dune::PDELab::BackendVectorSelector<SUB_GFS_POT,Real>::Type USUB_POT;

      typedef Acme1Output<GFS_CON,GFS_POT,U_CON,U_POT,PHYSICS> ACME1_OUTPUT;

#if USE_SUBGRID==2
      typedef SubGridCoefficientVectorRestrictor<GV,SubGV,GFS_POT,SUB_GFS_POT,U_POT,USUB_POT> UPOT_Restrictor;
      typedef SubGridCoefficientVectorRestrictor<GV,SubGV,GFS_CON,SUB_GFS_CON,U_CON,USUB_CON> UCON_Restrictor;
      typedef SubGridCoefficientVectorInterpolator<SubGV,GV,SUB_GFS_POT,GFS_POT,USUB_POT,U_POT> UPOT_Interpolator;
      typedef SubGridCoefficientVectorInterpolator<SubGV,GV,SUB_GFS_CON,GFS_CON,USUB_CON,U_CON> UCON_Interpolator;
#endif
#if USE_SUBGRID==1
      typedef GridCoefficientVectorRestrictor<GV,SubGV,GFS_POT,SUB_GFS_POT,U_POT,USUB_POT> UPOT_Restrictor;
      typedef GridCoefficientVectorRestrictor<GV,SubGV,GFS_CON,SUB_GFS_CON,U_CON,USUB_CON> UCON_Restrictor;
      typedef GridCoefficientVectorInterpolator<SubGV,GV,SUB_GFS_POT,GFS_POT,USUB_POT,U_POT> UPOT_Interpolator;
      typedef GridCoefficientVectorInterpolator<SubGV,GV,SUB_GFS_CON,GFS_CON,USUB_CON,U_CON> UCON_Interpolator;
#endif
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
    typedef Acme1Traits<NUMBER_OF_SPECIES> Traits;

    static const bool useRowPreconditioner = false;
    static const bool writeIfNotConverged = true;
    static const bool resetTime = true;

    Acme1Setup(const GV& gv_, PHYSICS& physics_,
        const SubGV& subGridView_Inside_, const SubGV& subGridView_Outside_)
    : gv(gv_),
      physics(physics_),
      subGridView_Inside(subGridView_Inside_),
      subGridView_Outside(subGridView_Outside_)
    {}

    void setup (double dtstart, double tend)
    {
      Real time = 0.0;
      Real dt = dtstart;
      Real tEquilibrium = 1000.0;

      const int dim = GV::dimension;
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
      // Full GFS on host grid
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM_CON,CONSTRAINTS,VBE> GFS_SINGLE_CON;
      GFS_SINGLE_CON gfsSingleCon(gv,femCon);
      // Power grid function space for all ion species
      typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper PowerGFSOrdering;
      typedef Dune::PDELab::PowerGridFunctionSpace<
        GFS_SINGLE_CON, NUMBER_OF_SPECIES, PowerGFSOrdering> GFS_CON;
      GFS_CON gfsCon(gfsSingleCon);

      // Actual GFS for calculations on subgrid
      typedef Dune::PDELab::GridFunctionSpace<SubGV,FEM_CON,CONSTRAINTS,VBE> SUB_GFS_SINGLE_CON;
      SUB_GFS_SINGLE_CON subGfsSingleCon_Inside(subGridView_Inside,femCon);
      SUB_GFS_SINGLE_CON subGfsSingleCon_Outside(subGridView_Outside,femCon);
      // Power grid function space for all ion species
      typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper PowerGFSOrdering;
      typedef Dune::PDELab::PowerGridFunctionSpace<
        SUB_GFS_SINGLE_CON, NUMBER_OF_SPECIES, PowerGFSOrdering> SUB_GFS_CON;
      SUB_GFS_CON subGfsCon_Inside(subGfsSingleCon_Inside);
      SUB_GFS_CON subGfsCon_Outside(subGfsSingleCon_Outside);

      // Grid function space for potential on host grid
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM_POT,CONSTRAINTS,VBE> GFS_POT;
      GFS_POT gfsPot(gv,femPot);

      // Grid function space for potential on subgrid
      typedef Dune::PDELab::GridFunctionSpace<SubGV,FEM_POT,CONSTRAINTS,VBE> SUB_GFS_POT;
      SUB_GFS_POT subGfsPot_Inside(subGridView_Inside,femPot);
      SUB_GFS_POT subGfsPot_Outside(subGridView_Outside,femPot);

      // Extract solution vector types
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CON,Real>::Type U_CON;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_POT,Real>::Type U_POT;

      typedef typename Dune::PDELab::BackendVectorSelector<SUB_GFS_CON,Real>::Type USUB_CON;
      typedef typename Dune::PDELab::BackendVectorSelector<SUB_GFS_POT,Real>::Type USUB_POT;
      
      // Types for different grid functions
      typedef Ax1VectorDiscreteGridFunction <PHYSICS,GFS_CON,U_CON> DGF_CON;
      typedef Dune::PDELab::DiscreteGridFunction<GFS_POT,U_POT> DGF_POT;

      // Helper grid functions for concentrations
      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_POT,U_POT> DGF_POT_GRAD;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<SUB_GFS_POT,USUB_POT> DGF_SUB_POT_GRAD;

      // Helper grid function for potential
      typedef ChargeDensityGridFunction<DGF_CON, PHYSICS> GF_CD;
      typedef IonicStrengthGridFunction<DGF_CON, PHYSICS> GF_IS;

      // Grid functions for initial values
      typedef typename PHYSICS::Traits::INITIAL_CON INITIAL_GF;
      typedef InitialCon<PHYSICS,INITIAL_GF,GV,Real,NUMBER_OF_SPECIES> INITIAL_CON;
      typedef InitialPot<GV,Real,PHYSICS,GF_CD> INITIAL_POT;

      // Helper grid function to calculate equilibirum concentrations from stationary Poisson-Boltzmann potential
      typedef PoissonBoltzmannRHSGridFunction<DGF_POT, INITIAL_CON, PHYSICS> GF_PB_RHS;
      typedef PoissonBoltzmannConcentrationGridFunction<INITIAL_CON, DGF_POT, SubGV, PHYSICS>
        GF_PB_CON;
      // ==============================================================================================

      // =========== Define solution vectors ==========================================================
      // Coefficient vectors on subgrids
      USUB_CON uoldCon_Inside(subGfsCon_Inside,0.0);
      USUB_CON unewCon_Inside(subGfsCon_Inside,0.0);
      USUB_POT uPot_Inside(subGfsPot_Inside,0.0);

      USUB_CON uoldCon_Outside(subGfsCon_Outside,0.0);
      USUB_CON unewCon_Outside(subGfsCon_Outside,0.0);
      USUB_POT uPot_Outside(subGfsPot_Outside,0.0);

      // Coefficient vectors on host grid
      U_CON uCon(gfsCon,0.0);
      U_CON uConPrevious(gfsCon,0.0);
      U_POT uPot(gfsPot,0.0);

      // Now let's put 'em all into one structure!
      Acme1SolutionVectors<Traits> solutionVectors(uCon, uPot, uoldCon_Inside, unewCon_Inside, uPot_Inside,
          uoldCon_Outside, unewCon_Outside, uPot_Outside);


      // Classes for transfer of coefficient vectors between host grid and subgrid
#if USE_SUBGRID==2
      typedef SubGridCoefficientVectorRestrictor<GV,SubGV,GFS_POT,SUB_GFS_POT,U_POT,USUB_POT> UPOT_Restrictor;
      typedef SubGridCoefficientVectorRestrictor<GV,SubGV,GFS_CON,SUB_GFS_CON,U_CON,USUB_CON> UCON_Restrictor;
      typedef SubGridCoefficientVectorInterpolator<SubGV,GV,SUB_GFS_POT,GFS_POT,USUB_POT,U_POT> UPOT_Interpolator;
      typedef SubGridCoefficientVectorInterpolator<SubGV,GV,SUB_GFS_CON,GFS_CON,USUB_CON,U_CON> UCON_Interpolator;
#endif
#if USE_SUBGRID==1
      typedef GridCoefficientVectorRestrictor<GV,SubGV,GFS_POT,SUB_GFS_POT,U_POT,USUB_POT> UPOT_Restrictor;
      typedef GridCoefficientVectorRestrictor<GV,SubGV,GFS_CON,SUB_GFS_CON,U_CON,USUB_CON> UCON_Restrictor;
      typedef GridCoefficientVectorInterpolator<SubGV,GV,SUB_GFS_POT,GFS_POT,USUB_POT,U_POT> UPOT_Interpolator;
      typedef GridCoefficientVectorInterpolator<SubGV,GV,SUB_GFS_CON,GFS_CON,USUB_CON,U_CON> UCON_Interpolator;
#endif
      // ==============================================================================================

      // ====== Create grid functions and interpolate initial values onto coefficient vectors =========
      // Initial concentrations
      INITIAL_GF initialGF(gv,physics.getParams());

      // Generic wrapper class for initial grid function
      INITIAL_CON initialCon(physics,initialGF,gv,physics.getParams());
      Dune::PDELab::interpolate(initialCon,gfsCon,uCon);
      uConPrevious = uCon;

      DGF_CON dgfCon(physics,gfsCon, uCon);
      DGF_CON dgfConPrevious(physics,gfsCon, uConPrevious);
      GF_CD gfChargeDensity(dgfCon, dgfConPrevious, physics);
      GF_IS gfIonicStrength(dgfCon, physics);

      // Initial potential
      INITIAL_POT initialPot(gv,physics,gfChargeDensity,2*degreePot);
      initialPot.setTime(time);
      Dune::PDELab::interpolate(initialPot,gfsPot,uPot);
      //Output::printSingleCoefficientVector(uPot, "uPot initial");

      DGF_POT dgfPot(gfsPot, uPot);
      
      UPOT_Restrictor::restrict(gv, subGridView_Inside, gfsPot, subGfsPot_Inside, uPot, uPot_Inside);
      UPOT_Restrictor::restrict(gv, subGridView_Outside, gfsPot, subGfsPot_Outside, uPot, uPot_Outside);

      UCON_Restrictor::restrict(gv, subGridView_Inside, gfsCon, subGfsCon_Inside, uCon, uoldCon_Inside);
      UCON_Restrictor::restrict(gv, subGridView_Outside, gfsCon, subGfsCon_Outside, uCon, uoldCon_Outside);
      unewCon_Inside = uoldCon_Inside;
      unewCon_Outside = uoldCon_Outside;

      DGF_POT_GRAD dgfGradPot(gfsPot, uPot);
      DGF_SUB_POT_GRAD dgfGradPot_Inside(subGfsPot_Inside, uPot_Inside);
      DGF_SUB_POT_GRAD dgfGradPot_Outside(subGfsPot_Outside, uPot_Outside);

      // Define and instantiate output class
      typedef Acme1Output<GFS_CON,GFS_POT,U_CON,U_POT,PHYSICS> Acme1Output;
      Acme1Output acme1Output(gfsCon, gfsPot, uCon, uPot, physics, 2*degreePot, 2*degreeCon);
      // ==============================================================================================

      // ========== Define Parameter classes containing the model problem =============================
      typedef typename PHYSICS::Traits::NERNST_PLANCK_BOUNDARY BOUNDARY_CON;
      BOUNDARY_CON boundaryCon(initialGF);
      typedef NernstPlanckParameters<SubGV,Real,PHYSICS,DGF_SUB_POT_GRAD,BOUNDARY_CON> PARAMETERS_CON;
      PARAMETERS_CON parametersCon_Inside(subGridView_Inside,physics,dgfGradPot_Inside,boundaryCon,tEquilibrium);
      PARAMETERS_CON parametersCon_Outside(subGridView_Outside,physics,dgfGradPot_Outside,boundaryCon,tEquilibrium);

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
      BCType_CON bctypeCon_Inside;
      BCType_CON bctypeCon_Outside;
      for(int i=0; i<NUMBER_OF_SPECIES; ++i)
      {
        BCType_SINGLE_CON* bctypeSingleCon_Inside = new BCType_SINGLE_CON(subGridView_Inside,parametersCon_Inside,i);
        bctypeCon_Inside.setChild(i,*bctypeSingleCon_Inside);
        BCType_SINGLE_CON* bctypeSingleCon_Outside = new BCType_SINGLE_CON(subGridView_Outside,parametersCon_Outside,i);
        bctypeCon_Outside.setChild(i,*bctypeSingleCon_Outside);
      }
      BCType_POT bctypePot(gv,parametersPot);

      typedef typename SUB_GFS_CON::template ConstraintsContainer<Real>::Type CC_CON;
      CC_CON ccCon_Inside;
      CC_CON ccCon_Outside;
      typedef typename GFS_POT::template ConstraintsContainer<Real>::Type CC_POT;
      CC_POT ccPot;

      Dune::PDELab::constraints( bctypePot, gfsPot, ccPot); // assemble potential constraints
      Dune::PDELab::constraints( bctypeCon_Inside, subGfsCon_Inside, ccCon_Inside); // assemble concentration constraints
      Dune::PDELab::constraints( bctypeCon_Outside, subGfsCon_Outside, ccCon_Outside); // assemble concentration constraints
      debug_info << "constrained potential dofs="     << ccPot.size() << " of " << gfsPot.globalSize() << std::endl;
      debug_info << "constrained cytosol concentration dofs=" << ccCon_Inside.size() << " of "
          << subGfsCon_Inside.globalSize() << std::endl;
      debug_info << "constrained extracellular concentration dofs=" << ccCon_Outside.size() << " of "
          << subGfsCon_Outside.globalSize() << std::endl;


      // ###### Dirichlet values ########
      typedef DirichletValuesSingleCon<PARAMETERS_CON> DirichletValuesSingleCon;

      // Create NUMBER_OF_SPECIES Dirichlet value classes
      typedef Dune::PDELab::PowerGridFunction<DirichletValuesSingleCon,NUMBER_OF_SPECIES> DirichletValuesCon;
      DirichletValuesCon dirichletValuesCon_Inside;
      DirichletValuesCon dirichletValuesCon_Outside;
      for(int i=0; i<NUMBER_OF_SPECIES; ++i)
      {
        DirichletValuesSingleCon* dirichletValuesSingleCon_Inside
          = new DirichletValuesSingleCon(subGridView_Inside,parametersCon_Inside,i);
        DirichletValuesSingleCon* dirichletValuesSingleCon_Outside
            = new DirichletValuesSingleCon(subGridView_Outside,parametersCon_Outside,i);
        dirichletValuesCon_Inside.setChild(i,*dirichletValuesSingleCon_Inside);
        dirichletValuesCon_Outside.setChild(i,*dirichletValuesSingleCon_Outside);
      }
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PARAMETERS_POT> DirichletValuesPot;
      DirichletValuesPot dirichletValuesPot(gv,parametersPot);

      // Is this necessary? Initial conditions should fulfill the boundary conditions anyways,
      // otherwise we have an ill-posed problem!
      //Dune::PDELab::interpolate(dirichletValuesCon_Inside,subGfsCon_Inside,unewCon_Inside);
      //Dune::PDELab::copy_nonconstrained_dofs(ccCon_Inside,uoldCon_Inside,unewCon_Inside);
      //Dune::PDELab::interpolate(dirichletValuesCon_Outside,subGfsCon_Outside,unewCon_Outside);
      //Dune::PDELab::copy_nonconstrained_dofs(ccCon_Outside,uoldCon_Outside,unewCon_Outside);
      Dune::PDELab::interpolate(dirichletValuesPot,gfsPot,uPot);
      //Output::printSingleCoefficientVector(uPot, "uPot t=0");
      // ==============================================================================================

      // ============== Make grid operator space =========================================
      // ##### Stationary operator #####
      //typedef PoissonBoltzmannLocalOperator<PHYSICS,PARAMETERS_POT,INITIAL_CON> STATIONARY_LOP_POT;
      //STATIONARY_LOP_POT stationaryLopPot(physics, parametersPot, initialCon);

      GF_PB_RHS gfPoissonBoltzmannRHS(dgfPot, initialCon, physics);
      typedef PoissonBoltzmannParameters<GV,Real,PHYSICS,GF_PB_RHS,BOUNDARY_POT> PARAMETERS_PB;
      PARAMETERS_PB parametersPB(physics,gfPoissonBoltzmannRHS,boundaryPot);

      typedef Dune::PDELab::ConvectionDiffusionFEM<PARAMETERS_PB,FEM_POT> STATIONARY_LOP_POT;
      STATIONARY_LOP_POT stationaryLopPot(parametersPB);

      // ##### Spatial part #####
      // ### Standard FEM (CG) Nernst-Planck local operator
      typedef Dune::PDELab::ConvectionDiffusionFEM<PARAMETERS_CON,FEM_CON> LOP_SINGLE_CON;
      LOP_SINGLE_CON lopSingleCon_Inside(parametersCon_Inside);
      LOP_SINGLE_CON lopSingleCon_Outside(parametersCon_Outside);
      typedef NernstPlanckDGPowerOperator<PARAMETERS_CON,FEM_CON,LOP_SINGLE_CON> LOP_CON;
      LOP_CON lopCon_Inside(parametersCon_Inside,lopSingleCon_Inside);
      LOP_CON lopCon_Outside(parametersCon_Outside,lopSingleCon_Outside);
	    
      // ### Standard FEM (CG) Poisson local operator
      typedef Dune::PDELab::ConvectionDiffusionFEM<PARAMETERS_POT,FEM_POT> LOP_POT;
      LOP_POT lopPot(parametersPot);
      
      // ##### Temporal part #####
      typedef NernstPlanckTimeLocalOperator<PHYSICS> TLOP_CON;
      TLOP_CON tlopCon_Inside(physics);
      TLOP_CON tlopCon_Outside(physics);

      typedef VBE::MatrixBackend MBE;
      typedef Dune::PDELab::GridOperator<GFS_POT,GFS_POT,STATIONARY_LOP_POT,MBE,Real,Real,Real,CC_POT,CC_POT>
        STATIONARY_GO_POT;
      STATIONARY_GO_POT stationaryGoPot(gfsPot,ccPot,gfsPot,ccPot,stationaryLopPot);

      typedef Dune::PDELab::GridOperator<SUB_GFS_CON,SUB_GFS_CON,LOP_CON,MBE,Real,Real,Real,CC_CON,CC_CON> GO_CON0;
      typedef Dune::PDELab::GridOperator<SUB_GFS_CON,SUB_GFS_CON,TLOP_CON,MBE,Real,Real,Real,CC_CON,CC_CON> GO_CON1;
      typedef Dune::PDELab::OneStepGridOperator<GO_CON0,GO_CON1> GO_CON;
      GO_CON0 goCon0_Inside(subGfsCon_Inside,ccCon_Inside,subGfsCon_Inside,ccCon_Inside,lopCon_Inside);
      GO_CON1 goCon1_Inside(subGfsCon_Inside,ccCon_Inside,subGfsCon_Inside,ccCon_Inside,tlopCon_Inside);
      GO_CON goCon_Inside(goCon0_Inside,goCon1_Inside);

      GO_CON0 goCon0_Outside(subGfsCon_Outside,ccCon_Outside,subGfsCon_Outside,ccCon_Outside,lopCon_Outside);
      GO_CON1 goCon1_Outside(subGfsCon_Outside,ccCon_Outside,subGfsCon_Outside,ccCon_Outside,tlopCon_Outside);
      GO_CON goCon_Outside(goCon0_Outside,goCon1_Outside);

      typedef Dune::PDELab::GridOperator<GFS_POT,GFS_POT,LOP_POT,MBE,Real,Real,Real,CC_POT,CC_POT> GO_POT;
      GO_POT goPot(gfsPot,ccPot,gfsPot,ccPot,lopPot);
      // ==============================================================================================


      // ========== Select a linear solver backend ====================================================
      //typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS_CON;

      //LS_CON lsCon(5000,false);
      //LS_POT lsPot(5000,false);
      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
      //LS ls(5000,false);
      typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS_CON;
      LS_CON lsCon_Inside(0); // verbose = 1
      LS_CON lsCon_Outside(0); // verbose = 1
      typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS_POT;
      //typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS_POT;
      //typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS_POT;
      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS_POT;
      LS_POT lsPot(0); // verbose = 1
      // ==============================================================================================

      // ========= Solver for (non-)linear problem per stage ==========================================
      // Stationary Poisson-Boltzmann solver
      typedef Ax1Newton<Acme1Output,STATIONARY_GO_POT,LS_POT,U_POT> STATIONARY_SOLVER_POT;
      std::string stationaryOutputPrefix = "octave/";
      stationaryOutputPrefix += stationaryLopPot.getName();
      STATIONARY_SOLVER_POT stationarySolverPot(acme1Output,stationaryGoPot,lsPot,
          stationaryOutputPrefix.c_str());
      stationarySolverPot.setReassembleThreshold(0.0);
      stationarySolverPot.setVerbosityLevel(0); // (2)
      stationarySolverPot.setAbsoluteLimit(1e-8);
      stationarySolverPot.setReduction(1e-6);
      stationarySolverPot.setMinLinearReduction(1e-4);
      stationarySolverPot.setMaxIterations(25);
      stationarySolverPot.setLineSearchMaxIterations(10);
      stationarySolverPot.setPrintMatrix(false);
      stationarySolverPot.setPrintRhs(false);

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
      //Output::printSingleCoefficientVector(uPot, "uPot initial Poisson");

      try {
        //uPot = 0;
        //stationarySolverPot.apply(uPot);

        // Calculate equilibirum concentrations from Poisson-Boltzmann potential in concentration domains
        //GF_PB_CON gfPoissonBoltzmannCon_Inside(initialCon, dgfPot, subGridView_Inside, physics);
        //GF_PB_CON gfPoissonBoltzmannCon_Outside(initialCon, dgfPot, subGridView_Outside, physics);

        //Dune::PDELab::interpolate(gfPoissonBoltzmannCon_Inside, subGfsCon_Inside, uoldCon_Inside);
        //Dune::PDELab::interpolate(gfPoissonBoltzmannCon_Outside, subGfsCon_Outside, uoldCon_Outside);

        // Transfer subgrid equilibrium concentrations to host grid and update solution vectors
        //UCON_Interpolator::interpolate(subGridView_Inside, gv, subGfsCon_Inside, gfsCon, uoldCon_Inside, uCon);
        //UCON_Interpolator::interpolate(subGridView_Outside, gv, subGfsCon_Outside, gfsCon, uoldCon_Outside, uCon);
        //unewCon_Inside = uoldCon_Inside;
        //unewCon_Outside = uoldCon_Outside;
      } catch (Dune::Exception& e) {
        debug_warn << e.what() << std::endl;
        //Output::printSingleCoefficientVector(uPot, "uPot Poisson-Boltzmann");
        throw e;
      }
      debug_verb << "============= SOLVED initial Poisson problem" << std::endl;

      // Transfer potential vector to subgrid
      UPOT_Restrictor::restrict(gv, subGridView_Inside, gfsPot, subGfsPot_Inside, uPot, uPot_Inside);
      UPOT_Restrictor::restrict(gv, subGridView_Outside, gfsPot, subGfsPot_Outside, uPot, uPot_Outside);


//      Output::printSingleCoefficientVector(uPot,"uPot");
//      Output::printSingleCoefficientVector(uPot_Inside,"uPot_Inside");
//      Output::printSingleCoefficientVector(uPot_Outside,"uPot_Outside");
//
//      Output::printMultipleComponentCoefficientVector(uCon,NUMBER_OF_SPECIES);
//      Output::printMultipleComponentCoefficientVector(unewCon_Inside,NUMBER_OF_SPECIES);
//      Output::printMultipleComponentCoefficientVector(unewCon_Outside,NUMBER_OF_SPECIES);

      typedef Ax1StationaryLinearProblemSolver<GO_CON,LS_CON,USUB_CON> PDESOLVER_CON;
      PDESOLVER_CON pdesolverCon_Inside(goCon_Inside,lsCon_Inside,1e-10,"octave/nernst_planck_inside");
      pdesolverCon_Inside.setPrintMatrix(false);
      pdesolverCon_Inside.setPrintRhs(false);
      pdesolverCon_Inside.setRowPreconditioner(useRowPreconditioner);
      pdesolverCon_Inside.setVerbosity(0); //2
      PDESOLVER_CON pdesolverCon_Outside(goCon_Outside,lsCon_Outside,1e-10,"octave/nernst_planck_outside");
      pdesolverCon_Outside.setPrintMatrix(false);
      pdesolverCon_Outside.setPrintRhs(false);
      pdesolverCon_Outside.setRowPreconditioner(useRowPreconditioner);
      pdesolverCon_Outside.setVerbosity(0); //2
      debug_verb << "Nernst-Planck PDE solvers set up" << std::endl;
      // ==============================================================================================

      // ========== time-stepper ======================================================================
      //Dune::PDELab::Alexander3Parameter<Real> timeStepperCon_Inside;
      //Dune::PDELab::Alexander3Parameter<Real> timeStepperCon_Outside;
      Dune::PDELab::ImplicitEulerParameter<Real> timeStepperCon_Inside;
      Dune::PDELab::ImplicitEulerParameter<Real> timeStepperCon_Outside;
      //Dune::PDELab::ExplicitEulerParameter<Real> timeStepperCon;
      //Dune::PDELab::RK4Parameter<Real> timeStepperCon;
      //Dune::PDELab::HeunParameter<Real> timeStepperCon;
      typedef Dune::PDELab::OneStepMethod<Real,GO_CON,PDESOLVER_CON,USUB_CON,USUB_CON> SOLVER_CON;
      SOLVER_CON solverCon_Inside(timeStepperCon_Inside,goCon_Inside,pdesolverCon_Inside);
      solverCon_Inside.setVerbosityLevel(0); // 2
      SOLVER_CON solverCon_Outside(timeStepperCon_Outside,goCon_Outside,pdesolverCon_Outside);
      solverCon_Outside.setVerbosityLevel(0); // 2
      
      // In case we need access to the pdesolver later on
      const PDESOLVER_CON& pdesolverConRef_Inside = solverCon_Inside.getPDESolver();
      const PDESOLVER_CON& pdesolverConRef_Outside = solverCon_Outside.getPDESolver();
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
        //solutionVectors.printHostGridVectors();
      }
      // ==============================================================================================

      // ========== Initial output ====================================================================
      double debyeLength = physics.getDebyeLength(gfIonicStrength);
      double h_elec = (physics.getParams().xMax() - 0.5 * physics.getParams().dMemb()) / (subGridView_Inside.size(0));

      debug_verb << "" << std::endl;
      debug_verb << "@@@@@@@ minimum Debye length / length scale " << debyeLength / physics.getLengthScale()<< std::endl;
      if(h_elec >= (debyeLength/physics.getLengthScale()))
      {
        debug_warn << "WARNING: Grid does not resolve Debye length [h_elec = " << h_elec << "]!" << std::endl;
      }
      debug_verb << "" << std::endl;

      typename Acme1Output::DiagnosticInfo& diagInfo = acme1Output.getDiagInfo();
      diagInfo.tEquilibrium = tEquilibrium;

      acme1Output.writeStep(time);
      Tools::pecletNumber(subGridView_Inside, parametersCon_Inside, dgfGradPot); // TODO dgfGradPot_Inside?
      Tools::pecletNumber(subGridView_Outside, parametersCon_Outside, dgfGradPot); // TODO dgfGradPot_Outside?
      debug_info  << std::endl << "########## initial output done" << std::endl << std::endl;
      // ==============================================================================================

      // ========= time loop ==========================================================================
      // Setup solution method
      typedef Acme1OperatorSplit<Traits,DirichletValuesCon,DirichletValuesPot,GO_CON,GO_POT,
                                 SOLVER_CON,SOLVER_POT> Acme1OperatorSplit;
      //typedef Acme1Iteration<Traits,DirichletValuesCon,DirichletValuesPot,GO_CON,GO_POT,
      //                                 SOLVER_CON,SOLVER_POT> Acme1OperatorSplit;
      Acme1OperatorSplit acme1OperatorSplit(uCon, uConPrevious, uPot, uoldCon_Inside, unewCon_Inside,
          uoldCon_Outside, unewCon_Outside, uPot_Inside, uPot_Outside,
          dirichletValuesCon_Inside, dirichletValuesCon_Outside, dirichletValuesPot,
          goCon_Inside, goCon_Outside, goPot, solverCon_Inside, solverCon_Outside,
          solverPot, gfsCon, gfsPot, subGfsCon_Inside, subGfsCon_Outside,
          subGfsPot_Inside, subGfsPot_Outside, physics, gv,
          subGridView_Inside, subGridView_Outside, acme1Output);
      acme1OperatorSplit.setReduction(physics.getParams().getReduction());
      acme1OperatorSplit.setAbsLimit(physics.getParams().getAbsLimit());
      acme1OperatorSplit.setMaxIterations(500);
      acme1OperatorSplit.setUseDefect(false);

      Real outputTimeInterval = physics.getParams().getOutputTimeInterval();
      int outputCounter = 1;
      bool printEveryTimeStep = (dt > outputTimeInterval && !physics.getParams().useAdaptiveTimeStep());

      // Number of iterations/time steps since last output
      int iterations = 0;
      int timeSteps = 0;

      typename PHYSICS::ChannelSet& channels = physics.getMembrane().getChannelSet();
      bool channelsInitialized = false;

      while (time<tend-1e-8)
      {
        //debug_verb << "TIME= " << time << std::endl;
        //debug_verb << "outputCounter= " << outputCounter << std::endl;

        // Initialize channels
        if(! channelsInitialized /*&& std::abs(time-tEquilibrium) < 1e-8*/)
        {
          updateChannels(channels, parametersCon_Inside, parametersCon_Outside,
              dgfCon, dgfPot, time, dt, channelsInitialized);
          channelsInitialized = true;
        }

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

        if(physics.getParams().useAdaptiveTimeStep())
        {
          Real safety_factor = 0.1;
          dt = safety_factor * Tools::getTimeStep(gv, physics, dgfGradPot);
          debug_verb << "Calculated time step: " << dt << std::endl;
        }
        physics.setTimeStep(dt);
        diagInfo.dt = dt;
        
        U_POT uoldPot = uPot;

        try{

          // ========= DO ONE TIME STEP ===========
          dt = acme1OperatorSplit.timeStep(time,dt);

        } catch(Dune::Exception& e)
        {
          // Write debug output to file for later inspection
          if(writeIfNotConverged)
          {
            debug_verb << "write debug output" << std::endl;
            for(int i=0; i<diagInfo.maxDiffCon.size(); i++)
            {
              std::stringstream infoStream;
              Real iteration = i;
              infoStream << "iteration: " << std::fixed << std::setw(3) << i;
              Output::gnuplotDoubleAppend("operator_split_debug.dat", iteration, diagInfo.maxDiffCon[i], diagInfo.maxDiffPot[i],
                  physics.getParams().getOutputPrefix(), infoStream.str());

            }
          }
          time += dt;
          // Write out non-converged solution
          acme1Output.writeStep(time);
          Tools::pecletNumber(subGridView_Inside, parametersCon_Inside, dgfGradPot); // TODO dgfGradPot_Inside?
          Tools::pecletNumber(subGridView_Outside, parametersCon_Outside, dgfGradPot); // TODO dgfGradPot_Outside?

          //Output::printMultipleComponentCoefficientVector(uCon, NUMBER_OF_SPECIES);
          throw e;
        }

        USUB_CON uConChange_Inside(unewCon_Inside);
        uConChange_Inside -= uoldCon_Inside;
        USUB_CON uConChange_Outside(unewCon_Outside);
        uConChange_Outside -= uoldCon_Outside;
        U_POT uPotChange(uPot);
        uPotChange -= uoldPot;
        {
          Dune::ios_base_all_saver hackbraten(std::cout);
          debug_info << std::scientific << std::setprecision(16);
          debug_info << "L2  con change inside  = " << uConChange_Inside.base().two_norm();
          debug_info << ", outside = " << uConChange_Outside.base().two_norm() << std::endl;
          debug_info << "MAX con change inside  = " << uConChange_Inside.base().infinity_norm();
          debug_info << ", outside = " << uConChange_Outside.base().infinity_norm() << std::endl;
          debug_info << "L2  pot change = " << uPotChange.base().two_norm() << std::endl;
          debug_info << "MAX pot change = " << uPotChange.base().infinity_norm() << std::endl;
        }

        uoldCon_Inside = unewCon_Inside;
        uoldCon_Outside = unewCon_Outside;

        time += dt;
        iterations += diagInfo.iterations;
        timeSteps++;

        // Update channels
        updateChannels(channels, parametersCon_Inside, parametersCon_Outside,
                      dgfCon, dgfPot, time, dt, true);

        double diffToNextTimeStep = time - outputCounter * outputTimeInterval;
        if (printEveryTimeStep || std::abs(diffToNextTimeStep) < 1e-8 || diffToNextTimeStep > 0)
        {
        	++outputCounter;

        	//Output::printSingleCoefficientVectorDG(uPot, "pot");
          //Output::printMultipleComponentCoefficientVectorDG(unewCon, NUMBER_OF_SPECIES);

        	acme1Output.writeStep(time);
        	Tools::pecletNumber(subGridView_Inside, parametersCon_Inside, dgfGradPot); // TODO dgfGradPot_Inside?
          Tools::pecletNumber(subGridView_Outside, parametersCon_Outside, dgfGradPot); // TODO dgfGradPot_Outside?

        	debug_info << std::endl << "########## output done, time: " << time
        	    << " [" << timeSteps << " steps, average #iterations: " << (iterations/timeSteps) << "] ###########"
        	    << std::endl << std::endl;

        	iterations = 0;
        	timeSteps = 0;
        }
      }
      // ==============================================================================================

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

    //! \brief Prepare channel flux calculation for next time step
    template<typename PARAMETERS_CON, typename DGF_CON, typename DGF_POT>
    void updateChannels(typename PHYSICS::ChannelSet& channels,
        PARAMETERS_CON& parametersCon_Inside, PARAMETERS_CON& parametersCon_Outside,
        DGF_CON& dgfCon, DGF_POT& dgfPot,
        Real time, Real dt, bool channelsInitialized)
    {
      // =========  =========================
      const std::map<int,int>& membraneElements = channels.getMembraneElements();

      std::valarray<Real> potJumps(membraneElements.size());
      std::valarray<typename DGF_CON::Traits::RangeType> conJumps(membraneElements.size());
      std::valarray<typename DGF_CON::Traits::RangeType> conUps(membraneElements.size());

      // Make one time step for channel conductances
      for(ElementIterator eit=gv.template begin<0>(); eit != gv.template end<0>(); ++eit)
      {
        if(physics.isMembrane(*eit))
        {
          const typename ElementIterator::Entity& e = *eit;

          //Membrane interfaces
          int elemIndex = physics.getElementIndex(e);
          int membraneElemIndex = physics.getLocalMembraneElementIndex(elemIndex);

          typename DGF_POT::Traits::RangeType potJump(0.0);
          physics.getMembranePotentialJump(e, dgfPot, potJump);
          potJumps[membraneElemIndex] = potJump;

          typename DGF_CON::Traits::RangeType conJump(0.0), conUp(0.0);
          physics.getMembraneConcentrationJump(e, dgfCon, conJump, conUp);
          conJumps[membraneElemIndex] = conJump;
          conUps[membraneElemIndex] = conUp;

          potJump = physics.convertTo_mV(potJump);

          if(not channelsInitialized)
          {
            // Init channels
            channels.initChannels(physics.getElementIndex(e), potJump);
          } else {
            // Do one time step for channels
            double real_dt = dt * physics.getTimeScale();
            debug_jochen << "Calling channels.timeStep()!" << std::endl;
            channels.timeStep(physics.getElementIndex(e), real_dt, potJump);
          }
        }
      }
      parametersCon_Inside.updatePotJump(potJumps);
      parametersCon_Outside.updatePotJump(potJumps);
      parametersCon_Inside.updateConJump(conJumps);
      parametersCon_Outside.updateConJump(conJumps);
      parametersCon_Inside.updateConUp(conUps);
      parametersCon_Outside.updateConUp(conUps);

      // Permanently use new vector of gating particles
      physics.getMembrane().updateState();

      debug_jochen << "================= Channel DEBUG Kacke ============================================" << std::endl;
      for(std::map<int,int>::const_iterator it = membraneElements.begin();
          it != membraneElements.end(); ++it)
      {
        int i = it->first;
        debug_jochen << "Membrane element #" << i << std::endl;

        for(int k=0; k<channels.size(); ++k)
        {
          debug_jochen << channels.getChannel(k).getName();
          debug_jochen << " - " << channels.getConductance(k, i) << std::endl;
          for(int j=0; j<channels.getChannel(k).numGatingParticles(); j++)
          {
            debug_jochen << "^" << channels.getGatingParticle(k,j,i) << std::endl;
          }
        }
        debug_jochen << std::endl;
      }
      debug_jochen << "=================================================================================" << std::endl;

  // =========================================================================================
    }


    void saveState(double time, double dt, Acme1SolutionVectors<Traits>& solutionVectors, std::string filename)
    {
      Ax1SimulationData<Real> simulationData(gv.size(0), time, dt, filename);

      std::vector<Real> stdVec;

      // Host grid
      solutionVectors.uCon.std_copy_to(stdVec);
      simulationData.addVector("uCon", stdVec);

      solutionVectors.uPot.std_copy_to(stdVec);
      simulationData.addVector("uPot", stdVec);

      // Subrid 'inside'
      solutionVectors.uoldCon_Inside.std_copy_to(stdVec);
      simulationData.addVector("uoldCon_Inside", stdVec);

      solutionVectors.unewCon_Inside.std_copy_to(stdVec);
      simulationData.addVector("unewCon_Inside", stdVec);

      solutionVectors.uPot_Inside.std_copy_to(stdVec);
      simulationData.addVector("uPot_Inside", stdVec);

      // Subgrid 'outside'
      solutionVectors.uoldCon_Outside.std_copy_to(stdVec);
      simulationData.addVector("uoldCon_Outside", stdVec);

      solutionVectors.unewCon_Outside.std_copy_to(stdVec);
      simulationData.addVector("unewCon_Outside", stdVec);

      solutionVectors.uPot_Outside.std_copy_to(stdVec);
      simulationData.addVector("uPot_Outside", stdVec);

      Ax1SimulationState<Real> state;
      state.setData(simulationData);
      state.saveState();
    }

    void loadState(double& time, double& dt, Acme1SolutionVectors<Traits>& solutionVectors, std::string filename)
    {
      Ax1SimulationState<Real> state;
      state.loadState(filename);

      Ax1SimulationData<Real>& simulationData = state.getData();

      assert(gv.size(0) == simulationData.getNElements());

      time = simulationData.getTime();
      dt = simulationData.getDt();

      std::vector<Real> stdVec;

      // Host grid
      stdVec = simulationData.getVector("uCon").data;
      solutionVectors.uCon.std_copy_from(stdVec);

      stdVec = simulationData.getVector("uPot").data;
      solutionVectors.uPot.std_copy_from(stdVec);

      // Subrid 'inside'
      stdVec = simulationData.getVector("uoldCon_Inside").data;
      solutionVectors.uoldCon_Inside.std_copy_from(stdVec);

      stdVec = simulationData.getVector("unewCon_Inside").data;
      solutionVectors.unewCon_Inside.std_copy_from(stdVec);

      stdVec = simulationData.getVector("uPot_Inside").data;
      solutionVectors.uPot_Inside.std_copy_from(stdVec);

      // Subgrid 'outside'
      stdVec = simulationData.getVector("uoldCon_Outside").data;
      solutionVectors.uoldCon_Outside.std_copy_from(stdVec);

      stdVec = simulationData.getVector("unewCon_Outside").data;
      solutionVectors.unewCon_Outside.std_copy_from(stdVec);

      stdVec = simulationData.getVector("uPot_Outside").data;
      solutionVectors.uPot_Outside.std_copy_from(stdVec);
    }

    const GV& gv;
    PHYSICS& physics;

    const SubGV& subGridView_Inside;
    const SubGV& subGridView_Outside;
};

#endif /* DUNE_AX1_ACME1_SETUP_HH */
