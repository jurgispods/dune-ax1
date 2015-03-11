#ifndef DUNE_AX1_ACME2_SETUP_HH
#define DUNE_AX1_ACME2_SETUP_HH

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
#include <dune/ax1/acme2/common/acme2_boundary.hh>
#include <dune/ax1/acme2/common/acme2_initial.hh>
#include <dune/ax1/acme2/common/acme2_output.hh>
#include <dune/ax1/acme2/common/acme2_solutionvectors.hh>
#include <dune/ax1/acme2/common/nernst_planck_parameters.hh>
#include <dune/ax1/acme2/common/poisson_parameters.hh>
#include <dune/ax1/acme2/common/poisson_boltzmann_parameters.hh>
#include <dune/ax1/acme2/common/acme2_simulation.hh>
#include <dune/ax1/acme2/operator/convectiondiffusionfem.hh>
#include <dune/ax1/acme2/operator/poisson_boltzmann_operator.hh>
#include <dune/ax1/acme2/operator/acme2_operator_fully_coupled.hh>

#include <dune/ax1/acme2/operator/acme2_toperator.hh>

template<class Grid, class GV, class PHYSICS, class SubGV>
class Acme2Setup
{
  public:
    
    //! Traits class providing function spaces and related data types for the acme2 use case
    template<int NUMBER_OF_SPECIES>
    struct Acme2Traits
    {
      typedef Acme2Setup<Grid,GV,PHYSICS,SubGV> SETUP;

      typedef PHYSICS Physics;

      typedef GV GridView;
      typedef Grid GridType;
      typedef typename GV::Grid::ctype Coord;
      typedef double Real;

      typedef SubGV SubGridView;
      typedef typename SubGV::Grid SubGridType;

#if USE_PARALLEL==1
      typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints CONSTRAINTS;
#else
      typedef Dune::PDELab::ConformingDirichletConstraints CONSTRAINTS;
#endif
      typedef Dune::PDELab::ISTLVectorBackend<AX1_BLOCKSIZE> VBE;

      typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,GV::dimension> FEM_CON;
      typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,GV::dimension> FEM_POT;

      //typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM_CON;
      //typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM_POT;

      // Nernst-Planck GFS (on electrolyte subdomain)
      typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper LexicographicOrdering;
      typedef Dune::PDELab::GridFunctionSpaceComponentBlockwiseMapper<1,1,1> PowerGFSOrdering;
      typedef Dune::PDELab::GridFunctionSpace<SubGV,FEM_CON,CONSTRAINTS,VBE> GFS_SINGLE_CON;
      typedef Dune::PDELab::PowerGridFunctionSpace<
        GFS_SINGLE_CON, NUMBER_OF_SPECIES, PowerGFSOrdering> GFS_CON;

      // Poisson GFS (on electrolyte or membrane subdomain)
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM_POT,CONSTRAINTS,VBE> GFS_POT;

      // Composite GFS for electrolyte
      //typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper CompositeGFSOrdering;
      // This has no effect when using dune-multidomain => ordering always lexicographic!
      typedef Dune::PDELab::GridFunctionSpaceComponentBlockwiseMapper<3,1> CompositeBlockedOrdering;

      // Multidomain GFS
      typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<GridType,VBE,GFS_CON,GFS_POT> MultiGFS;

      // Extract SubGFS for use in gridfunctions; use these after creating MultiGFS!
      typedef typename Dune::PDELab::GridFunctionSubSpace<MultiGFS,0> GFS_CON_SUB;
      typedef typename Dune::PDELab::GridFunctionSubSpace<MultiGFS,1> GFS_POT_SUB;

      // Extract solution vector type
      //typedef typename Dune::PDELab::BackendVectorSelector<GFS_CON,Real>::Type U_CON;
      //typedef typename Dune::PDELab::BackendVectorSelector<GFS_POT,Real>::Type U_POT;
      typedef typename Dune::PDELab::BackendVectorSelector<MultiGFS,Real>::Type U;

      // Various gridfunctions
      typedef Dune::PDELab::VectorDiscreteGridFunction <GFS_CON_SUB,U> DGF_CON_SUB;
      typedef Dune::PDELab::VectorDiscreteGridFunctionGradient <GFS_CON_SUB,U> DGF_CON_GRAD_SUB;

      typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_CON_SUB, PHYSICS> DGF_CON;
      typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_CON_GRAD_SUB, PHYSICS> DGF_CON_GRAD;

      typedef Dune::PDELab::DiscreteGridFunction<GFS_POT_SUB,U> DGF_POT;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_POT_SUB,U> DGF_POT_GRAD;

      typedef ChargeDensityGridFunction<DGF_CON, PHYSICS> GF_CD;
      typedef IonicStrengthGridFunction<DGF_CON, PHYSICS> GF_IS;

      typedef MembraneFluxGridFunction<DGF_CON,DGF_POT,PHYSICS> GF_MEMB_FLUX;

      // Grid functions for initial values (subdomains)
      typedef typename PHYSICS::Traits::INITIAL_CON INITIAL_GF_CON;
      typedef InitialCon<GV,Real,NUMBER_OF_SPECIES,PHYSICS,INITIAL_GF_CON> INITIAL_CON;
      typedef typename PHYSICS::Traits::INITIAL_POT INITIAL_GF_POT;
      typedef InitialPot<GV,Real,PHYSICS,GF_CD,INITIAL_GF_POT> INITIAL_POT;

      typedef Dune::PDELab::CompositeGridFunction<INITIAL_CON,INITIAL_POT> INITIAL_ELEC;

      // Helper grid function to calculate equilibirum concentrations from stationary Poisson-Boltzmann potential
      //typedef PoissonBoltzmannRHSGridFunction<DGF_POT, INITIAL_CON, PHYSICS> GF_PB_RHS;
      //typedef PoissonBoltzmannConcentrationGridFunction<INITIAL_CON, DGF_POT, SubGV, PHYSICS> GF_PB_CON;

      // Boundary values
      typedef typename PHYSICS::Traits::NERNST_PLANCK_BOUNDARY BOUNDARY_CON;
      typedef typename PHYSICS::Traits::POISSON_BOUNDARY BOUNDARY_POT;

      // Parameter classes
      typedef NernstPlanckParameters<GV,Real,PHYSICS,BOUNDARY_CON,GF_MEMB_FLUX> PARAMETERS_CON;
      typedef PoissonParameters<GV,Real,PHYSICS,BOUNDARY_POT> PARAMETERS_POT;
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
    typedef Acme2Traits<NUMBER_OF_SPECIES> Traits;

    static const bool writeIfNotConverged = true;

    Acme2Setup(Grid& grid_, GV& gv_, PHYSICS& physics_,
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
      const Acme2Parameters& params = physics.getParams();

      Real tEquilibrium = physics.getParams().tEquilibrium();
      if(tEquilibrium > time)
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
      // Single-component GFS for one ion species
      typename Traits::GFS_SINGLE_CON gfsSingleCon(elecGV,femCon);
      // Power grid function space for all ion species
      typename Traits::GFS_CON gfsCon(gfsSingleCon);

      // Full GFS for potential on electrolyte and membrane subdomain
      typename Traits::GFS_POT gfsPot(gv,femPot);

      // Multidomain GFS
      //typename Traits::MultiGFS multigfs(grid,gfsCon,gfsPot);
      typename Traits::MultiGFS multigfs(grid,gfsCon,gfsPot);

      // SubGridFunctionSpaces for use in gridfunctions
      typename Traits::GFS_CON_SUB gfsConSub(multigfs);
      typename Traits::GFS_POT_SUB gfsPotSub(multigfs);
      // ==============================================================================================


      // =========== Define solution vectors ==========================================================
      // Coefficient vectors on subdomains
      //typename Traits::U_CON uoldCon(gfsCon,0.0);
      //typename Traits::U_CON unewCon(gfsCon,0.0);
      //debug_jochen << "U_CON.N(): " << unewCon.N() << std::endl;
      //debug_jochen << "U_CON.flatsize(): " << unewCon.flatsize() << std::endl;

      //typename Traits::U_POT uoldPot(gfsPot,0.0);
      //typename Traits::U_POT unewPot(gfsPot,0.0);
      //debug_jochen << "U_POT.N(): " << unewPot.N() << std::endl;
      //debug_jochen << "U_POT.flatsize(): " << unewPot.flatsize() << std::endl;

      typename Traits::U uold(multigfs,0.0);
      typename Traits::U unew = uold;
      debug_jochen << "U.N(): " << unew.N() << std::endl;
      debug_jochen << "U.flatsize(): " << unew.flatsize() << std::endl;
      debug_jochen << "con/pot GFS size: " << gfsConSub.size() << " / " << gfsPotSub.size() << std::endl;

      // Now let's put 'em all into one structure!
      //Acme2SolutionVectors<Traits> solutionVectors(uold, unew, uoldCon, unewCon, uoldPot, unewPot);
      Acme2SolutionVectors<Traits> solutionVectors(uold, unew);
      // ==============================================================================================

      // ====== Create grid functions and interpolate initial values onto coefficient vectors =========
      // Initial concentrations
      typename Traits::INITIAL_GF_CON initialGFCon(gv,physics.getParams());

      // Generic wrapper class for initial grid function
      typename Traits::INITIAL_CON initialCon(gv,physics,initialGFCon);
      initialCon.setTime(time);

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
      typename Traits::INITIAL_POT initialPot(gv,physics,gfChargeDensity,initialGFPot,2*degreePot);
      initialPot.setTime(time);

      typename Traits::DGF_POT dgfPot(gfsPotSub, unew);
      typename Traits::DGF_POT dgfOldPot(gfsPotSub, uold);
      typename Traits::DGF_POT_GRAD dgfGradPot(gfsPotSub, unew);

      // Flag 'true' for updating channel states in each time step
      typename Traits::GF_MEMB_FLUX gfMembFlux(dgfOldCon, dgfOldPot, physics, true, tEquilibrium);

      // Composite initial grid function for the electrolyte subdomain
      typename Traits::INITIAL_ELEC initialElec(initialCon,initialPot);

      // Define and instantiate output class
      //typedef Acme2Output<Traits,GV,SubGV,typename Traits::GFS_CON,
      //    typename Traits::GFS_POT,typename Traits::U_CON,typename Traits::U_POT,PHYSICS> ACME2_OUTPUT;
      typedef Acme2Output<Traits> ACME2_OUTPUT;

      ACME2_OUTPUT acme2Output(gv,elecGV,membGV,multigfs,gfsConSub,gfsPotSub,solutionVectors,
          gfMembFlux,physics,2*degreePot,2*degreeCon);
      // ==============================================================================================


      // ========== Define Parameter classes containing the model problem =============================
      typename Traits::BOUNDARY_CON boundaryCon(physics.getParams(), initialGFCon);
      typename Traits::PARAMETERS_CON parametersCon(gv,physics,boundaryCon,gfMembFlux,tEquilibrium);

      typename Traits::BOUNDARY_POT boundaryPot(physics.getParams(), initialGFPot);
      typename Traits::PARAMETERS_POT parametersPot(physics,boundaryPot);
      // ==============================================================================================

      // ============== BOUNDARY CONDITION TYPES ======================================================
      typedef BCTypeSingleCon<typename Traits::PARAMETERS_CON,PHYSICS> BCType_SINGLE_CON;
      typedef Dune::PDELab::PowerConstraintsParameters<BCType_SINGLE_CON,NUMBER_OF_SPECIES> BCType_CON;
      typedef BCTypePot<typename Traits::PARAMETERS_POT,PHYSICS> BCType_POT;

      typedef Dune::PDELab::CompositeConstraintsParameters<BCType_CON, BCType_POT> BCType_ELEC;

      // Create NUMBER_OF_SPECIES boundary condition type classes
      BCType_CON bctypeCon;
      for(int i=0; i<NUMBER_OF_SPECIES; ++i)
      {
        BCType_SINGLE_CON* bctypeSingleCon = new BCType_SINGLE_CON(gv,parametersCon,physics,i);
        bctypeCon.setChild(i,*bctypeSingleCon);
      }
      BCType_POT bctypePot(gv,parametersPot,physics);

      BCType_ELEC bctypeElec(bctypeCon, bctypePot);
      // ==============================================================================================


      // ============== Define local operators ========================================================
      // ##### Spatial part #####
      // ### Standard FEM (CG) fully-coupled Poisson-Nernst-Planck local operator on electrolyte subdomain
      typedef Acme2OperatorFullyCoupled<typename Traits::PARAMETERS_CON, typename Traits::PARAMETERS_POT,
          typename Traits::FEM_CON, typename Traits::FEM_POT> LOP_ELEC;
      LOP_ELEC lopElec(parametersCon, parametersPot);

      // ### Standard FEM (CG) Poisson local operator on membrane subdomain
      typedef Dune::PDELab::ConvectionDiffusionFEM<typename Traits::PARAMETERS_POT,typename Traits::FEM_POT> LOP_MEMB;
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

      // Here, initial values coming from the both subdomain initial GFs are written to uold
      Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,uold,initialElec,elecSubProblem);
      Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,uold,initialPot,membSubProblem);
      unew = uold;
      // ==============================================================================================

      // ========== Select a linear solver backend ====================================================
#if USE_PARALLEL==1
      //typedef Dune::PDELab::ISTLBackend_NOVLP_CG_NOPREC<typename Traits::MultiGFS> LS;
      typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<IGO> LS;
      LS ls(multigfs);
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
      LS ls(1,1.0,5000,5);

      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<IGO,true> LS; // nö / ja
      //LS ls(5000,5);
      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SOR<IGO,true> LS; // ? / ja
      //LS ls(5000,5);
      //typedef Dune::PDELab::ISTLBackend_SEQ_LS_AMG_SSOR<IGO> LS; // nö / ja
      //LS ls(5000,5);
      //typedef Dune::PDELab::ISTLBackend_SEQ_LS_AMG_SOR<IGO> LS; //nö / ja
      //LS ls(5000,5);

      //typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
      //LS ls(5); // verbose = 1

      debug_info << "Using linear solver " << Tools::getTypeName(ls) << std::endl;

      //typedef typename FAKE_IGO::Traits::Jacobian M;
      //typename M::BaseT MatrixType ;
      //typedef typename FAKE_IGO::Traits::Domain FAKE_V;
      //typedef typename IGO::Traits::Domain REAL_V;
      //REAL_V v1(unew);
      //FAKE_V v2(unew);
      //debug_jochen << "M: " << Tools::getTypeName(bla) << std::endl;
      //debug_jochen << "V: " << Tools::getTypeName(v1) << std::endl;
      //debug_jochen << "V: " << Tools::getTypeName(v2) << std::endl;
      //debug_jochen << "VectorType: " << Tools::getTypeName(v2) << std::endl;

#endif
#if 1
      // ==============================================================================================

      // ========= Solver for nonlinear problem per stage ==========================================
      typedef Ax1Newton<ACME2_OUTPUT,IGO,LS,typename Traits::U> PDESOLVER;
      PDESOLVER pdesolver(acme2Output,igo,ls,"octave/acme2");
      pdesolver.setReassembleThreshold(0.0);
      pdesolver.setVerbosityLevel(5); // 2
      pdesolver.setReduction(params.getReduction()); // 1e-10
      pdesolver.setAbsoluteLimit(params.getAbsLimit()); // 1e-10
      pdesolver.setMinLinearReduction(1e-5); // has no effect when using direct solver like SuperLU
      pdesolver.setMaxIterations(50);
      pdesolver.setLineSearchStrategy(PDESOLVER::hackbuschReuskenAcceptBest);
      pdesolver.setLineSearchMaxIterations(50); // 10
      pdesolver.setPrintMatrix(params.doPrintMatrix());
      pdesolver.setPrintRhs(params.doPrintRhs());
      pdesolver.setRowPreconditioner(params.useRowNormPreconditioner());
      pdesolver.setReorderMatrix(params.doReorderMatrix());
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
        acme2Output.loadState(time,dt,loadFilename);
        debug_info << "======================================================================================"
            << std::endl;
        debug_info << "Loaded simulation state from file " << loadFilename << std::endl;
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


        if(not params.doContinueSimulation())
        {
          time = 0.0;
        } else {
          if(tend < time)
            DUNE_THROW(Dune::Exception, "Start time (" << time << ") > end time (" << tend << ")!");
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
      //Tools::compositeToChildCoefficientVector(multigfs, unew, unewCon, 0);
      //uoldCon = unewCon;
      //Tools::compositeToChildCoefficientVector(multigfs, unew, unewPot, 1);
      //uoldPot = unewPot;

      //Output::printSingleCoefficientVector(unewCon, "uCon");
      //Output::printSingleCoefficientVector(unewPot, "uPot");
      //Output::printSingleCoefficientVector(unew, "u");

      double debyeLength, h_x, h_y;
      physics.getDebyeLength(gfIonicStrength, debyeLength, h_x, h_y);

      debug_info << std::endl;
      debug_info << "@@@@@@@ minimum Debye length / length scale " << debyeLength / physics.getLengthScale()
          << std::endl;
      debug_info << "@@@@@@@ h_x = " << h_x << ", h_y = " << h_y << std::endl;
      if(h_x >= (debyeLength/physics.getLengthScale()) || h_y >= (debyeLength/physics.getLengthScale()))
      {
        debug_warn << "WARNING: Grid does not resolve Debye length [h_x = "
            << h_x << ", h_y = " << h_y << "]!" << std::endl;
      }
      debug_verb << "" << std::endl;

      debug_jochen << "Diffusion length for this initial dt: "
          << 2*std::sqrt(physics.getDiffCoeff(Na, 0) * dt)
          << std::endl;
      // ==============================================================================================

      // ========== Run simulation ====================================================================
      Acme2Simulation<Traits,ACME2_OUTPUT>::run(time, dt, dtstart, tend, tEquilibrium,
          physics, gv, membGV,
          solver,
          parametersCon,
          multigfs,
          uold, unew,
          //uoldCon, unewCon, uoldPot, unewPot,
          gfMembFlux, dgfCon, dgfPot, dgfGradPot,
          acme2Output, solutionVectors);
      // ==============================================================================================

      // ======= Save simulation state ================================================================
      if(physics.getParams().doSaveState())
      {
        std::string saveFilename = physics.getParams().getSaveFilename();
        acme2Output.saveState(time,dt,saveFilename);
        debug_info << "=============================================================" << std::endl;
        debug_info << "Saved simulation state to file " << saveFilename << std::endl;
        debug_info << "=============================================================" << std::endl;
      }
      // ==============================================================================================

      debug_info << "Total Acme2Output time: " << acme2Output.getTotalOutputTime() << "s" << std::endl;
#endif
    }

  private:


    Grid& grid;

    GV& gv;
    PHYSICS& physics;

    SubGV& elecGV;
    SubGV& membGV;
};

#endif /* DUNE_AX1_ACME2_SETUP_HH */
