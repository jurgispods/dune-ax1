#ifndef DUNE_AX1_ACME0_OPERATOR_FULLY_IMPLICIT_OLD_HH
#define DUNE_AX1_ACME0_OPERATOR_FULLY_IMPLICIT_OLD_HH

#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

#include <dune/ax1/common/constants.hh>

template<class BCT, class Physics>
class Acme0LocalOperator :
  public Dune::PDELab::NumericalJacobianApplyVolume<Acme0LocalOperator<BCT,Physics> >,
  public Dune::PDELab::NumericalJacobianVolume<Acme0LocalOperator<BCT,Physics> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<Acme0LocalOperator<BCT,Physics> >,
  public Dune::PDELab::NumericalJacobianBoundary<Acme0LocalOperator<BCT,Physics> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  
  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doAlphaBoundary = true };
  
  // constructor stores parameters
  Acme0LocalOperator (const BCT& bct_,
                      Physics& physics_,
                      unsigned int intorder_=2)
    : bct(bct_),
      physics(physics_),
      intorder(intorder_)
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // select the two components (assume Galerkin scheme U=V)
    typedef typename LFSU::template Child<0>::Type LFSU_CON;
    typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
    typedef typename LFSU::template Child<1>::Type LFSU_POT;
    
    const LFSU_CON& lfsuCon = lfsu.template getChild<0>();
    const LFSU_POT& lfsuPot = lfsu.template getChild<1>();
    
    // Numer of local power function space must equal number of ion species
    assert(NUMBER_OF_SPECIES == physics.numOfSpecies());

    // domain and range field type (assume both components have same RF)
    typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianType;
    typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename LFSU::Traits::SizeType size_type;
        
    // dimensions
    const int dim = EG::Geometry::dimension;
    const int dimw = EG::Geometry::dimensionworld;
    
    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
    
    // get element index
    int elemIndex = physics.getElementIndex(eg.entity());
    
    // get group index
    int subdomainIndex = physics.getSubdomainIndex( elemIndex );
    
    // get permittivity of the element
    RF permittivity = physics.getPermittivity( subdomainIndex );
    
    // get potential jump through membrane
    RF deltaPot = 0.0;
    if ( physics.isMembrane(subdomainIndex) ) // membrane
    {
      RF potMin = 0.0;
      RF potMax = 0.0;
      for (typename Dune::QuadratureRule<DF,dim>::const_iterator 
           it=rule.begin(); it!=rule.end(); ++it)
      {
        std::vector<RangeType> phiPot(lfsuPot.size());
        lfsuPot.finiteElement().localBasis().evaluateFunction(it->position(),phiPot);
        RF Pot = 0.0;
        for (size_type i=0; i<lfsuPot.size(); i++) Pot += x(lfsuPot,i)*phiPot[i];
        if ( Pot < potMin ) potMin = Pot;
        if ( Pot > potMax ) potMax = Pot;
      }
      deltaPot = std::abs( potMax - potMin );
      debug_verb << "deltaPot = " << deltaPot << std::endl;
    }
    
    
    // current source test
    RF current = 0.0;
    if ( elemIndex == 32 ) current = 0.0;
    
    // loop over quadrature points #################################################################
    
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator 
           it=rule.begin(); it!=rule.end(); ++it)
      {
        
        const Dune::FieldMatrix<DF,dimw,dim>
          jac = eg.geometry().jacobianInverseTransposed(it->position());

        RF factor = it->weight()*eg.geometry().integrationElement(it->position());

        // *************** Conentration part *************************
        
        std::vector<LFSU_SINGLE_CON> v_lfsuSingleCon;
        std::vector<RF> v_singleCon;
        std::vector<Dune::FieldVector<RF,dim> > v_gradSingleCon;
        std::vector<std::vector<RangeType> > v_phiSingleCon;
        std::vector<std::vector<Dune::FieldVector<RF,dim> > > v_gradPhiSingleCon;

        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          const LFSU_SINGLE_CON& lfsuSingleCon = lfsuCon.child(j);

          // evaluate basis functions on reference element
          std::vector<RangeType> phiSingleCon(lfsuSingleCon.size());
          lfsuSingleCon.finiteElement().localBasis().evaluateFunction(it->position(),phiSingleCon);

          // compute Con at integration point
          RF singleCon=0.0;
          // localIndex() maps dof within leaf space to all dofs within given element
          for (size_type i=0; i<lfsuSingleCon.size(); i++) singleCon += x(lfsuSingleCon,i)*phiSingleCon[i];

          // evaluate gradient of basis functions on reference element
          std::vector<JacobianType> jsSingleCon(lfsuSingleCon.size());
          lfsuSingleCon.finiteElement().localBasis().evaluateJacobian(it->position(),jsSingleCon);

          // transform gradients from reference element to real element
          std::vector<Dune::FieldVector<RF,dim> > gradPhiSingleCon(lfsuSingleCon.size());
          for (size_type i=0; i<lfsuSingleCon.size(); i++) jac.mv(jsSingleCon[i][0],gradPhiSingleCon[i]);

          // compute gradient of Con
          Dune::FieldVector<RF,dim> gradSingleCon(0.0);
          for (size_type i=0; i<lfsuSingleCon.size(); i++)
            gradSingleCon.axpy(x(lfsuSingleCon,i),gradPhiSingleCon[i]);
          
          v_lfsuSingleCon.push_back(lfsuSingleCon);
          v_singleCon.push_back(singleCon);
          v_gradSingleCon.push_back(gradSingleCon);
          v_phiSingleCon.push_back(phiSingleCon);
          v_gradPhiSingleCon.push_back(gradPhiSingleCon);
          
        }
        
        // *************** Potential part *************************
        
        // evaluate basis functions on reference element
        std::vector<RangeType> phiPot(lfsuPot.size());
        lfsuPot.finiteElement().localBasis().evaluateFunction(it->position(),phiPot);
        
        // compute Pot at integration point
        RF Pot=0.0;
        // localIndex() maps dof within leaf space to all dofs within given element
        for (size_type i=0; i<lfsuPot.size(); i++) Pot += x(lfsuPot,i)*phiPot[i];
        
        // evaluate gradient of basis functions on reference element
        std::vector<JacobianType> jsPot(lfsuPot.size());
        lfsuPot.finiteElement().localBasis().evaluateJacobian(it->position(),jsPot);

        // transform gradients from reference element to real element
        std::vector<Dune::FieldVector<RF,dim> > gradPhiPot(lfsuPot.size());
        for (size_type i=0; i<lfsuPot.size(); i++) jac.mv(jsPot[i][0],gradPhiPot[i]);
        
        // compute gradient of Pot
        Dune::FieldVector<RF,dim> gradPot(0.0);
        for (size_type i=0; i<lfsuPot.size(); i++) gradPot.axpy(x(lfsuPot,i),gradPhiPot[i]);
        
        // integrate components ****************************************************************
        
        // eq. 0: Nernst-Planck
        
        // To get the diffusion coefficients of the membrane, we first need to solve
        // the coupled system of ODEs for all the channels on this element

        double dt = physics.getTimeStep();
        physics.updateDiffCoeffs(subdomainIndex, elemIndex, dt, deltaPot);


        RF chargeDensity = 0.0;
        
        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          RF diffCoeff = physics.getDiffCoeff( subdomainIndex, j, elemIndex );
          RF valence = physics.getValence( j );
          chargeDensity += valence*v_singleCon[j];
          
          //if ( j == 1 ) current = 0.0; // turn current off if cl
          
          for (size_type i=0; i<v_lfsuSingleCon[j].size(); ++i)
          {
            r.accumulate(v_lfsuSingleCon[j],i,(
              diffCoeff*(v_gradSingleCon[j]+valence*v_singleCon[j]*gradPot)*v_gradPhiSingleCon[j][i]
              -current*v_phiSingleCon[j][i]
            )*factor);
          }
        }
        
        // eq. 1: Poisson
        
        RF poissonConstant = physics.getPoissonConstant();
        
        /*
        Dune::FieldVector<RF,dim> globalpos = eg.geometry().global(it->position());
        if ( std::abs(globalpos[0])<0.12 )
        {
          permittivity=2;
          std::cout << "### " << elemIndex << " " << globalpos[0] << std::endl;
        }
        */
        
        for (size_type i=0; i<lfsuPot.size(); ++i)
        {
          r.accumulate(lfsuPot,i,(
          -gradPot*gradPhiPot[i]+
                    chargeDensity/permittivity*poissonConstant*phiPot[i]
          )*factor);
        }
      }
  }

  
  // boundary integral #############################################################################
  
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // select the two components (assume Galerkin scheme U=V)
    typedef typename LFSU::template Child<0>::Type LFSU_CON;
    typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
    typedef typename LFSU::template Child<1>::Type LFSU_POT;
    
    const LFSU_CON& lfsuCon = lfsu.template getChild<0>();
    const LFSU_POT& lfsuPot = lfsu.template getChild<1>();
    
    // select components of constraints parameters
    typedef typename BCT::template Child<0>::Type BCT_CON;
    typedef typename BCT_CON::template Child<0>::Type BCT_SINGLE_CON;
    typedef typename BCT::template Child<1>::Type BCT_POT;
    
    const BCT_CON& bctCon = bct.template getChild<0>();
    const BCT_POT& bctPot = bct.template getChild<1>();
    
    // some types
    typedef typename LFSU_POT::Traits::FiniteElementType::
    Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU_POT::Traits::FiniteElementType::
    Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU_POT::Traits::FiniteElementType::
    Traits::LocalBasisType::Traits::JacobianType JacobianType;
    typedef typename LFSU_POT::Traits::FiniteElementType::
    Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename LFSU_POT::Traits::SizeType size_type;
    
    // dimensions
    const int dim = IG::dimension;
    const int dimw = IG::Geometry::dimensionworld;
    
    // select quadrature rule for face
    Dune::GeometryType gtface = ig.geometryInInside().type();
    const Dune::QuadratureRule<DF,dim-1>& 
    rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);
    
    // loop over quadrature points and integrate normal flux
    for(typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin();it!=rule.end();++it)
    {
      // Potential /////////////////////////////////////////////////////////////////////////////////
      
      if ( not bctPot.isDirichlet( ig, it->position() ) )
      {
        // position of quadrature point in local coordinates of element 
        Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());
        
        // evaluate basis functions at integration point
        std::vector<RangeType> phiPot(lfsuPot.size());
        lfsuPot.finiteElement().localBasis().evaluateFunction(local,phiPot);
        
        // evaluate flux boundary condition
        //Dune::FieldVector<RF,dim> globalpos = ig.geometry().global(it->position());
        RF J = 0.0;
        
        // integrate J
        RF factor = it->weight()*ig.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsuPot.size(); i++)
          r.accumulate(lfsuPot,i, J*phiPot[i] * factor);
      }
      
      // Concentrations ////////////////////////////////////////////////////////////////////////////
      
      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        const BCT_SINGLE_CON& bctSingleCon = bctCon.child(j);
        
        if ( not bctSingleCon.isDirichlet( ig, it->position() ) )
        {
          if ( bctSingleCon.nonZeroBoundaryFlux( ig, it->position(), j ) )
          {
            const LFSU_SINGLE_CON& lfsuSingleCon = lfsuCon.child(j);
            
            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());
            
            // evaluate basis functions on reference element
            std::vector<RangeType> phiSingleCon(lfsuSingleCon.size());
            lfsuSingleCon.finiteElement().localBasis().evaluateFunction(local,phiSingleCon);
            
            // evaluate flux boundary condition
            //Dune::FieldVector<RF,dim> globalpos = ig.geometry().global(it->position());
            RF J = 0.1;
            
            // integrate J
            RF factor = it->weight()*ig.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsuSingleCon.size(); i++)
              r.accumulate(lfsuSingleCon,i, -J*phiSingleCon[i] * factor);
          }
        }
      }
      
    }
  }
  
private:
  const BCT& bct;
  Physics& physics;
  unsigned int intorder;
};

#endif /* DUNE_AX1_ACME0_OPERATOR_FULLY_IMPLICIT_OLD_HH */
