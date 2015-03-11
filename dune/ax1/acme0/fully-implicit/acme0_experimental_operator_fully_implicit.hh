/*
 * acme0_operator_fully_implicit_old.hh
 *
 *  Created on: Nov 24, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME0_EXPERIMENTAL_OPERATOR_FULLY_IMPLICIT_HH
#define DUNE_AX1_ACME0_EXPERIMENTAL_OPERATOR_FULLY_IMPLICIT_HH

#include <dune/ax1/common/constants.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

/**
 * \brief This local operator is designed for the fully-coupled solution of the Nernst-Planck/Poisson system
 * using standard finite elements
 */
template<typename ParamCon, typename ParamPot, typename FiniteElementMapCon, typename FiniteElementMapPot>
class Acme0ExperimentalOperatorFullyImplicit
      : public Dune::PDELab::NumericalJacobianVolume<Acme0ExperimentalOperatorFullyImplicit
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot> >,
        public Dune::PDELab::NumericalJacobianApplyVolume<Acme0ExperimentalOperatorFullyImplicit
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot> >,
        public Dune::PDELab::NumericalJacobianBoundary<Acme0ExperimentalOperatorFullyImplicit
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot> >,
        public Dune::PDELab::NumericalJacobianApplyBoundary<Acme0ExperimentalOperatorFullyImplicit
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot> >,
        public Dune::PDELab::FullSkeletonPattern,
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::LocalOperatorDefaultFlags,
        public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename ParamPot::Traits::RangeFieldType>
{
    enum { dim = ParamPot::Traits::GridViewType::dimension };

    typedef typename ParamPot::Traits::RangeFieldType Real;

  public:
    // pattern assembly flags
    enum { doPatternVolume = true};
    enum { doPatternSkeleton = false };

    // residual assembly flags
    enum { doAlphaVolume  = true };
    enum { doAlphaSkeleton  = false};
    enum { doAlphaBoundary  = true};
    enum { doLambdaVolume  = false};

    //! constructor: pass parameter object
    Acme0ExperimentalOperatorFullyImplicit (ParamCon& paramCon_, ParamPot& paramPot_, int intorderadd_=0)
      : paramCon(paramCon_),
        paramPot(paramPot_),
        intorderadd(intorderadd_)
    {}


    // volume integral depending on test and ansatz functions
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {
      // Extract leaf local function spaces from type tree
      typedef typename LFSU::template Child<0>::Type LFSU_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_CON;
      typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV_CON::template Child<0>::Type LFSV_SINGLE_CON;
      typedef typename LFSU::template Child<1>::Type LFSU_POT;
      typedef typename LFSV::template Child<1>::Type LFSV_POT;

      const LFSU_CON& lfsuCon = lfsu.template getChild<0>();
      const LFSV_CON& lfsvCon = lfsv.template getChild<0>();
      const LFSU_POT& lfsuPot = lfsu.template getChild<1>();
      const LFSV_POT& lfsvPot = lfsv.template getChild<1>();

      // This operator works for identical finite element spaces for con and pot only!
      assert(lfsuPot.finiteElement().localBasis().order()
          == lfsuCon.child(0).finiteElement().localBasis().order());

      //Extract all types from LFSU_POT, assert local function spaces for concentrations are identical!
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
      const int dim = EG::Geometry::dimension;
      const int dimw = EG::Geometry::dimensionworld;

      // select quadrature rule
      Dune::GeometryType gt = eg.geometry().type();
      const int intorder = intorderadd + 2*lfsuPot.finiteElement().localBasis().order();
      const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

      // evaluate diffusion tensor at cell center, assume it is constant over elements
      typename ParamPot::Traits::PermTensorType tensorPot;
      std::vector<typename ParamCon::Traits::PermTensorType> tensorCon;
      Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);

      tensorPot = paramPot.A(eg.entity(),localcenter);
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        paramCon.setIonSpecies(j);
        typename ParamCon::Traits::PermTensorType tensorSingleCon;
        tensorSingleCon = paramCon.A(eg.entity(),localcenter);
        tensorCon.push_back(tensorSingleCon);
      }

      // loop over quadrature points
      for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate basis functions
        std::vector<RangeType> phi(lfsuPot.size());
        lfsuPot.finiteElement().localBasis().evaluateFunction(it->position(),phi);
        //const std::vector<RangeType>& phi = cache.evaluateFunction(it->position(),lfsuPot.finiteElement().localBasis());

        // evaluate u
        // Poisson part: potential
        RF uPot=0.0;
        for (size_type i=0; i<lfsuPot.size(); i++)
          uPot += x(lfsuPot,i)*phi[i];

        //debug_jochen << "uPot = " << uPot << std::endl;

        // Nernst-Planck part: concentrations
        std::vector<RF> uCon;
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          RF uSingleCon=0.0;
          for (size_type i=0; i<lfsuCon.child(j).size(); i++)
            uSingleCon += x(lfsuCon.child(j),i)*phi[i];
          uCon.push_back(uSingleCon);
          //debug_jochen << "uCon[ = " << j << "] = "<< uSingleCon << std::endl;
        }


        // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
        std::vector<JacobianType> js(lfsuPot.size());
        lfsuPot.finiteElement().localBasis().evaluateJacobian(it->position(),js);
        //const std::vector<JacobianType>& js = cache.evaluateJacobian(it->position(),lfsuPot.finiteElement().localBasis());

        // transform gradients of shape functions to real element
        const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsuPot.size());
        for (size_type i=0; i<lfsuPot.size(); i++)
          jac.mv(js[i][0],gradphi[i]);

        // compute gradient of u
        // Poisson part: potential
        Dune::FieldVector<RF,dim> graduPot(0.0);
        for (size_type i=0; i<lfsuPot.size(); i++)
          graduPot.axpy(x(lfsuPot,i),gradphi[i]);

        //debug_jochen << "graduPot = " << graduPot << std::endl;

        // Nernst-Planck part: concentrations
        std::vector<Dune::FieldVector<RF,dim> > graduCon;
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          Dune::FieldVector<RF,dim> graduSingleCon(0.0);
          for (size_type i=0; i<lfsuCon.child(j).size(); i++)
            graduSingleCon.axpy(x(lfsuCon.child(j),i),gradphi[i]);
          graduCon.push_back(graduSingleCon);
          //debug_jochen << "graduCon[" << j << "] = " << graduCon[j] << std::endl;
        }

        // compute A * gradient of u
        // Poisson part: potential
        Dune::FieldVector<RF,dim> AgraduPot(0.0);
        tensorPot.umv(graduPot,AgraduPot);

        // Nernst-Planck part: concentrations
        std::vector<Dune::FieldVector<RF,dim> > AgraduCon;
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          Dune::FieldVector<RF,dim> AgraduSingleCon(0.0);
          tensorCon[j].umv(graduCon[j],AgraduSingleCon);
          AgraduCon.push_back(AgraduSingleCon);
        }

        // evaluate velocity field, sink term and source term
        typename ParamPot::Traits::RangeType bPot = paramPot.b(eg.entity(),it->position());
        typename ParamPot::Traits::RangeFieldType cPot = paramPot.c(eg.entity(),it->position());
        typename ParamPot::Traits::RangeFieldType fPot = paramPot.f(eg.entity(),it->position(),uCon); // depending on cd

        std::vector<typename ParamCon::Traits::RangeType> bCon;
        std::vector<typename ParamCon::Traits::RangeFieldType> cCon;
        std::vector<typename ParamCon::Traits::RangeFieldType> fCon;
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          paramCon.setIonSpecies(j);
          typename ParamCon::Traits::RangeType bSingleCon = paramCon.b(eg.entity(),it->position(),graduPot); // depending on gradPot
          typename ParamCon::Traits::RangeFieldType cSingleCon = paramCon.c(eg.entity(),it->position());
          typename ParamCon::Traits::RangeFieldType fSingleCon = paramCon.f(eg.entity(),it->position());
          bCon.push_back(bSingleCon);
          cCon.push_back(cSingleCon);
          fCon.push_back(fSingleCon);
        }

        // integrate (A grad u)*grad phi_i - u b*grad phi_i + c*u*phi_i
        RF factor = it->weight() * eg.geometry().integrationElement(it->position());

        for (size_type i=0; i<lfsuPot.size(); i++)
         r.accumulate(lfsuPot,i,( AgraduPot*gradphi[i] - uPot*(bPot*gradphi[i]) + (cPot*uPot-fPot)*phi[i] )*factor);

        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          for (size_type i=0; i<lfsuCon.child(j).size(); i++)
           r.accumulate(lfsuCon.child(j),i,( AgraduCon[j]*gradphi[i] - uCon[j]*(bCon[j]*gradphi[i])
               + (cCon[j]*uCon[j]-fCon[j])*phi[i] )*factor);
        }

      }
    }


    /*

    // jacobian of volume term
    template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          M& mat) const
    {
      // Extract leaf local function spaces from type tree
      typedef typename LFSU::template Child<0>::Type LFSU_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_CON;
      typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV_CON::template Child<0>::Type LFSV_SINGLE_CON;
      typedef typename LFSU::template Child<1>::Type LFSU_POT;
      typedef typename LFSV::template Child<1>::Type LFSV_POT;

      const LFSU_CON& lfsuCon = lfsu.template getChild<0>();
      const LFSV_CON& lfsvCon = lfsv.template getChild<0>();
      const LFSU_POT& lfsuPot = lfsu.template getChild<1>();
      const LFSV_POT& lfsvPot = lfsv.template getChild<1>();

      // Nernst-Planck part: concentrations
      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        const LFSU_SINGLE_CON& lfsuSingleCon = lfsuCon.child(j);
        const LFSV_SINGLE_CON& lfsvSingleCon = lfsvCon.child(j);
        paramCon.setIonSpecies(j);
        opCon.jacobian_volume(eg, lfsuSingleCon, x, lfsvSingleCon, mat);
      }

      // Poisson part: potential
      opPot.jacobian_volume(eg, lfsuPot, x, lfsvPot, mat);
    }

    // skeleton integral depending on test and ansatz functions
    // each face is only visited ONCE!
    template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_skeleton (const IG& ig,
                         const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                         const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                         R& r_s, R& r_n) const
    {
      // Extract leaf local function spaces from type tree
      typedef typename LFSU::template Child<0>::Type LFSU_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_CON;
      typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV_CON::template Child<0>::Type LFSV_SINGLE_CON;
      typedef typename LFSU::template Child<1>::Type LFSU_POT;
      typedef typename LFSV::template Child<1>::Type LFSV_POT;

      const LFSU_CON& lfsuCon_s = lfsu_s.template getChild<0>();
      const LFSV_CON& lfsvCon_s = lfsv_s.template getChild<0>();
      const LFSU_CON& lfsuCon_n = lfsu_n.template getChild<0>();
      const LFSV_CON& lfsvCon_n = lfsv_n.template getChild<0>();
      const LFSU_POT& lfsuPot_s = lfsu_s.template getChild<1>();
      const LFSV_POT& lfsvPot_s = lfsv_s.template getChild<1>();
      const LFSU_POT& lfsuPot_n = lfsu_n.template getChild<1>();
      const LFSV_POT& lfsvPot_n = lfsv_n.template getChild<1>();

      // Check if we are on the membrane boundary
      bool isMembrane_s = paramCon.getPhysics().isMembrane(*ig.inside());
      bool isMembrane_n = paramCon.getPhysics().isMembrane(*ig.outside());
      bool membrane = (isMembrane_s || isMembrane_n);
      bool membraneBoundary = (isMembrane_s || isMembrane_n) && not (isMembrane_s && isMembrane_n);

      // No flux on the membrane boundary or within membrane => skip alpha_skeleton!
      if(not membrane)
      {
        // Nernst-Planck part: concentrations
        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          const LFSU_SINGLE_CON& lfsuSingleCon_s = lfsuCon_s.child(j);
          const LFSV_SINGLE_CON& lfsvSingleCon_s = lfsvCon_s.child(j);
          const LFSU_SINGLE_CON& lfsuSingleCon_n = lfsuCon_n.child(j);
          const LFSV_SINGLE_CON& lfsvSingleCon_n = lfsvCon_n.child(j);
          paramCon.setIonSpecies(j);
          opCon.alpha_skeleton(ig, lfsuSingleCon_s, x_s, lfsvSingleCon_s,
                                lfsuSingleCon_n, x_n, lfsvSingleCon_n, r_s, r_n);
        }
      }

      // Poisson part: potential
      opPot.alpha_skeleton(ig, lfsuPot_s, x_s, lfsvPot_s,
                               lfsuPot_n, x_n, lfsvPot_n, r_s, r_n);
    }

    template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_skeleton (const IG& ig,
                            const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                            const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                            M& mat_ss, M& mat_sn,
                            M& mat_ns, M& mat_nn) const
    {
      // Extract leaf local function spaces from type tree
      typedef typename LFSU::template Child<0>::Type LFSU_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_CON;
      typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV_CON::template Child<0>::Type LFSV_SINGLE_CON;
      typedef typename LFSU::template Child<1>::Type LFSU_POT;
      typedef typename LFSV::template Child<1>::Type LFSV_POT;

      const LFSU_CON& lfsuCon_s = lfsu_s.template getChild<0>();
      const LFSV_CON& lfsvCon_s = lfsv_s.template getChild<0>();
      const LFSU_CON& lfsuCon_n = lfsu_n.template getChild<0>();
      const LFSV_CON& lfsvCon_n = lfsv_n.template getChild<0>();
      const LFSU_POT& lfsuPot_s = lfsu_s.template getChild<1>();
      const LFSV_POT& lfsvPot_s = lfsv_s.template getChild<1>();
      const LFSU_POT& lfsuPot_n = lfsu_n.template getChild<1>();
      const LFSV_POT& lfsvPot_n = lfsv_n.template getChild<1>();

      // Check if we are on the membrane boundary
      bool isMembrane_s = paramCon.getPhysics().isMembrane(*ig.inside());
      bool isMembrane_n = paramCon.getPhysics().isMembrane(*ig.outside());
      bool membrane = (isMembrane_s || isMembrane_n);
      bool membraneBoundary = (isMembrane_s || isMembrane_n) && not (isMembrane_s && isMembrane_n);

      // No flux on the membrane boundary or within membrane => skip alpha_skeleton!
      if(not membrane)
      {
        // Nernst-Planck part: concentrations
        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          const LFSU_SINGLE_CON& lfsuSingleCon_s = lfsuCon_s.child(j);
          const LFSV_SINGLE_CON& lfsvSingleCon_s = lfsvCon_s.child(j);
          const LFSU_SINGLE_CON& lfsuSingleCon_n = lfsuCon_n.child(j);
          const LFSV_SINGLE_CON& lfsvSingleCon_n = lfsvCon_n.child(j);
          paramCon.setIonSpecies(j);
          opCon.jacobian_skeleton(ig, lfsuSingleCon_s, x_s, lfsvSingleCon_s,
                                lfsuSingleCon_n, x_n, lfsvSingleCon_n,
                                mat_ss, mat_sn, mat_ns, mat_nn);
        }
      }

      // Poisson part: potential
      opPot.jacobian_skeleton(ig, lfsuPot_s, x_s, lfsvPot_s,
                               lfsuPot_n, x_n, lfsvPot_n,
                               mat_ss, mat_sn, mat_ns, mat_nn);
    }
    */

    // boundary integral depending on test and ansatz functions
    // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
    template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_boundary (const IG& ig,
                         const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                         R& r_s) const
    {
      // Extract leaf local function spaces from type tree
      typedef typename LFSU::template Child<0>::Type LFSU_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_CON;
      typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV_CON::template Child<0>::Type LFSV_SINGLE_CON;
      typedef typename LFSU::template Child<1>::Type LFSU_POT;
      typedef typename LFSV::template Child<1>::Type LFSV_POT;

      const LFSU_CON& lfsuCon_s = lfsu_s.template getChild<0>();
      const LFSV_CON& lfsvCon_s = lfsv_s.template getChild<0>();
      const LFSU_POT& lfsuPot_s = lfsu_s.template getChild<1>();
      const LFSV_POT& lfsvPot_s = lfsv_s.template getChild<1>();

      // This operator works for identical finite element spaces for con and pot only!
      assert(lfsuPot_s.finiteElement().localBasis().order()
          == lfsuCon_s.child(0).finiteElement().localBasis().order());

      //Extract all types from LFSU_POT, assert local function spaces for concentrations are identical!
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
      const int dimw = IG::dimensionworld;

      // This operator works for identical finite element spaces for con and pot only!
      assert(lfsuPot_s.finiteElement().localBasis().order()
          == lfsuCon_s.child(0).finiteElement().localBasis().order());

      // select quadrature rule
      Dune::GeometryType gtface = ig.geometryInInside().type();
      const int intorder = intorderadd+2*lfsuPot_s.finiteElement().localBasis().order();
      const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

      // evaluate boundary condition type
      Dune::FieldVector<DF,dim-1> facecenterlocal = Dune::GenericReferenceElements<DF,dim-1>::general(gtface).position(0,0);
      Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type bctypePot;
      bctypePot = paramPot.bctype(ig.intersection(),facecenterlocal);
      std::vector<Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type> bctypeCon;
      bool allDirichlet = (bctypePot==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet);
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        paramCon.setIonSpecies(j);
        Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type bctypeSingleCon
          = paramCon.bctype(ig.intersection(),facecenterlocal);
        allDirichlet = (allDirichlet && bctypeSingleCon==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet);
        bctypeCon.push_back(bctypeSingleCon);
      }

      // skip rest if we have only Dirichlet boundarys
      if (allDirichlet) return;

      // loop over quadrature points and integrate normal flux
      for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        {
          // position of quadrature point in local coordinates of element
          Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

          // evaluate shape functions (assume Galerkin method)
          std::vector<RangeType> phi(lfsuPot_s.size());
          lfsuPot_s.finiteElement().localBasis().evaluateFunction(local,phi);
          //const std::vector<RangeType>& phi = cache.evaluateFunction(local,lfsu_s.finiteElement().localBasis());

          if (bctypePot==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann)
            {
              // evaluate flux boundary condition
              typename ParamPot::Traits::RangeFieldType jPot = paramPot.j(ig.intersection(),it->position());

              // integrate j
              RF factor = it->weight()*ig.geometry().integrationElement(it->position());
              for (size_type i=0; i<lfsuPot_s.size(); i++)
                r_s.accumulate(lfsuPot_s,i,jPot*phi[i]*factor);
            }

          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            paramCon.setIonSpecies(j);
            if (bctypeCon[j]==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann)
            {
              // evaluate flux boundary condition
              typename ParamCon::Traits::RangeFieldType jSingleCon = paramCon.j(ig.intersection(),it->position());

              // integrate j
              RF factor = it->weight()*ig.geometry().integrationElement(it->position());
              for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
                r_s.accumulate(lfsuCon_s.child(j),i,jSingleCon*phi[i]*factor);
            }
          }

          if (bctypePot==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow)
            {
            // Throw exception for outflow BC, we don't want to use it in our application
            DUNE_THROW(Dune::Exception, "Outflow boundary not implemented!");
              // evaluate u
              RF uPot=0.0;
              for (size_type i=0; i<lfsuPot_s.size(); i++)
                uPot += x_s(lfsuPot_s,i)*phi[i];

              // evaluate velocity field and outer unit normal
              typename ParamPot::Traits::RangeType bPot = paramPot.b(*(ig.inside()),local);
              const Dune::FieldVector<DF,dim> n = ig.unitOuterNormal(it->position());

              // evaluate outflow boundary condition
              typename ParamPot::Traits::RangeFieldType oPot = paramPot.o(ig.intersection(),it->position());

              // integrate o
              RF factor = it->weight()*ig.geometry().integrationElement(it->position());
              for (size_type i=0; i<lfsu_s.size(); i++)
                r_s.accumulate(lfsuPot_s,i,( (bPot*n)*uPot + oPot)*phi[i]*factor);
            }
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            paramCon.setIonSpecies(j);
            if (bctypeCon[j]==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow)
              {
              // Throw exception for outflow BC, we don't want to use it in our application
              DUNE_THROW(Dune::Exception, "Outflow boundary not implemented!");


                // evaluate u
                RF uCon=0.0;
                for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
                  uCon += x_s(lfsuCon_s.child(j),i)*phi[i];

                // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
                std::vector<JacobianType> js(lfsuPot_s.size());
                lfsuPot_s.finiteElement().localBasis().evaluateJacobian(local,js);
                //const std::vector<JacobianType>& js = cache.evaluateJacobian(it->position(),lfsuPot.finiteElement().localBasis());

                // transform gradients of shape functions to real element
                const Dune::FieldMatrix<DF,dimw,dim> jac = ig.inside()->geometry().jacobianInverseTransposed(local);
                std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsuPot_s.size());
                for (size_type i=0; i<lfsuPot_s.size(); i++)
                  jac.mv(js[i][0],gradphi[i]);

                // compute gradient of u
                // Poisson part: potential
                Dune::FieldVector<RF,dim> graduPot(0.0);
                for (size_type i=0; i<lfsuPot_s.size(); i++)
                  graduPot.axpy(x_s(lfsuPot_s,i),gradphi[i]);

                // evaluate velocity field and outer unit normal
                typename ParamCon::Traits::RangeType bCon = paramCon.b(*(ig.inside()),local,graduPot);
                const Dune::FieldVector<DF,dim> n = ig.unitOuterNormal(it->position());

                // evaluate outflow boundary condition
                typename ParamCon::Traits::RangeFieldType oCon = paramCon.o(ig.intersection(),it->position());

                // integrate o
                RF factor = it->weight()*ig.geometry().integrationElement(it->position());
                for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
                  r_s.accumulate(lfsuCon_s.child(j),i,( (bCon*n)*uCon + oCon)*phi[i]*factor);
              }
          }
        }
    }


/*
    template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_boundary (const IG& ig,
                            const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                            M& mat_ss) const
    {
      // Extract leaf local function spaces from type tree
      typedef typename LFSU::template Child<0>::Type LFSU_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_CON;
      typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV_CON::template Child<0>::Type LFSV_SINGLE_CON;
      typedef typename LFSU::template Child<1>::Type LFSU_POT;
      typedef typename LFSV::template Child<1>::Type LFSV_POT;

      const LFSU_CON& lfsuCon_s = lfsu_s.template getChild<0>();
      const LFSV_CON& lfsvCon_s = lfsv_s.template getChild<0>();
      const LFSU_POT& lfsuPot_s = lfsu_s.template getChild<1>();
      const LFSV_POT& lfsvPot_s = lfsv_s.template getChild<1>();

      // Nernst-Planck part: concentrations
      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        const LFSU_SINGLE_CON& lfsuSingleCon_s = lfsuCon_s.child(j);
        const LFSV_SINGLE_CON& lfsvSingleCon_s = lfsvCon_s.child(j);
        paramCon.setIonSpecies(j);
        opCon.jacobian_boundary(ig, lfsuSingleCon_s, x_s, lfsvSingleCon_s, mat_ss);
      }

      // Poisson part: potential
      opPot.jacobian_boundary(ig, lfsuPot_s, x_s, lfsvPot_s, mat_ss);
    }

    // volume integral depending only on test functions
    template<typename EG, typename LFSV, typename R>
    void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
    {
      // Extract leaf local function spaces from type tree
      typedef typename LFSV::template Child<0>::Type LFSV_CON;
      typedef typename LFSV_CON::template Child<0>::Type LFSV_SINGLE_CON;
      typedef typename LFSV::template Child<1>::Type LFSV_POT;

      const LFSV_CON& lfsvCon = lfsv.template getChild<0>();
      const LFSV_POT& lfsvPot = lfsv.template getChild<1>();

      // Nernst-Planck part: concentrations
      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        const LFSV_SINGLE_CON& lfsvSingleCon = lfsvCon.child(j);
        paramCon.setIonSpecies(j);
        opCon.lambda_volume(eg, lfsvSingleCon, r);
      }

      // Poisson part: potential
      opPot.lambda_volume(eg, lfsvPot, r);
    }
    */

    //! set time in parameter class
    void setTime (double t)
    {
      paramCon.setTime(t);
      paramPot.setTime(t);
    }

  private:
    ParamCon& paramCon;  // Nernst-Planck parameter class
    ParamPot& paramPot;  // Poisson parameter class
    int intorderadd;
    typedef typename FiniteElementMapPot::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;
    Dune::PDELab::LocalBasisCache<LocalBasisType> cache;

};



#endif /* DUNE_AX1_ACME0_EXPERIMENTAL_OPERATOR_FULLY_IMPLICIT_HH */
