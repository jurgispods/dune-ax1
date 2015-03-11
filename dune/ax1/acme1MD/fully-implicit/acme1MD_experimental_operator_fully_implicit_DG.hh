/*
 * acme1MD_operator_fully_implicit_old.hh
 *
 *  Created on: Nov 24, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME1MD_EXPERIMENTAL_OPERATOR_FULLY_IMPLICIT_DG_HH
#define DUNE_AX1_ACME1MD_EXPERIMENTAL_OPERATOR_FULLY_IMPLICIT_DG_HH

#include <dune/ax1/common/constants.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#define USECACHE 0


/**
 * \brief This local operator is designed for the fully-coupled solution of the Nernst-Planck/Poisson system
 * using standard finite elements
 */
template<typename ParamCon, typename ParamPot, typename FiniteElementMapCon, typename FiniteElementMapPot>
class Acme1MDExperimentalOperatorFullyImplicitDG
      : public Dune::PDELab::NumericalJacobianVolume<Acme1MDExperimentalOperatorFullyImplicitDG
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot> >,
        public Dune::PDELab::NumericalJacobianApplyVolume<Acme1MDExperimentalOperatorFullyImplicitDG
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot> >,
        public Dune::PDELab::NumericalJacobianBoundary<Acme1MDExperimentalOperatorFullyImplicitDG
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot> >,
        public Dune::PDELab::NumericalJacobianApplyBoundary<Acme1MDExperimentalOperatorFullyImplicitDG
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

    enum { doPatternVolume = true };
    enum { doPatternSkeleton = true };

    // residual assembly flags
    enum { doAlphaVolume  = true };
    enum { doAlphaSkeleton  = true };
    enum { doAlphaBoundary  = true };
    enum { doLambdaVolume  = true };

    //! constructor: pass parameter object
    Acme1MDExperimentalOperatorFullyImplicitDG (ParamCon& paramCon_, ParamPot& paramPot_, int intorderadd_=0)
      : paramCon(paramCon_),
        paramPot(paramPot_),
        intorderadd(intorderadd_)
    {
      DUNE_THROW(Dune::Exception, "This fully-coupled DG operator is not usable yet, some methods are not implemented!");
    }


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
      // domain and range field type
      typedef typename LFSU::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DF;
      typedef typename LFSU::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType RF;
      typedef typename LFSU::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::JacobianType JacobianType;
      typedef typename LFSU::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeType RangeType;
      typedef typename LFSU::Traits::SizeType size_type;

      bool membrane = param.getPhysics().isMembrane(eg.entity());

      // dimensions
      const int dim = EG::Geometry::dimension;
      const int dimw = EG::Geometry::dimensionworld;
      const int order = lfsuPot.finiteElement().localBasis().order();
      const int intorder = intorderadd + quadrature_factor * order;

      // select quadrature rule
      Dune::GeometryType gt = eg.geometry().type();
      const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

      // evaluate diffusion tensor at cell center, assume it is constant over elements
      typename ParamPot::Traits::PermTensorType APot;
      std::vector<typename ParamCon::Traits::PermTensorType> ACon;
      Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        paramCon.setIonSpecies(j);
        typename ParamCon::Traits::PermTensorType ASingleCon(0.0);
        ASingleCon = paramCon.A(eg.entity(),localcenter);
        ACon.push_back(ASingleCon);
      }

      // transformation
      Dune::FieldMatrix<DF,dimw,dim> jac;

      // loop over quadrature points
      for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        {
          // evaluate basis functions
#if USECACHE==0
          std::vector<RangeType> phi(lfsuPot.size());
          lfsuPot.finiteElement().localBasis().evaluateFunction(it->nodePositions(),phi);
#else
          const std::vector<RangeType>& phi = cache[order].evaluateFunction(it->position(),lfsuPot.finiteElement().localBasis());
#endif

          // evaluate u
          RF uPot=0.0;
          for (size_type i=0; i<lfsuPot.size(); i++)
            uPot += x(lfsuPot,i)*phi[i];

          std::vector<RF> uCon;
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            RF uSingleCon = 0.0;
            for (size_type i=0; i<lfsuCon.child(j).size(); i++)
              uSingleCon += x(lfsuCon.child(j),i)*phi[i];
            uCon.push_back(uSingleCon);
          }

          // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
#if USECACHE==0
          std::vector<JacobianType> js(lfsuPot.size());
          lfsuPot.finiteElement().localBasis().evaluateJacobian(it->nodePositions(),js);
#else
          const std::vector<JacobianType>& js = cache[order].evaluateJacobian(it->position(),lfsuPot.finiteElement().localBasis());
#endif

          // transform gradients of shape functions to real element
          jac = eg.geometry().jacobianInverseTransposed(it->position());
          std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsuPot.size());
          for (size_type i=0; i<lfsuPot.size(); i++)
            jac.mv(js[i][0],gradphi[i]);

          // compute gradient of u
          Dune::FieldVector<RF,dim> graduPot(0.0);
          for (size_type i=0; i<lfsuPot.size(); i++)
            graduPot.axpy(x(lfsuPot,i),gradphi[i]);

          std::vector<Dune::FieldVector<RF,dim> > graduCon;
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            Dune::FieldVector<RF,dim> graduSingleCon(0.0);
            for (size_type i=0; i<lfsuCon.child(j).size(); i++)
              graduSingleCon.axpy(x(lfsuCon.child(j),i),gradphi[i]);
            graduCon.push_back(graduSingleCon);
          }

          // compute K * gradient of u
          Dune::FieldVector<RF,dim> AgraduPot(0.0);
          APot.umv(graduPot,AgraduPot);

          std::vector<Dune::FieldVector<RF,dim> > AgraduCon;
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            Dune::FieldVector<RF,dim> AgraduSingleCon(0.0);
            ACon.umv(graduCon[j],AgraduSingleCon);
            AgraduCon.push_back(AgraduSingleCon);
          }

          // evaluate velocity field
          typename ParamPot::Traits::RangeType bPot = paramPot.b(eg.entity(),it->position());
          // evaluate reaction term
          typename ParamPot::Traits::RangeFieldType cPot = paramPot.c(eg.entity(),it->position());

          std::vector<typename ParamCon::Traits::RangeType> bCon;
          std::vector<typename ParamCon::Traits::RangeFieldType> cCon;
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            paramCon.setIonSpecies(j);
            // evaluate velocity field
            typename ParamCon::Traits::RangeType bSingleCon = paramCon.b(eg.entity(),it->position(), graduPot);
            // evaluate reaction term
            typename ParamCon::Traits::RangeFieldType cSingleCon = paramCon.c(eg.entity(),it->position());

            bCon.push_back(bSingleCon);
            cCon.push_back(cSingleCon);
          }

          // integrate (K grad u - bu)*grad phi_i + a*u*phi_i
          RF factor = it->weight() * eg.geometry().integrationElement(it->position());
          for (size_type i=0; i<lfsvPot.size(); i++)
            r.accumulate(lfsvPot,i,( AgraduPot*gradphi[i] - uPot*(bPot*gradphi[i]) + cPot*uPot*phi[i] )*factor);

          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            for (size_type i=0; i<lfsvCon.child(j).size(); i++)
              r.accumulate(lfsvCon.child(j),i,( AgraduCon[j]*gradphi[i] - uCon[j]*(bCon[j]*gradphi[i])
                  + cCon[j]*uCon[j]*phi[i] )*factor);
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
    */

    /*
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
    */

    /*
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

      // domain and range field type
      typedef typename LFSV::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DF;
      typedef typename LFSV::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType RF;
      typedef typename LFSV::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeType RangeType;
      typedef typename LFSU::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::JacobianType JacobianType;
      typedef typename LFSV::Traits::SizeType size_type;

      // dimensions
      const int dim = IG::dimension;
      const int intorder = intorderadd+2*lfsuPot_s.finiteElement().localBasis().order();

      // evaluate permeability tensors
      const Dune::FieldVector<DF,dim>&
        inside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
      typename ParamPot::Traits::PermTensorType APot_s;
      APot_s = paramPot.A(*(ig.inside()),inside_local);

      std::vector<typename ParamCon::Traits::PermTensorType> ACon_s;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        paramCon.setIonSpecies(j);
        typename ParamPot::Traits::PermTensorType ASingleCon_s(0.0);
        ASinglePot_s = paramCon.A(*(ig.inside()),inside_local);
        ACon_s.push_back(ASinglePot_s);
      }

      // select quadrature rule
      Dune::GeometryType gtface = ig.geometryInInside().type();
      const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

      // face diameter
      DF h_s;
      DF hmax_s = 0.;
      element_size(ig.inside()->geometry(),h_s,hmax_s);
      RF h_F = h_s;
      h_F = ig.inside()->geometry().volume()/ig.geometry().volume(); // Houston!

      // transformation
      Dune::FieldMatrix<DF,dim,dim> jac;

      // evaluate boundary condition
      const Dune::FieldVector<DF,dim-1>
        face_local = Dune::GenericReferenceElements<DF,dim-1>::general(gtface).position(0,0);
      BCType bctypePot = paramPot.bctype(ig.intersection(),face_local);
      bool allNeumann = (bctypePot == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
      bool hasDirichlet = (bctypePot == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet);
      std::vector<BCType> bctypeCon;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        paramCon.setIonSpecies(j);
        BCType bctypeSingleCon = paramCon.bctype(ig.intersection(),face_local);
        bctypeCon.push_back(bctypeSingleCon);
        allNeumann = allNeumann && (bctypeSingleCon == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        hasDirichlet = hasDirichlet || (bctypePot == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet);
      }

      // compute weights
      const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
      Dune::FieldVector<RF,dim> An_FPot_s;
      APot_s.mv(n_F,An_FPot_s);

      std::vector<Dune::FieldVector<RF,dim> > An_FCon_s;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        Dune::FieldVector<RF,dim> An_FSingleCon_s;
        ACon_s[j].mv(n_F,An_FSingleCon_s);
        An_FCon_s.push_back(An_FSingleCon_s);
      }

      RF harmonic_averagePot;
      if (weights==ConvectionDiffusionDGWeights::weightsOn)
        harmonic_averagePot = An_FPot_s*n_F;
      else
        harmonic_averagePot = 1.0;

      std::vector<RF> harmonic_averageCon;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        RF harmonic_averageSingleCon = 1.0;
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          harmonic_averageSingleCon = An_FCon_s[j]*n_F;

        harmonic_averageCon.push_back(harmonic_averageSingleCon);
      }

      // get polynomial degree
      const int order_s = lfsuPot_s.finiteElement().localBasis().order();
      int degree = order_s;

      // penalty factor
      RF penalty_factorPot = (alpha/h_F) * harmonic_averagePot * degree*(degree+dim-1);
      std::vector<RF> penalty_factorCon;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        RF penalty_factorSingleCon = (alpha/h_F) * harmonic_averageCon[j] * degree*(degree+dim-1);
        penalty_factorCon.push_back(penalty_factorSingleCon);
      }

      // loop over quadrature points and integrate normal flux
      for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        {
          // position of quadrature point in local coordinates of elements
          Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());

          // local normal
          const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(it->position());

          // evaluate basis functions
#if USECACHE==0
          std::vector<RangeType> phi_s(lfsuPot_s.size());
          lfsuPot_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);
#else
          const std::vector<RangeType>& phi_s = cache[order_s].evaluateFunction(iplocal_s,lfsuPot_s.finiteElement().localBasis());
#endif

          // integration factor
          RF factor = it->weight() * ig.geometry().integrationElement(it->position());


          // Handle Neumann boundaries
          if (bctypePot == ConvectionDiffusionBoundaryConditions::Neumann)
            {
              // evaluate flux boundary condition
              RF jPot = paramPot.j(ig.intersection(),it->position());

              // integrate
              for (size_type i=0; i<lfsvPot_s.size(); i++)
                r_s.accumulate(lfsuPot_s,i,j * phi_s[i] * factor);
            }
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            paramCon.setIonSpecies(j);
            if (bctypeCon[j] == ConvectionDiffusionBoundaryConditions::Neumann)
              {
                // evaluate flux boundary condition
                RF jCon = paramCon.j(ig.intersection(),it->position());

                // integrate
                for (size_type i=0; i<lfsvCon_s.child(j).size(); i++)
                  r_s.accumulate(lfsvCon_s.child(j),i,jCon * phi_s[i] * factor);
              }
          }


          if(allNeumann) continue;

          // evaluate u
          RF uPot_s=0.0;
          for (size_type i=0; i<lfsuPot_s.size(); i++)
            uPot_s += x_s(lfsuPot_s,i)*phi_s[i];

          std::vector<RF> uCon_s;
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            RF uSingleCon_s = 0.0;
            for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
              uSingleCon_s += x_s(lfsuCon_s.child(j),i)*phi_s[i];
            uCon_s.push_back(uSingleCon_s);
          }

          // transform gradients of shape functions to real element
          jac = ig.inside()->geometry().jacobianInverseTransposed(iplocal_s);
          std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsuPot_s.size());
          for (size_type i=0; i<lfsuPot_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);

          // compute gradient of u
          Dune::FieldVector<RF,dim> graduPot_s(0.0);
          for (size_type i=0; i<lfsuPot_s.size(); i++)
            graduPot_s.axpy(x_s(lfsuPot_s,i),tgradphi_s[i]);

          std::vector<Dune::FieldVector<RF,dim> > graduCon_s;
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            Dune::FieldVector<RF,dim> graduSingleCon_s(0.0);
            for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
              graduSingleCon_s.axpy(x_s(lfsuCon_s.child(j),i),tgradphi_s[i]);
            graduCon_s.push_back(graduSingleCon_s);
          }

          // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
          //typename ParamPot::Traits::RangeType bPot = paramPot.b(*(ig.inside()),iplocal_s);
          typename ParamPot::Traits::RangeType bPot = 0.5*(paramPot.b(*(ig.inside()),iplocal_s)
              +paramPot.b(*(ig.outside()),iplocal_n));

          RF normalfluxPot = bPot*n_F_local;
          if (bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow)
            {
              if (normalflux<-1e-30)
                DUNE_THROW(Dune::Exception,"Outflow boundary condition on inflow!");

              // convection term
              RF term1 = u_s * normalfluxPot *factor;
              for (size_type i=0; i<lfsuPot_s.size(); i++)
                r_s.accumulate(lfsuPot_s,i,term1 * phi_s[i]);

              // evaluate flux boundary condition
              RF oPot = paramPot.o(ig.intersection(),it->position());

              // integrate
              for (size_type i=0; i<lfsvPot_s.size(); i++)
                r_s.accumulate(lfsuPot_s,i,oPot * phi_s[i] * factor);
            }

          std::vector<RF> normalfluxCon;
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            paramCon.setIonSpecies(j);

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            //typename ParamPot::Traits::RangeType bPot = paramPot.b(*(ig.inside()),iplocal_s);
            typename ParamCon::Traits::RangeType bCon = 0.5*(paramCon.b(*(ig.inside()),iplocal_s)
                +paramCon.b(*(ig.outside()),iplocal_n));

            RF normalfluxSingleCon = bCon[j]*n_F_local;
            normalfluxCon.push_back(normalfluxSingleCon);

            if (bctype[j] == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow)
            {
              if (normalfluxCon[j]<-1e-30)
                DUNE_THROW(Dune::Exception,"Outflow boundary condition on inflow!");

              // convection term
              RF term1 = uCon_s[j] * normalfluxCon[j] *factor;
              for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
                r_s.accumulate(lfsuCon_s.child(j),i,term1 * phi_s[i]);

              // evaluate flux boundary condition
              RF oCon = paramCon.o(ig.intersection(),it->position());

              // integrate
              for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
                r_s.accumulate(lfsuCon_s.child(j),i,oCon * phi_s[i] * factor);
            }
          }

          if(not hasDirichlet) continue;

          // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
#if USECACHE==0
          std::vector<JacobianType> gradphi_s(lfsuPot_s.size());
          lfsuPot_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
#else
          const std::vector<JacobianType>& gradphi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsuPot_s.finiteElement().localBasis());
#endif

          if(bctypePot == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet)
            {
              // evaluate Dirichlet boundary condition
              RF gPot = paramPot.g(*(ig.inside()),iplocal_s);

              // upwind
              RF omegaupPot_s, omegaupPot_n;
              if (normalfluxPot>=0.0)
                {
                  omegaupPot_s = 1.0;
                  omegaupPot_n = 0.0;
                }
              else
                {
                  omegaupPot_s = 0.0;
                  omegaupPot_n = 1.0;
                }

              // convection term
              RF term1 = (omegaupPot_s*uPot_s + omegaupPot_n*gPot) * normalfluxPot *factor;
              for (size_type i=0; i<lfsuPot_s.size(); i++)
                r_s.accumulate(lfsuPot_s,i,term1 * phi_s[i]);

              // diffusion term
              RF term2 =  (An_FPot_s*graduPot_s) * factor;
              for (size_type i=0; i<lfsuPot_s.size(); i++)
                r_s.accumulate(lfsuPot_s,i,-term2 * phi_s[i]);

              // (non-)symmetric IP term
              RF term3 = (uPot_s-gPot) * factor;
              for (size_type i=0; i<lfsuPot_s.size(); i++)
                r_s.accumulate(lfsuPot_s,i,term3 * theta * (An_FPot_s*tgradphi_s[i]));

              // standard IP term
              RF term4 = penalty_factorPot * (uPot_s-gPot) * factor;
              for (size_type i=0; i<lfsuPot_s.size(); i++)
                r_s.accumulate(lfsuPot_s,i,term4 * phi_s[i]);

            } else {

              // One of the concentrations has Dirichlet boundary
              for(int j=0; j<NUMBER_OF_SPECIES; j++)
              {
                paramCon.setIonSpecies(j);
                if(bctypeCon[j] == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet)
                {

                RF gCon = paramCon.g(*(ig.inside()),iplocal_s);

                RF omegaupCon_s, omegaupCon_n;
                if (normalfluxCon[j]>=0.0)
                  {
                    omegaupCon_s = 1.0;
                    omegaupCon_n = 0.0;
                  }
                else
                  {
                    omegaupCon_s = 0.0;
                    omegaupCon_n = 1.0;
                  }


                  // convection term
                  RF term1 = (omegaupCon_s*uCon_s[j] + omegaupCon_n*gCon) * normalfluxCon[j] *factor;
                  for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
                    r_s.accumulate(lfsuCon_s.child(j),i,term1 * phi_s[i]);

                  // diffusion term
                  RF term2 =  (An_FCon_s[j]*graduCon_s[j]) * factor;
                  for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
                    r_s.accumulate(lfsuCon_s.child(j),i,-term2 * phi_s[i]);

                  // (non-)symmetric IP term
                  RF term3 = (uCon_s[j]-gCon) * factor;
                  for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
                    r_s.accumulate(lfsuCon_s.child(j),i,term3 * theta * (An_FCon_s[j]*tgradphi_s[i]));

                  // standard IP term
                  RF term4 = penalty_factorCon[j] * (uCon_s[j]-gCon) * factor;
                  for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
                    r_s.accumulate(lfsuCon_s.child(j),i,term4 * phi_s[i]);
                }
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

    ConvectionDiffusionDGMethod::Type method;
    ConvectionDiffusionDGWeights::Type weights;
    Real alpha, beta;
    int intorderadd;
    Real theta;
    typedef typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;

    typedef Dune::PDELab::LocalBasisCache<LocalBasisType> Cache;

    // In theory it is possible that one and the same local operator is
    // called first with a finite element of one type and later with a
    // finite element of another type.  Since finite elements of different
    // type will usually produce different results for the same local
    // coordinate they cannot share a cache.  Here we use a vector of caches
    // to allow for different orders of the shape functions, which should be
    // enough to support p-adaptivity.  (Another likely candidate would be
    // differing geometry types, i.e. hybrid meshes.)

    std::vector<Cache> cache;

    template<class GEO>
    void element_size (const GEO& geo, typename GEO::ctype& hmin, typename GEO::ctype hmax) const
    {
      typedef typename GEO::ctype DF;
      hmin = 1.0E100;
      hmax = -1.0E00;
      const int dim = GEO::coorddimension;
      if (dim==1)
        {
          Dune::FieldVector<DF,dim> x = geo.corner(0);
          x -= geo.corner(1);
          hmin = hmax = x.two_norm();
          return;
        }
      else
        {
          Dune::GeometryType gt = geo.type();
          for (int i=0; i<Dune::GenericReferenceElements<DF,dim>::general(gt).size(dim-1); i++)
            {
              Dune::FieldVector<DF,dim> x = geo.corner(Dune::GenericReferenceElements<DF,dim>::general(gt).subEntity(i,dim-1,0,dim));
              x -= geo.corner(Dune::GenericReferenceElements<DF,dim>::general(gt).subEntity(i,dim-1,1,dim));
              hmin = std::min(hmin,x.two_norm());
              hmax = std::max(hmax,x.two_norm());
            }
          return;
        }
    }

};



#endif /* DUNE_AX1_ACME1MD_EXPERIMENTAL_OPERATOR_FULLY_IMPLICIT_DG_HH */
