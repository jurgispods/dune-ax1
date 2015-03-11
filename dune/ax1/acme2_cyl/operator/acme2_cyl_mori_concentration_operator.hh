#ifndef DUNE_AX1_ACME2CYL_MORI_CONCENTRATION_OPERATOR_HH
#define DUNE_AX1_ACME2CYL_MORI_CONCENTRATION_OPERATOR_HH

#include <dune/pdelab/finiteelement/localbasiscache.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/ax1_lfs_tools.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_geometrywrapper.hh>

/**
 * \brief This local operator is designed for the fully-coupled solution of the Nernst-Planck/Poisson system
 * using standard finite elements
 */
template<typename ParamCon, typename ParamPot, typename FiniteElementMapCon, typename FiniteElementMapPot,
  typename SolutionContainer, bool useMembraneContributions=true>
class Acme2CylMoriConcentrationOperator
      : public Dune::PDELab::NumericalJacobianVolume<Acme2CylMoriConcentrationOperator
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot,SolutionContainer,useMembraneContributions> >,
        public Dune::PDELab::NumericalJacobianApplyVolume<Acme2CylMoriConcentrationOperator
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot,SolutionContainer,useMembraneContributions> >,
        public Dune::PDELab::NumericalJacobianBoundary<Acme2CylMoriConcentrationOperator
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot,SolutionContainer,useMembraneContributions> >,
        public Dune::PDELab::NumericalJacobianApplyBoundary<Acme2CylMoriConcentrationOperator
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot,SolutionContainer,useMembraneContributions> >,
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
    Acme2CylMoriConcentrationOperator (ParamCon& paramCon_, ParamPot& paramPot_, const int lfsPotSize,
        const SolutionContainer& solutionContainer_, int intorderadd_=0)
      : paramCon(paramCon_),
        paramPot(paramPot_),
        intorderadd(intorderadd_),
        solutionContainer(solutionContainer_),
        POTENTIAL_RESIDUAL_SCALE(paramPot.getPhysics().getParams().getPotentialResidualScalingFactor()),
        lfsuPotSize(lfsPotSize),
        doVolumeScaling(paramPot.getPhysics().getParams().doVolumeScaling()),
        volumeScale(lfsuPotSize, 1.0),
        local_position(lfsuPotSize),
        local_volume(lfsuPotSize),
        c(lfsuPotSize),
        uCon(0.0), // has FieldVector type because we want to stick it into flux classes
        uConOld(0.0), // has FieldVector type because we want to stick it into flux classes
        bCon(NUMBER_OF_SPECIES),
        cCon(NUMBER_OF_SPECIES),
        fCon(NUMBER_OF_SPECIES),
        jCon(NUMBER_OF_SPECIES),
        bctypeCon(NUMBER_OF_SPECIES),
        tensorCon(NUMBER_OF_SPECIES),
        graduCon(NUMBER_OF_SPECIES),
        AgraduCon(NUMBER_OF_SPECIES),
        useMori(paramPot.getPhysics().getParams().useMori())
    {
      debug_jochen << "Acme2CylMoriConcentrationOperator(), useMembraneContributions="
                   << useMembraneContributions << std::endl;
    }


    // volume integral depending on test and ansatz functions
    template<typename EG_ORIG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG_ORIG& eg_orig, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {
      //debug_jochen << "Acme2CylMoriConcentrationOperator::alpha_volume @ "<< eg_orig.geometry().center() << std::endl;
      if(!useMembraneContributions && paramPot.getPhysics().isMembrane(eg_orig.entity()))
      {
//          debug_jochen << "Skipping intersection coming from membrane element!" << std::endl;
        return;
      }

      typedef typename Acme2CylGeometrySwitch::ElementGeometrySwitch<EG_ORIG>::type EG;
      const EG& eg(eg_orig);

      typename EG::Geometry geo = eg.geometry();

      // Extract leaf local function spaces from type tree
      typedef LFSU LFSU_CON;
      typedef LFSV LFSV_CON;
      // Get MultiDomainLFS from lfsuPot, do binding and solution vector extraction on the root space,
      // then fetch concentrations power LFS and its scalar children from the root space
      typedef typename LFSU::MultiDomainLFS MultiLFS;
      typedef typename MultiLFS::template Child<1>::Type LFSU_POT;
      typedef LFSU_POT LFSV_POT;
      typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV_CON::template Child<0>::Type LFSV_SINGLE_CON;

      const LFSU_CON& lfsuCon = lfsu;
      const LFSV_CON& lfsvCon = lfsv;

      // TODO Create as class members
      const MultiLFS& multiLFS = lfsu.multiDomainLFS();

      // ========================= Create local coordinate vector for POTENTIAL ===========================
      typedef Dune::PDELab::LFSIndexCache<MultiLFS> MultiLFS_CACHE;
      //multiLFS.bind(eg.entity());
      MultiLFS_CACHE multiLFSCache(multiLFS);
      multiLFSCache.update();

      typedef typename SolutionContainer::U::template ConstLocalView<MultiLFS_CACHE> XView;
      // We use _last_ iteration's potential for evaluating the gradient appearing in the drift term
      XView xPotView(*solutionContainer.getSolutionPotNew());
      X xPot(multiLFS.maxSize());
      xPotView.bind(multiLFSCache);
      xPotView.read(xPot);
      xPotView.unbind();

      const LFSU_POT& lfsuPot = multiLFS.template child<1>();

      // ========================= Create local coordinate vector for CONCENTRATION =======================
      // We use _last_ iteration's concentration for evaluating the concentration appearing in the drift term
      XView xConView(*solutionContainer.getSolutionConNew());
      X xCon(multiLFS.maxSize());
      xConView.bind(multiLFSCache);
      xConView.read(xCon);
      xConView.unbind();
      // Now xPot and xCon contain the _old_ local coefficients from last iteration!

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
        Traits::LocalBasisType::Traits::DomainType DomainType;
      typedef typename LFSU_POT::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeType RangeType;

      typedef typename LFSU_POT::Traits::SizeType size_type;

      // dimensions
      const int dim = EG::Geometry::dimension;
      const int dimw = EG::Geometry::dimensionworld;

      // ============================================================================================
      // Determine Lagrange points for each basis (=test) function
      // => get node position the test functions belong to
      if(doVolumeScaling)
      {
        for (int k=0; k<dim; k++)
        {
          CoordinateEvaluation f(k);
          lfsuPot.finiteElement().localInterpolation().interpolate(f,c);
          for (size_type i=0; i<lfsuPot.size(); i++)
            local_position[i][k] = c[i];
        }
        for (size_type i=0; i<lfsuPot.size(); i++)
        {
          local_volume[i] = paramPot.getPhysics().getNodeVolume(geo.global(local_position[i]));
          volumeScale[i] = 1./ local_volume[i];
        }
      }
      // ============================================================================================

      // select quadrature rule
      Dune::GeometryType gt = geo.type();
      const int intorder = intorderadd + 2*lfsuPot.finiteElement().localBasis().order();
      const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
      
      // evaluate diffusion tensor at cell center, assume it is constant over elements
      Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        typename ParamCon::Traits::PermTensorType tensorSingleCon;
        tensorSingleCon = paramCon[j]->A(eg.entity(),localcenter);
        tensorCon[j] = tensorSingleCon;
      }

      // loop over quadrature points
      for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate basis functions
#if USE_CACHE
        const std::vector<RangeType>& phi = cache.evaluateFunction(it->position(),lfsuPot.finiteElement().localBasis());
#else
        std::vector<RangeType> phi(lfsuPot.size());
        lfsuPot.finiteElement().localBasis().evaluateFunction(it->position(),phi);
#endif


        // evaluate u
        // Poisson part: potential
        // TODO This is actually not needed and only calculated for debug output. Remove this!
        RF uPot=0.0;
        for (size_type i=0; i<lfsuPot.size(); i++)
          uPot += xPot(lfsuPot,i)*phi[i];
        //debug_jochen << "[Acme2CylMoriConcentrationOperator] uPot = " << uPot << std::endl;

        // Nernst-Planck part: concentrations
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          RF uSingleCon=0.0;
          for (size_type i=0; i<lfsuCon.child(j).size(); i++)
            uSingleCon += x(lfsuCon.child(j),i)*phi[i];
          uCon[j] = uSingleCon;
          //debug_jochen << "[Acme2CylMoriConcentrationOperator] uCon[" << j << "] = " << uSingleCon << std::endl;

          RF uSingleConOld=0.0;
          for (size_type i=0; i<lfsuCon.child(j).size(); i++)
            uSingleConOld += xCon(lfsuCon.child(j),i)*phi[i];
          uConOld[j] = uSingleConOld;
          //debug_jochen << "[Acme2CylMoriConcentrationOperator] uConOld[" << j << "] = " << uSingleConOld << std::endl;
        }
        
//        DomainType injectionPosition;
//        injectionPosition[0] = paramPot.getPhysics().getParams().stimulation.get("position_x",
//            (paramPot.getPhysics().getParams().xMax() - 3.0) / 2.0);
//        injectionPosition[1] = paramPot.getPhysics().getParams().stimulation.get("position_y",
//            (paramPot.getPhysics().getParams().yMemb()[0] - 3.0) / 2.0);
//
//        DomainType firstCorner = geo.corner(0);
//        DomainType lastCorner = geo.corner(geo.corners()-1);
//        bool elementContainsInjectionPosition = Tools::lessOrEqualThan(firstCorner, injectionPosition)
//                    && Tools::greaterOrEqualThan(lastCorner, injectionPosition);
//        //if(eg_orig.geometry().corner(0)[0]-1e-6 < 0.0 && std::abs(geo.corner(0)[1]-400.) < 1e-6)
//        if(elementContainsInjectionPosition)
//        {
//          debug_jochen << "Acme2CylMoriConcentrationOperator::alpha_volume @ "<< geo.center() << std::endl;
//          debug_jochen << "[Acme2CylMoriConcentrationOperator] uPot = " << uPot << std::endl;
//          for(int j=0; j<NUMBER_OF_SPECIES; j++)
//          {
//            debug_jochen << "[Acme2CylMoriConcentrationOperator] uCon[" << j << "] = " << uCon[j] << std::endl;
//          }
//          for(int j=0; j<NUMBER_OF_SPECIES; j++)
//          {
//            debug_jochen << "[Acme2CylMoriConcentrationOperator] uConOld[" << j << "] = " << uConOld[j] << std::endl;
//          }
//        }

        // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
#if USE_CACHE
        const std::vector<JacobianType>& js = cache.evaluateJacobian(it->position(),lfsuPot.finiteElement().localBasis());
#else
        std::vector<JacobianType> js(lfsuPot.size());
        lfsuPot.finiteElement().localBasis().evaluateJacobian(it->position(),js);
#endif

        // transform gradients of shape functions to real element
        const Dune::FieldMatrix<DF,dimw,dim> jac = geo.jacobianInverseTransposed(it->position());
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsuPot.size());
        for (size_type i=0; i<lfsuPot.size(); i++)
          jac.mv(js[i][0],gradphi[i]);

        // compute gradient of u
        // operator-split: Use _old_ potential gradient for drift term
        Dune::FieldVector<RF,dim> graduPot(0.0);
        for (size_type i=0; i<lfsuPot.size(); i++)
          graduPot.axpy(xPot(lfsuPot,i),gradphi[i]);

        //debug_jochen << "[PoissonNernstPlanckOperator] graduPot = " << graduPot << std::endl;

        // Nernst-Planck part: concentrations
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          // operator-split: Use _new_ concentration gradient for diffusion term
          Dune::FieldVector<RF,dim> graduSingleCon(0.0);
          for (size_type i=0; i<lfsuCon.child(j).size(); i++)
            graduSingleCon.axpy(x(lfsuCon.child(j),i),gradphi[i]);
          graduCon[j] = graduSingleCon;

          //debug_jochen << "[PoissonNernstPlanckOperator] graduCon[" << j << "] = " << graduSingleCon << std::endl;
        }

        // Nernst-Planck part: concentrations
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          Dune::FieldVector<RF,dim> AgraduSingleCon(0.0);
          tensorCon[j].umv(graduCon[j],AgraduSingleCon);
          AgraduCon[j] = AgraduSingleCon;
        }

        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          // paramCon::b depending on gradPot
          typename ParamCon::Traits::RangeType bSingleCon = paramCon[j]->b(eg.entity(),it->position(),graduPot);
          typename ParamCon::Traits::RangeFieldType cSingleCon = paramCon[j]->c(eg.entity(),it->position());
          typename ParamCon::Traits::RangeFieldType fSingleCon = paramCon[j]->f(eg.entity(),it->position());
          bCon[j] = bSingleCon;
          cCon[j] = cSingleCon;
          fCon[j] = fSingleCon;
        }

        // integrate (A grad u)*grad phi_i - u b*grad phi_i + c*u*phi_i
        RF factor = it->weight() * geo.integrationElement(it->position());
        //debug_jochen << "Factor = " << factor << std::endl;

        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          for (size_type i=0; i<lfsuCon.child(j).size(); i++)
          {

//            //if(true)
//            if(geo.center()[1] < 500.0)
//            {
//              debug_jochen << "[LOP] con @" << lfsuCon.child(j).dofIndex(i)
//                  //<< " -> [" << lfsu.gridFunctionSpace().ordering().mapIndex(lfsuPot.dofIndex(i)) << "] "
//                  << ": " << ((AgraduCon[j]*gradphi[i] - uConOld[j]*(bCon[j]*gradphi[i])
//                      + (cCon[j]*uCon[j]-fCon[j])*phi[i] )*factor*volumeScale[i]) << std::endl;
//
//              debug_jochen << "      AgraduCon[: " << j << "]: "<< AgraduCon[j] << std::endl;
//              debug_jochen << "      uConOld[: " << j << "]: "<< uConOld[j] << std::endl;
//              debug_jochen << "      bCon[: " << j << "]: "<< bCon[j] << std::endl;
//              debug_jochen << "      cCon[: " << j << "]: "<< cCon[j] << std::endl;
//              debug_jochen << "      uCon[: " << j << "]: "<< uCon[j] << std::endl;
//              debug_jochen << "      fCon[: " << j << "]: "<< fCon[j] << std::endl;
//              debug_jochen << "      gradphi[: " << i << "]: " << gradphi[i] << std::endl;
//              debug_jochen << "      phi[: " << i << "]: " << phi[i] << std::endl;
//              debug_jochen << "      volumeScale[: " << i << "]: " << volumeScale[i] << std::endl;
//            }

            // operator-split: Use _old_ concentration for drift term
            r.accumulate(lfsuCon.child(j),i,( AgraduCon[j]*gradphi[i] - uConOld[j]*(bCon[j]*gradphi[i])
               + (cCon[j]*uCon[j]-fCon[j])*phi[i] )*factor*volumeScale[i]);
          }
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
    template<typename IG_ORIG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_boundary (const IG_ORIG& ig_orig,
                         const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                         R& r_s) const
    {
      //debug_jochen << "Acme2CylMoriConcentrationOperator::alpha_boundary @ "<< ig_orig.geometry().center() << std::endl;
      if(!useMembraneContributions && paramPot.getPhysics().isMembraneInterface(ig_orig.intersection())
                  && paramPot.getPhysics().isMembrane(*ig_orig.inside()))
      {
        //debug_jochen << "Skipping intersection coming from membrane element!" << std::endl;
        return;
      }

      typedef typename Acme2CylGeometrySwitch::IntersectionGeometrySwitch<IG_ORIG>::type IG;
      const IG& ig(ig_orig);

      // Extract leaf local function spaces from type tree
      typedef LFSU LFSU_CON;
      typedef LFSV LFSV_CON;
      // Get MultiDomainLFS from lfsuPot, do binding and solution vector extraction on the root space,
      // then fetch concentrations power LFS and its scalar children from the root space
      typedef typename LFSU::MultiDomainLFS MultiLFS;
      typedef typename MultiLFS::template Child<1>::Type LFSU_POT;
      typedef LFSU_POT LFSV_POT;
      typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV_CON::template Child<0>::Type LFSV_SINGLE_CON;

      const LFSU_CON& lfsuCon_s = lfsu_s;
      const LFSV_CON& lfsvCon_s = lfsv_s;
      const LFSU_SINGLE_CON lfsuFirst_s = lfsu_s.child(0);

      // TODO Create as class members
      const MultiLFS& multiLFS = lfsu_s.multiDomainLFS();

      // ========================= Create local coordinate vector for POTENTIAL ===========================
      typedef Dune::PDELab::LFSIndexCache<MultiLFS> MultiLFS_CACHE;
      //multiLFS.bind(eg.entity());
      MultiLFS_CACHE multiLFSCache(multiLFS);
      multiLFSCache.update();

      typedef typename SolutionContainer::U::template ConstLocalView<MultiLFS_CACHE> XView;
      // Use _new_ potential values to calculate membrane (and Mori) fluxes
      XView xPotView(*solutionContainer.getSolutionPotNew());
      X xPot(multiLFS.maxSize());
      xPotView.bind(multiLFSCache);
      xPotView.read(xPot);
      xPotView.unbind();

      // We use _old_ concentration for evaluating the membrane fluxes (Mori fluxes don't care anyway)
      XView xConView(*solutionContainer.getSolutionConOld());
      X xCon(multiLFS.maxSize());
      xConView.bind(multiLFSCache);
      xConView.read(xCon);
      xConView.unbind();

      const LFSU_POT& lfsuPot_s = multiLFS.template child<1>();

      //Extract all types from LFSU_SINGLE_CON, assert local function spaces for concentrations are identical!
      typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DF;
      typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType RF;
      typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::JacobianType JacobianType;
      typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeType RangeType;

      typedef typename LFSU_SINGLE_CON::Traits::SizeType size_type;

      // dimensions
      const int dim = IG::dimension;
      const int dimw = IG::dimensionworld;

      // This operator works for identical finite element spaces for con and pot only!
//      assert(lfsuPot_s.finiteElement().localBasis().order()
//          == lfsuCon_s.child(0).finiteElement().localBasis().order());


      // ============================================================================================
      // Determine Lagrange points for each basis (=test) function
      // => get node position the test functions belong to
      if(doVolumeScaling)
      {
        for (int k=0; k<dim; k++)
        {
          CoordinateEvaluation f(k);
          lfsuFirst_s.finiteElement().localInterpolation().interpolate(f,c);
          for (size_type i=0; i<lfsuFirst_s.size(); i++)
            local_position[i][k] = c[i];
        }
        for (size_type i=0; i<lfsuFirst_s.size(); i++)
        {
          local_volume[i] = paramPot.getPhysics().getNodeVolume(
              (*ig.inside()).geometry().global(local_position[i]) );
          volumeScale[i] = 1./local_volume[i];
        }
      }
      // ============================================================================================

      // select quadrature rule
      Dune::GeometryType gtface = ig.geometryInInside().type();
      const int intorder = intorderadd+2*lfsuFirst_s.finiteElement().localBasis().order();
      const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

      // evaluate boundary condition type
      Dune::FieldVector<DF,dim-1> facecenterlocal = Dune::GenericReferenceElements<DF,dim-1>::general(gtface).position(0,0);

      bool allDirichlet = true;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type bctypeSingleCon
          = paramCon[j]->bctype(ig.intersection(),facecenterlocal);
        allDirichlet = (allDirichlet && bctypeSingleCon==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet);
        bctypeCon[j] = bctypeSingleCon;
      }

      // skip rest if we have only Dirichlet boundarys
      if (allDirichlet) return;

      // loop over quadrature points and integrate normal flux
      for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
          // position of quadrature point in local coordinates of element
          Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

          // evaluate shape functions (assume Galerkin method)
#if USE_CACHE
          const std::vector<RangeType>& phi = cache.evaluateFunction(local,lfsuFirst_s.finiteElement().localBasis());
#else
          std::vector<RangeType> phi(lfsuFirst_s.size());
          lfsuFirst_s.finiteElement().localBasis().evaluateFunction(local,phi);
#endif

          // evaluate potential
          RF uPot=0.0;
          for (size_type i=0; i<lfsuPot_s.size(); i++)
            uPot += xPot(lfsuPot_s,i)*phi[i];

          // evaluate concentrations
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            RF uSingleConOld=0.0;
            for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
              uSingleConOld += xCon(lfsuCon_s.child(j),i)*phi[i];
            uConOld[j] = uSingleConOld;
          }

          //debug_jochen << "Quadrature position " << ig.geometry().global(it->position()) << std::endl;
          //debug_jochen << "Factor: "
          //    << it->weight()*ig.geometry().integrationElement(it->position()) << std::endl;

          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            if (bctypeCon[j]==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann)
            {
              // evaluate flux boundary condition
              // New implicit version: Hand over current value of potential and concentrations to flux classes!
              jCon[j] = paramCon[j]->j(ig.intersection(),it->position(),true,true,uPot,uConOld);

//              debug_jochen << "    j" << ION_NAMES[j] << " @ "
//                  << ig.geometry().global(it->position()) << ": " << jCon[j] << std::endl;

              // integrate j
              RF factor = it->weight()*ig.geometry().integrationElement(it->position());

              for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
              {
//                //if(true)
//                if(ig.geometry().center()[1] < 500.0)
//                {
//                  debug_jochen << "[LOP] con @" << lfsuCon_s.child(j).dofIndex(i)
//                      //<< " -> [" << lfsu.gridFunctionSpace().ordering().mapIndex(lfsuPot.dofIndex(i)) << "] "
//                      << ": " << (jCon[j]*phi[i]*factor*volumeScale[i]) << std::endl;
//
//                  debug_jochen << "      jCon[" << j << "]: "<< jCon[j] << std::endl;
//                  debug_jochen << "      phi[: " << i << "]: " << phi[i] << std::endl;
//                  debug_jochen << "      volumeScale[: " << i << "]: " << volumeScale[i] << std::endl;
//                }

                r_s.accumulate(lfsuCon_s.child(j),i,jCon[j]*phi[i]*factor*volumeScale[i]);
              }
            }
          }

          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            if (bctypeCon[j]==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow)
            {
              // Throw exception for outflow BC, we don't want to use it in our application
              DUNE_THROW(Dune::Exception, "Outflow boundary not implemented!");
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

    template<typename U>
    void preAssembly(const U& u)
    {
      paramCon.preAssembly(u);
    }

    void preAssembly()
    {
      paramCon.preAssembly();
    }


  private:
    ParamCon& paramCon;  // Nernst-Planck parameter class
    ParamPot& paramPot;  // Poisson parameter class
    const SolutionContainer& solutionContainer;
    int intorderadd;
    const double POTENTIAL_RESIDUAL_SCALE;

    typedef typename FiniteElementMapPot::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;
    Dune::PDELab::LocalBasisCache<LocalBasisType> cache;

    const int lfsuPotSize;
    const bool doVolumeScaling;

    mutable std::vector<double> volumeScale;
    mutable std::vector<Dune::FieldVector<Real,dim> > local_position;
    mutable std::vector<double> local_volume;
    mutable std::vector<Real> c;

    mutable typename ParamCon::ParamT::ConRangeType uCon;
    mutable typename ParamCon::ParamT::ConRangeType uConOld;
    mutable std::vector<typename ParamCon::Traits::RangeType> bCon;
    mutable std::vector<typename ParamCon::Traits::RangeFieldType> cCon;
    mutable std::vector<typename ParamCon::Traits::RangeFieldType> fCon;
    mutable std::vector<typename ParamCon::Traits::RangeFieldType> jCon;

    mutable std::vector<Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type> bctypeCon;

    mutable std::vector<typename ParamCon::Traits::PermTensorType> tensorCon;
    mutable std::vector<Dune::FieldVector<Real,dim> > graduCon;
    mutable std::vector<Dune::FieldVector<Real,dim> > AgraduCon;
    const bool useMori;
};



#endif /* DUNE_AX1_ACME2CYL_MORI_CONCENTRATION_OPERATOR_HH */
