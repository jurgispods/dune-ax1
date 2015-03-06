/*
 * acme2_cyl_operator_fully_implicit_old.hh
 *
 *  Created on: Nov 24, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_OPERATOR_FULLY_COUPLED_HH
#define DUNE_AX1_ACME2CYL_OPERATOR_FULLY_COUPLED_HH

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
bool useMembraneContributions=true>
class Acme2CylOperatorFullyCoupled
      : public Dune::PDELab::NumericalJacobianVolume<Acme2CylOperatorFullyCoupled
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot,useMembraneContributions> >,
        public Dune::PDELab::NumericalJacobianApplyVolume<Acme2CylOperatorFullyCoupled
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot,useMembraneContributions> >,
        public Dune::PDELab::NumericalJacobianBoundary<Acme2CylOperatorFullyCoupled
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot,useMembraneContributions> >,
        public Dune::PDELab::NumericalJacobianApplyBoundary<Acme2CylOperatorFullyCoupled
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot,useMembraneContributions> >,
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
    Acme2CylOperatorFullyCoupled (ParamCon& paramCon_, ParamPot& paramPot_, const int lfsPotSize, int intorderadd_=0)
      : paramCon(paramCon_),
        paramPot(paramPot_),
        intorderadd(intorderadd_),
        POTENTIAL_RESIDUAL_SCALE(paramPot.getPhysics().getParams().getPotentialResidualScalingFactor()),
        lfsuPotSize(lfsPotSize),
        doVolumeScaling(paramPot.getPhysics().getParams().doVolumeScaling()),
        volumeScale(lfsuPotSize, 1.0),
        local_position(lfsuPotSize),
        local_volume(lfsuPotSize),
        c(lfsuPotSize),
        uCon(0.0),
        bCon(NUMBER_OF_SPECIES),
        cCon(NUMBER_OF_SPECIES),
        fCon(NUMBER_OF_SPECIES),
        bctypeCon(NUMBER_OF_SPECIES),
        tensorCon(NUMBER_OF_SPECIES),
        graduCon(NUMBER_OF_SPECIES),
        AgraduCon(NUMBER_OF_SPECIES),
        fullyImplicitMembraneFlux(paramPot.getPhysics().getParams().boundary.get("fullyImplicitMembraneFlux",false)),
        useMori(paramPot.getPhysics().getParams().useMori()),
        useJacobianVolume(paramPot.getPhysics().getParams().general.get("useElecOperatorJacobianVolume",false)),
        lengthScale(paramPot.getPhysics().getLengthScale()),
        timeScale(paramPot.getPhysics().getTimeScale())
    {
      debug_jochen << "Acme2CylOperatorFullyCoupled(), useMembraneContributions="
                   << useMembraneContributions << std::endl;

      assert(! useMori);
    }


    // volume integral depending on test and ansatz functions
    template<typename EG_ORIG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG_ORIG& eg_orig, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {
      //debug_jochen << "PoissonNernstPlanckOperator::alpha_volume @ "<< eg_orig.geometry().center() << std::endl;
//      if(!useMembraneContributions && paramPot.getPhysics().isMembrane(eg_orig.entity()))
//      {
//        //debug_jochen << "Skipping intersection coming from membrane element!" << std::endl;
//        return;
//      }

      typedef typename Acme2CylGeometrySwitch::ElementGeometrySwitch<EG_ORIG>::type EG;
      const EG& eg(eg_orig);

      // Extract leaf local function spaces from type tree
      typedef typename LFSU::template Child<0>::Type LFSU_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_CON;
      typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV_CON::template Child<0>::Type LFSV_SINGLE_CON;
      typedef typename LFSU::template Child<1>::Type LFSU_POT;
      typedef typename LFSV::template Child<1>::Type LFSV_POT;

      const LFSU_CON& lfsuCon = lfsu.template child<0>();
      const LFSV_CON& lfsvCon = lfsv.template child<0>();
      const LFSU_POT& lfsuPot = lfsu.template child<1>();
      const LFSV_POT& lfsvPot = lfsv.template child<1>();

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
          local_volume[i] = paramPot.getPhysics().getNodeVolume(eg.geometry().global(local_position[i]));
          volumeScale[i] = 1./ local_volume[i];
        }
      }
      // ============================================================================================

      // select quadrature rule
      Dune::GeometryType gt = eg.geometry().type();
      const int intorder = intorderadd + 2*lfsuPot.finiteElement().localBasis().order();
      const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
      
      // evaluate diffusion tensor at cell center, assume it is constant over elements
      Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
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
        RF uPot=0.0;
        for (size_type i=0; i<lfsuPot.size(); i++)
          uPot += x(lfsuPot,i)*phi[i];

        //debug_jochen << "[PoissonNernstPlanckOperator] uPot = " << uPot << std::endl;

        //debug_jochen << "uPot = " << uPot << std::endl;

        // Nernst-Planck part: concentrations
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          RF uSingleCon=0.0;
          for (size_type i=0; i<lfsuCon.child(j).size(); i++)
            uSingleCon += x(lfsuCon.child(j),i)*phi[i];
          uCon[j] = uSingleCon;

          //debug_jochen << "[PoissonNernstPlanckOperator] uCon[" << j << "] = " << uSingleCon << std::endl;
        }
        
        // evaluate diffusion tensor at quadrature point (Mori: not constant over elements!)
        typename ParamPot::Traits::PermTensorType tensorPot;
        // Mori parameters need concentrations
        tensorPot = paramPot.A(eg.entity(),it->position(),uCon);
        //tensorPot = paramPot.A(eg.entity(),it->position());

        // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
#if USE_CACHE
        const std::vector<JacobianType>& js = cache.evaluateJacobian(it->position(),lfsuPot.finiteElement().localBasis());
#else
        std::vector<JacobianType> js(lfsuPot.size());
        lfsuPot.finiteElement().localBasis().evaluateJacobian(it->position(),js);
#endif


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

        //debug_jochen << "[PoissonNernstPlanckOperator] graduPot = " << graduPot << std::endl;

        // Nernst-Planck part: concentrations
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          Dune::FieldVector<RF,dim> graduSingleCon(0.0);
          for (size_type i=0; i<lfsuCon.child(j).size(); i++)
            graduSingleCon.axpy(x(lfsuCon.child(j),i),gradphi[i]);
          graduCon[j] = graduSingleCon;

          //debug_jochen << "[PoissonNernstPlanckOperator] graduCon[" << j << "] = " << graduSingleCon << std::endl;
        }

        // compute A * gradient of u
        // Poisson part: potential
        Dune::FieldVector<RF,dim> AgraduPot(0.0);
        tensorPot.umv(graduPot,AgraduPot);

        // Nernst-Planck part: concentrations
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          Dune::FieldVector<RF,dim> AgraduSingleCon(0.0);
          tensorCon[j].umv(graduCon[j],AgraduSingleCon);
          AgraduCon[j] = AgraduSingleCon;
        }

        // evaluate velocity field, sink term and source term
        typename ParamPot::Traits::RangeType bPot = paramPot.b(eg.entity(),it->position());
        typename ParamPot::Traits::RangeFieldType cPot = paramPot.c(eg.entity(),it->position());
        typename ParamPot::Traits::RangeFieldType fPot = paramPot.f(eg.entity(),it->position(),uCon); // depending on cd

        typename ParamCon::Traits::RangeFieldType fConSum(0.0);
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          // paramCon::b depending on gradPot
          typename ParamCon::Traits::RangeType bSingleCon = paramCon[j]->b(eg.entity(),it->position(),graduPot);
          typename ParamCon::Traits::RangeFieldType cSingleCon = paramCon[j]->c(eg.entity(),it->position());
          typename ParamCon::Traits::RangeFieldType fSingleCon = paramCon[j]->f(eg.entity(),it->position());
          bCon[j] = bSingleCon;
          cCon[j] = cSingleCon;
          fCon[j] = fSingleCon;

          fConSum += paramPot.getPhysics().getValence(j) * fCon[j];
        }

        // integrate (A grad u)*grad phi_i - u b*grad phi_i + c*u*phi_i
        RF factor = it->weight() * eg.geometry().integrationElement(it->position());

        //debug_jochen << "Factor = " << factor << std::endl;


//        if(fCon[Na] > 0)
//        {
//          debug_jochen << eg.geometry().center() << ": fCon[Na] = " << fCon[Na]
//            << "<-> AgraduCon[Na] = " << AgraduCon[Na] << " <-> bCon[Na] = " << bCon[Na] << std::endl;
//          double termA = 0.0;
//          double termB = 0.0;
//          double termF = 0.0;
//          for (size_type i=0; i<lfsuCon.child(Na).size(); i++)
//          {
//            termA += AgraduCon[Na]*gradphi[i];
//            termB += - uCon[Na]*(bCon[Na]*gradphi[i]);
//            termF += -fCon[Na]*phi[i];
//          }
//          debug_jochen << "   termF = " << termF  << " <-> termA = " << termA << " <-> termB = "
//              << termB << std::endl;
//
//          double termAPot = 0.0;
//          double termBPot = 0.0;
//          double termFPot = 0.0;
//          for (size_type i=0; i<lfsuCon.child(0).size(); i++)
//          {
//            termAPot += AgraduPot*gradphi[i];
//            termBPot += bGrad*gradphi[i];
//            termFPot += -fPot*phi[i];
//          }
//          debug_jochen << "   termFPot = " << termFPot  << " <-> termAPot = " << termAPot << " <-> termBPot = "
//              << termBPot << std::endl;
//        }


        //debug_jochen << "QuadratureRule it->weight = " << it->weight()
        //    << ", geometry().integrationElement(" << eg.geometry().global(it->position()) << ") = "
        //    << it->weight() * eg.geometry().integrationElement(it->position()) << std::endl;

        for (size_type i=0; i<lfsuPot.size(); i++)
        {
//          debug_jochen << "[LOP] pot @" << lfsuPot.dofIndex(i) << " -> [" << lfsu.gridFunctionSpace().ordering().mapIndex(lfsuPot.dofIndex(i)) << "] "
//              << (( AgraduPot*gradphi[i] - uPot*(bPot*gradphi[i]) + (cPot*uPot-fPot)*phi[i] )
//              *factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE) << std::endl;
          r.accumulate(lfsuPot,i,( AgraduPot*gradphi[i] - uPot*(bPot*gradphi[i]) + (cPot*uPot-fPot)*phi[i] )
              *factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE);
        }

        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          for (size_type i=0; i<lfsuCon.child(j).size(); i++)
          {
//            debug_jochen << "[LOP] con @" << lfsuCon.child(j).dofIndex(i)
//                << " -> [" << lfsu.gridFunctionSpace().ordering().mapIndex(lfsuCon.child(j).dofIndex(i)) << "] "
//                << (( AgraduCon[j]*gradphi[i] - uCon[j]*(bCon[j]*gradphi[i]) + (cCon[j]*uCon[j]-fCon[j])*phi[i] )*factor*volumeScale[i])
//                << std::endl;

//            if(fCon[j] > 1e-100)
//            {
//              debug_jochen << "LOP1: " << ( AgraduCon[j]*gradphi[i]) << std::endl;
//              debug_jochen << "LOP2: " << (uCon[j]*(bCon[j]*gradphi[i])) << std::endl;
//              debug_jochen << "LOP3: " << ((cCon[j]*uCon[j]-fCon[j])*phi[i] ) << std::endl;
//              debug_jochen << "LOP4: " << (uCon[j] ) << std::endl;
//              debug_jochen << "LOP5: " << (fCon[j] ) << std::endl;
//            }

            r.accumulate(lfsuCon.child(j),i,( AgraduCon[j]*gradphi[i] - uCon[j]*(bCon[j]*gradphi[i])
               + (cCon[j]*uCon[j]-fCon[j])*phi[i] )*factor*volumeScale[i]);
          }
        }

      }
    }

#ifdef AX1_USE_JACOBIAN_METHODS_IN_OPERATOR
    // jacobian of volume term
    template<typename EG_ORIG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_volume (const EG_ORIG& eg_orig, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          M& mat) const
    {
      if(! useJacobianVolume)
      {
        Dune::PDELab::NumericalJacobianVolume<Acme2CylOperatorFullyCoupled<ParamCon,ParamPot,FiniteElementMapCon,
                FiniteElementMapPot,useMembraneContributions> >::jacobian_volume(eg_orig, lfsu, x, lfsv,mat);
        return;
      }

      typedef typename Acme2CylGeometrySwitch::ElementGeometrySwitch<EG_ORIG>::type EG;
      const EG& eg(eg_orig);

      // Extract leaf local function spaces from type tree
      typedef typename LFSU::template Child<0>::Type LFSU_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_CON;
      typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV_CON::template Child<0>::Type LFSV_SINGLE_CON;
      typedef typename LFSU::template Child<1>::Type LFSU_POT;
      typedef typename LFSV::template Child<1>::Type LFSV_POT;

      const LFSU_CON& lfsuCon = lfsu.template child<0>();
      const LFSV_CON& lfsvCon = lfsv.template child<0>();
      const LFSU_POT& lfsuPot = lfsu.template child<1>();
      const LFSV_POT& lfsvPot = lfsv.template child<1>();

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
          local_volume[i] = paramPot.getPhysics().getNodeVolume(eg.geometry().global(local_position[i]));
          volumeScale[i] = 1./ local_volume[i];
        }
      }
      // ============================================================================================

      // select quadrature rule
      Dune::GeometryType gt = eg.geometry().type();
      const int intorder = intorderadd + 2*lfsuPot.finiteElement().localBasis().order();
      const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

      // evaluate diffusion tensor at cell center, assume it is constant over elements
      Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
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

//        // evaluate u
//        // Poisson part: potential
//        RF uPot=0.0;
//        for (size_type i=0; i<lfsuPot.size(); i++)
//          uPot += x(lfsuPot,i)*phi[i];

        // Nernst-Planck part: concentrations
        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          RF uSingleCon=0.0;
          for (size_type i=0; i<lfsuCon.child(j).size(); i++)
            uSingleCon += x(lfsuCon.child(j),i)*phi[i];
          uCon[j] = uSingleCon;

          //debug_jochen << "[PoissonNernstPlanckOperator] uCon[" << j << "] = " << uSingleCon << std::endl;
        }

        // evaluate diffusion tensor at quadrature point (Mori: not constant over elements!)
        typename ParamPot::Traits::PermTensorType tensorPot;
        // Mori parameters need concentrations
        tensorPot = paramPot.A(eg.entity(),it->position(),uCon);

        // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
#if USE_CACHE
        const std::vector<JacobianType>& js = cache.evaluateJacobian(it->position(),lfsuPot.finiteElement().localBasis());
#else
        std::vector<JacobianType> js(lfsuPot.size());
        lfsuPot.finiteElement().localBasis().evaluateJacobian(it->position(),js);
#endif

        // transform gradients of shape functions to real element
        const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsuPot.size());
        Dune::FieldVector<RF,dim> graduPot(0.0);
        std::vector<Dune::FieldVector<RF,dim> > AgradphiPot(lfsuPot.size());
        for (size_type i=0; i<lfsuPot.size(); i++)
        {
          jac.mv(js[i][0],gradphi[i]);
          // Calculate potential gradient (needed for Acon)
          graduPot.axpy(x(lfsuPot,i),gradphi[i]);
          tensorPot.mv(gradphi[i],AgradphiPot[i]);
        }

        std::vector<std::vector<Dune::FieldVector<RF,dim> > > AgradphiCon(NUMBER_OF_SPECIES);
        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          AgradphiCon[j].resize(lfsuPot.size());
          for (size_type i=0; i<lfsuPot.size(); i++)
          {
            tensorCon[j].mv(gradphi[i],AgradphiCon[j][i]);
          }
        }

        // evaluate velocity field, sink term and source term (potential)
        typename ParamPot::Traits::RangeType bPot = paramPot.b(eg.entity(),it->position());
        typename ParamPot::Traits::RangeFieldType cPot = paramPot.c(eg.entity(),it->position());

        // evaluate velocity field, sink term and source term (concentrations)
        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          bCon[j] = paramCon[j]->b(eg.entity(),it->position(),graduPot);
          cCon[j] = paramCon[j]->c(eg.entity(),it->position());
        }

        RF factor = it->weight() * eg.geometry().integrationElement(it->position());

        // integrate (A grad phi_j)*grad phi_i - phi_j b*grad phi_i + c*phi_j*phi_i (concentrations)
        for(int k=0; k<NUMBER_OF_SPECIES; ++k)
        {
          // DOF j is the changing variable
          for (size_type j=0; j<lfsuCon.child(k).size(); j++)
          {
            // concentration-concentration coupling: influence of j on concentration DOF i
            for (size_type i=0; i<lfsuCon.child(k).size(); i++)
            {
              mat.accumulate(lfsuCon.child(k),i,lfsuCon.child(k),j,( AgradphiCon[k][j]*gradphi[i]-phi[j]*(bCon[k]*gradphi[i])
                  +cCon[k]*phi[j]*phi[i] )*factor*volumeScale[i]);
            }
            // concentration-potential coupling: influence of j on potential DOF i
            for (size_type i=0; i<lfsuPot.size(); i++)
            {
              // For the influence of a change in n_j on the potenial, we need the derivative of the
              // concentration-dependent right hand side function f() of the Poisson equation
              typename ParamPot::Traits::RangeFieldType dfPot_du =
                  paramPot.dfdu(eg.entity(),it->position(),k,phi[j]);

              mat.accumulate(lfsuPot,i,lfsuCon.child(k),j,( -dfPot_du*phi[i] )
                  *factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE);
            }
          }
        }

        // integrate (A grad phi_j)*grad phi_i - phi_j b*grad phi_i + c*phi_j*phi_i (potential)
        for (size_type j=0; j<lfsuPot.size(); j++)
        {
          // potential-potential coupling
          for (size_type i=0; i<lfsuPot.size(); i++)
          {
            mat.accumulate(lfsuPot,i,lfsuPot,j,( AgradphiPot[j]*gradphi[i]-phi[j]*(bPot*gradphi[i])+cPot*phi[j]*phi[i] )
                *factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE);
          }
          // potential-concentration couplings
          for(int k=0; k<NUMBER_OF_SPECIES; ++k)
          {
            for (size_type i=0; i<lfsuCon.child(k).size(); i++)
            {
              typename ParamCon::Traits::RangeType dbCon_du =
                    paramCon[k]->dbdu(eg.entity(),it->position(),gradphi[j]);

              mat.accumulate(lfsuCon.child(k),i,lfsuPot,j,( -uCon[k]*(dbCon_du*gradphi[i]) )*factor*volumeScale[i]);
            }
          }
        }
      }
    }
#endif

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
      //debug_jochen << "PoissonNernstPlanckOperator::alpha_boundary @ "<< ig_orig.geometry().center() << std::endl;
      if(!useMembraneContributions && paramPot.getPhysics().isMembraneInterface(ig_orig.intersection())
                  && paramPot.getPhysics().isMembrane(*ig_orig.inside()))
      {
        //debug_jochen << "Skipping intersection coming from membrane element!" << std::endl;
        return;
      }

      typedef typename Acme2CylGeometrySwitch::IntersectionGeometrySwitch<IG_ORIG>::type IG;
      const IG& ig(ig_orig);

      // Extract leaf local function spaces from type tree
      typedef typename LFSU::template Child<0>::Type LFSU_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_CON;
      typedef typename LFSU_CON::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV_CON::template Child<0>::Type LFSV_SINGLE_CON;
      typedef typename LFSU::template Child<1>::Type LFSU_POT;
      typedef typename LFSV::template Child<1>::Type LFSV_POT;

      const LFSU_CON& lfsuCon_s = lfsu_s.template child<0>();
      const LFSV_CON& lfsvCon_s = lfsv_s.template child<0>();
      const LFSU_POT& lfsuPot_s = lfsu_s.template child<1>();
      const LFSV_POT& lfsvPot_s = lfsv_s.template child<1>();

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


      // ============================================================================================
      // Determine Lagrange points for each basis (=test) function
      // => get node position the test functions belong to
      if(doVolumeScaling)
      {
        for (int k=0; k<dim; k++)
        {
          CoordinateEvaluation f(k);
          lfsuPot_s.finiteElement().localInterpolation().interpolate(f,c);
          for (size_type i=0; i<lfsuPot_s.size(); i++)
            local_position[i][k] = c[i];
        }
        for (size_type i=0; i<lfsuPot_s.size(); i++)
        {
          local_volume[i] = paramPot.getPhysics().getNodeVolume(
              (*ig.inside()).geometry().global(local_position[i]) );
          volumeScale[i] = 1./local_volume[i];
        }
      }
      // ============================================================================================

      // select quadrature rule
      Dune::GeometryType gtface = ig.geometryInInside().type();
      const int intorder = intorderadd+2*lfsuPot_s.finiteElement().localBasis().order();
      const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

      // evaluate boundary condition type
      Dune::FieldVector<DF,dim-1> facecenterlocal = Dune::ReferenceElements<DF,dim-1>::general(gtface).position(0,0);
      Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type bctypePot;
      bctypePot = paramPot.bctype(ig.intersection(),facecenterlocal);

      bool allDirichlet = (bctypePot==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet);
      bool allConcentrationsNeumann = true;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type bctypeSingleCon
          = paramCon[j]->bctype(ig.intersection(),facecenterlocal);
        allDirichlet = (allDirichlet && bctypeSingleCon==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet);
        allConcentrationsNeumann = allConcentrationsNeumann &&
            (bctypeSingleCon == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
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
          const std::vector<RangeType>& phi = cache.evaluateFunction(local,lfsuPot_s.finiteElement().localBasis());
#else
          std::vector<RangeType> phi(lfsuPot_s.size());
          lfsuPot_s.finiteElement().localBasis().evaluateFunction(local,phi);
#endif

          // evaluate potential
          RF uPot=0.0;
          for (size_type i=0; i<lfsuPot_s.size(); i++)
            uPot += x_s(lfsuPot_s,i)*phi[i];

          // evaluate concentrations
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            RF uSingleCon=0.0;
            for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
              uSingleCon += x_s(lfsuCon_s.child(j),i)*phi[i];
            uCon[j] = uSingleCon;
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
              typename ParamCon::Traits::RangeFieldType jSingleCon =
                  paramCon[j]->j(ig.intersection(),it->position(),false,fullyImplicitMembraneFlux,uPot,uCon);

//              debug_jochen << "    j" << ION_NAMES[j] << " @ "
//                  << ig.geometry().global(it->position()) << ": " << jSingleCon << std::endl;

              // integrate j
              RF factor = it->weight()*ig.geometry().integrationElement(it->position());

              for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
              {
//                debug_jochen << "[LOP] con @" << lfsuCon_s.child(j).dofIndex(i)
//                                << " -> [" << lfsu_s.gridFunctionSpace().ordering().mapIndex(lfsuCon_s.child(j).dofIndex(i)) << "] "
//                                << (jSingleCon*phi[i]*factor*volumeScale[i]) << std::endl;
                r_s.accumulate(lfsuCon_s.child(j),i,jSingleCon*phi[i]*factor*volumeScale[i]);
              }
            }
          }

          // Only do something if there is a Neumann boundary for the potential; do nothing when bctype
          // is 'Dirichlet' or 'None'
          if (bctypePot==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann)
          {
            // evaluate flux boundary condition
            typename ParamPot::Traits::RangeFieldType jPot(0.0);

            jPot = paramPot.j(ig.intersection(),it->position());

            // integrate j
            RF factor = it->weight()*ig.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsuPot_s.size(); i++)
            {
//                debug_jochen << "[LOP] pot @" << lfsuPot_s.dofIndex(i) << " -> [" << lfsu_s.gridFunctionSpace().ordering().mapIndex(lfsuPot_s.dofIndex(i))
//                    << "] " << (jPot*phi[i]*factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE) << std::endl;
              r_s.accumulate(lfsuPot_s,i,jPot*phi[i]*factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE);
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
            {
              r_s.accumulate(lfsuPot_s,i,( (bPot*n)*uPot + oPot)*phi[i]*factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE);
            }
          }
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            if (bctypeCon[j]==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow)
            {
              // Throw exception for outflow BC, we don't want to use it in our application
              DUNE_THROW(Dune::Exception, "Outflow boundary not implemented!");

              // evaluate u
              RF uCon=0.0;
              for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
                uCon += x_s(lfsuCon_s.child(j),i)*phi[i];

              // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
#if USE_CACHE
              const std::vector<JacobianType>& js = cache.evaluateJacobian(local,lfsuPot_s.finiteElement().localBasis());
#else
              std::vector<JacobianType> js(lfsuPot_s.size());
              lfsuPot_s.finiteElement().localBasis().evaluateJacobian(local,js);
#endif


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
              typename ParamCon::Traits::RangeType bCon = paramCon[j]->b(*(ig.inside()),local,graduPot);
              const Dune::FieldVector<DF,dim> n = ig.unitOuterNormal(it->position());

              // evaluate outflow boundary condition
              typename ParamCon::Traits::RangeFieldType oCon = paramCon[j]->o(ig.intersection(),it->position());

              // integrate o
              RF factor = it->weight()*ig.geometry().integrationElement(it->position());
              for (size_type i=0; i<lfsuCon_s.child(j).size(); i++)
              {
                r_s.accumulate(lfsuCon_s.child(j),i,( (bCon*n)*uCon + oCon)*phi[i]*factor*volumeScale[i]);
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

    void preAssembly()
    {
      paramCon.preAssembly(fullyImplicitMembraneFlux);
    }

  private:
    ParamCon& paramCon;  // Nernst-Planck parameter class
    ParamPot& paramPot;  // Poisson parameter class
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
    mutable std::vector<typename ParamCon::Traits::RangeType> bCon;
    mutable std::vector<typename ParamCon::Traits::RangeFieldType> cCon;
    mutable std::vector<typename ParamCon::Traits::RangeFieldType> fCon;

    mutable std::vector<Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type> bctypeCon;

    mutable std::vector<typename ParamCon::Traits::PermTensorType> tensorCon;
    mutable std::vector<Dune::FieldVector<Real,dim> > graduCon;
    mutable std::vector<Dune::FieldVector<Real,dim> > AgraduCon;
    const bool fullyImplicitMembraneFlux;
    const bool useMori;
    const bool useJacobianVolume;

    double lengthScale;
    double timeScale;

};



#endif /* DUNE_AX1_ACME2CYL_OPERATOR_FULLY_COUPLED_HH */
