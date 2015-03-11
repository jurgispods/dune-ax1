/*
 * acme0_operator_fully_implicit_old.hh
 *
 *  Created on: Nov 24, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME0_OPERATOR_FULLY_IMPLICIT_HH
#define DUNE_AX1_ACME0_OPERATOR_FULLY_IMPLICIT_HH

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/acme0/common/convectiondiffusiondg.hh>

template<typename ParamCon, typename ParamPot, typename FiniteElementMapCon, typename FiniteElementMapPot,
         typename NernstPlanckOperator, typename PoissonOperator>
class Acme0OperatorFullyImplicit
      : public Dune::PDELab::NumericalJacobianApplyVolume<Acme0OperatorFullyImplicit
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot,NernstPlanckOperator,PoissonOperator> >,
        public Dune::PDELab::NumericalJacobianApplySkeleton<Acme0OperatorFullyImplicit
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot,NernstPlanckOperator,PoissonOperator> >,
        public Dune::PDELab::NumericalJacobianApplyBoundary<Acme0OperatorFullyImplicit
        <ParamCon,ParamPot,FiniteElementMapCon,FiniteElementMapPot,NernstPlanckOperator,PoissonOperator> >,
        public Dune::PDELab::FullSkeletonPattern,
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::LocalOperatorDefaultFlags,
        public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename ParamCon::Traits::RangeFieldType>
{
    enum { dim = ParamCon::Traits::GridViewType::dimension };

    typedef NernstPlanckOperator OPCon;
    typedef PoissonOperator      OPPot;

    typedef typename ParamCon::Traits::RangeFieldType Real;

  public:
    // pattern assembly flags
    enum { doPatternVolume = (OPCon::doPatternVolume || OPPot::doPatternVolume)};
    enum { doPatternSkeleton = (OPCon::doPatternSkeleton || OPPot::doPatternSkeleton) };

    // residual assembly flags
    enum { doAlphaVolume  = (OPCon::doAlphaVolume || OPPot::doAlphaVolume) };
    enum { doAlphaSkeleton  = (OPCon::doAlphaSkeleton || OPPot::doAlphaSkeleton) };
    enum { doAlphaBoundary  = (OPCon::doAlphaBoundary || OPPot::doAlphaBoundary) };
    enum { doLambdaVolume  = (OPCon::doLambdaVolume || OPPot::doLambdaVolume) };

    //! constructor: pass parameter object
    Acme0OperatorFullyImplicit (ParamCon& paramCon_, ParamPot& paramPot_,
        NernstPlanckOperator& opCon_, PoissonOperator& opPot_)
      : paramCon(paramCon_),
        paramPot(paramPot_),
        opCon(opCon_),
        opPot(opPot_)
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

      // Nernst-Planck part: concentrations
      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        const LFSU_SINGLE_CON& lfsuSingleCon = lfsuCon.child(j);
        const LFSV_SINGLE_CON& lfsvSingleCon = lfsvCon.child(j);
        paramCon.setIonSpecies(j);
        opCon.alpha_volume(eg, lfsuSingleCon, x, lfsvSingleCon, r);
      }

      // Poisson part: potential
      opPot.alpha_volume(eg, lfsuPot, x, lfsvPot, r);
    }

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

      // Nernst-Planck part: concentrations
      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        const LFSU_SINGLE_CON& lfsuSingleCon_s = lfsuCon_s.child(j);
        const LFSV_SINGLE_CON& lfsvSingleCon_s = lfsvCon_s.child(j);
        paramCon.setIonSpecies(j);
        opCon.alpha_boundary(ig, lfsuSingleCon_s, x_s, lfsvSingleCon_s, r_s);
      }

      // Poisson part: potential
      opPot.alpha_boundary(ig, lfsuPot_s, x_s, lfsvPot_s, r_s);
    }

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

    //! set time in parameter class
    void setTime (double t)
    {
      paramCon.setTime(t);
      paramPot.setTime(t);
    }

  private:
    ParamCon& paramCon;  // Nernst-Planck parameter class
    ParamPot& paramPot;  // Poisson parameter class
    OPCon& opCon;
    OPPot& opPot;
};


#endif /* DUNE_AX1_ACME0_OPERATOR_FULLY_IMPLICIT_HH */
