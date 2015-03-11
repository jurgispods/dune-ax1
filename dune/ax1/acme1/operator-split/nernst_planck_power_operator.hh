/*
 * nernst_planck_power_operator.hh
 *
 *  Created on: Sep 2, 2011
 *      Author: jpods
 */

#ifndef NERNST_PLANCK_POWER_OPERATOR_HH
#define NERNST_PLANCK_POWER_OPERATOR_HH

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/pattern.hh>

#include <dune/ax1/common/constants.hh>

template<typename T, typename FiniteElementMap, typename NernstPlanckOperator>
class NernstPlanckDGPowerOperator
      : public Dune::PDELab::NumericalJacobianApplyVolume<
          NernstPlanckDGPowerOperator<T,FiniteElementMap,NernstPlanckOperator> >,
        //public Dune::PDELab::NumericalJacobianApplySkeleton<
        //  NernstPlanckDGPowerOperator<T,FiniteElementMap,NernstPlanckOperator> >,
        public Dune::PDELab::NumericalJacobianApplyBoundary<
          NernstPlanckDGPowerOperator<T,FiniteElementMap,NernstPlanckOperator> >,
        //public Dune::PDELab::FullSkeletonPattern,
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::LocalOperatorDefaultFlags,
        public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>
{
    enum { dim = T::Traits::GridViewType::dimension };

    typedef NernstPlanckOperator OP;
    typedef typename T::Traits::RangeFieldType Real;

  public:
    // pattern assembly flags
    enum { doPatternVolume = OP::doPatternVolume };
    enum { doPatternSkeleton = OP::doPatternSkeleton };

    // residual assembly flags
    enum { doAlphaVolume  = OP::doAlphaVolume };
    enum { doAlphaSkeleton  = OP::doAlphaSkeleton };
    enum { doAlphaBoundary  = OP::doAlphaBoundary };
    enum { doLambdaVolume  = OP::doLambdaVolume };

    //! constructor: pass parameter object
    NernstPlanckDGPowerOperator (T& param_, OP& op_)
      : param(param_), op(op_)
    {}


    // volume integral depending on test and ansatz functions
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {
      // local function space for the NUMBER_OF_SPECIES components (assume Galerkin scheme U=V)
      typedef typename LFSU::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_SINGLE_CON;

      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        const LFSU_SINGLE_CON& lfsuSingleCon = lfsu.child(j);
        const LFSV_SINGLE_CON& lfsvSingleCon = lfsv.child(j);

        param.setIonSpecies(j);

        op.alpha_volume(eg, lfsuSingleCon, x, lfsvSingleCon, r);
      }
    }

    // jacobian of volume term
    template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          M& mat) const
    {
      // local function space for the NUMBER_OF_SPECIES components (assume Galerkin scheme U=V)
      typedef typename LFSU::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_SINGLE_CON;

      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        const LFSU_SINGLE_CON& lfsuSingleCon = lfsu.child(j);
        const LFSV_SINGLE_CON& lfsvSingleCon = lfsv.child(j);

        param.setIonSpecies(j);

        op.jacobian_volume(eg, lfsuSingleCon, x, lfsvSingleCon, mat);
      }
    }

    // skeleton integral depending on test and ansatz functions
    // each face is only visited ONCE!
    template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_skeleton (const IG& ig,
                         const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                         const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                         R& r_s, R& r_n) const
    {
      // local function space for the NUMBER_OF_SPECIES components (assume Galerkin scheme U=V)
      typedef typename LFSU::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_SINGLE_CON;

      // Check if we are on the membrane boundary
      bool isMembrane_s = param.getPhysics().isMembrane(*ig.inside());
      bool isMembrane_n = param.getPhysics().isMembrane(*ig.outside());
      bool membrane = (isMembrane_s || isMembrane_n);
      bool membraneBoundary = (isMembrane_s || isMembrane_n) && not (isMembrane_s && isMembrane_n);

      // We are on the membrane boundary
      if(membraneBoundary)
      {
        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          const LFSU_SINGLE_CON& lfsuSingleCon_s = lfsu_s.child(j);
          const LFSV_SINGLE_CON& lfsvSingleCon_s = lfsv_s.child(j);
          const LFSU_SINGLE_CON& lfsuSingleCon_n = lfsu_n.child(j);
          const LFSV_SINGLE_CON& lfsvSingleCon_n = lfsv_n.child(j);

          param.setIonSpecies(j);

          // Treat the membrane as a boundary
          op.alpha_boundary(ig, lfsuSingleCon_s, x_s, lfsvSingleCon_s, r_s);
          op.alpha_boundary(ig, lfsuSingleCon_n, x_n, lfsvSingleCon_n, r_n);
        }
      }

      /*
       * Legacy code used with ConvectionDiffusionDG operator
       *
      // No flux on the membrane boundary or within membrane => skip alpha_skeleton!
      if(not membrane)
      {
        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          const LFSU_SINGLE_CON& lfsuSingleCon_s = lfsu_s.child(j);
          const LFSV_SINGLE_CON& lfsvSingleCon_s = lfsv_s.child(j);
          const LFSU_SINGLE_CON& lfsuSingleCon_n = lfsu_n.child(j);
          const LFSV_SINGLE_CON& lfsvSingleCon_n = lfsv_n.child(j);

          param.setIonSpecies(j);

          op.alpha_skeleton(ig, lfsuSingleCon_s, x_s, lfsvSingleCon_s,
                                lfsuSingleCon_n, x_n, lfsvSingleCon_n, r_s, r_n);
        }
      }
      */
    }

    template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_skeleton (const IG& ig,
                            const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                            const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                            M& mat_ss, M& mat_sn,
                            M& mat_ns, M& mat_nn) const
    {
      // local function space for the NUMBER_OF_SPECIES components (assume Galerkin scheme U=V)
      typedef typename LFSU::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_SINGLE_CON;

      // Check if we are on the membrane boundary
      bool isMembrane_s = param.getPhysics().isMembrane(*ig.inside());
      bool isMembrane_n = param.getPhysics().isMembrane(*ig.outside());
      bool membrane = (isMembrane_s || isMembrane_n);
      bool membraneBoundary = (isMembrane_s || isMembrane_n) && not (isMembrane_s && isMembrane_n);

      // We are on the membrane boundary
      if(membraneBoundary)
      {
        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          const LFSU_SINGLE_CON& lfsuSingleCon_s = lfsu_s.child(j);
          const LFSV_SINGLE_CON& lfsvSingleCon_s = lfsv_s.child(j);
          const LFSU_SINGLE_CON& lfsuSingleCon_n = lfsu_n.child(j);
          const LFSV_SINGLE_CON& lfsvSingleCon_n = lfsv_n.child(j);

          param.setIonSpecies(j);

          // Treat the membrane as a boundary
          op.jacobian_boundary(ig, lfsuSingleCon_s, x_s, lfsvSingleCon_s, mat_ss); // TODO Is this correct?
          op.jacobian_boundary(ig, lfsuSingleCon_n, x_n, lfsvSingleCon_n, mat_nn); // TODO Is this correct?
        }
      }

      /*
       * Legacy code used with ConvectionDiffusionDG operator
       *
      // No flux on the membrane boundary or within membrane => skip jacobian_skeleton!
      if(not membrane)
      {
        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          const LFSU_SINGLE_CON& lfsuSingleCon_s = lfsu_s.child(j);
          const LFSV_SINGLE_CON& lfsvSingleCon_s = lfsv_s.child(j);
          const LFSU_SINGLE_CON& lfsuSingleCon_n = lfsu_n.child(j);
          const LFSV_SINGLE_CON& lfsvSingleCon_n = lfsv_n.child(j);

          param.setIonSpecies(j);

          op.jacobian_skeleton(ig, lfsuSingleCon_s, x_s, lfsvSingleCon_s,
                                lfsuSingleCon_n, x_n, lfsvSingleCon_n,
                                mat_ss, mat_sn, mat_ns, mat_nn);
        }
      }
      */
    }

    // boundary integral depending on test and ansatz functions
    // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
    template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_boundary (const IG& ig,
                         const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                         R& r_s) const
    {
      // local function space for the NUMBER_OF_SPECIES components (assume Galerkin scheme U=V)
      typedef typename LFSU::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_SINGLE_CON;

      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        const LFSU_SINGLE_CON& lfsuSingleCon_s = lfsu_s.child(j);
        const LFSV_SINGLE_CON& lfsvSingleCon_s = lfsv_s.child(j);

        param.setIonSpecies(j);

        op.alpha_boundary(ig, lfsuSingleCon_s, x_s, lfsvSingleCon_s, r_s);
      }
    }

    template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_boundary (const IG& ig,
                            const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                            M& mat_ss) const
    {
      // local function space for the NUMBER_OF_SPECIES components (assume Galerkin scheme U=V)
      typedef typename LFSU::template Child<0>::Type LFSU_SINGLE_CON;
      typedef typename LFSV::template Child<0>::Type LFSV_SINGLE_CON;

      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        const LFSU_SINGLE_CON& lfsuSingleCon_s = lfsu_s.child(j);
        const LFSV_SINGLE_CON& lfsvSingleCon_s = lfsv_s.child(j);

        param.setIonSpecies(j);

        op.jacobian_boundary(ig, lfsuSingleCon_s, x_s, lfsvSingleCon_s, mat_ss);
      }
    }

    // volume integral depending only on test functions
    template<typename EG, typename LFSV, typename R>
    void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
    {
      // local function space for the NUMBER_OF_SPECIES components (assume Galerkin scheme U=V)
      typedef typename LFSV::template Child<0>::Type LFSV_SINGLE_CON;

      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        const LFSV_SINGLE_CON& lfsvSingleCon = lfsv.child(j);

        param.setIonSpecies(j);

        op.lambda_volume(eg, lfsvSingleCon, r);
      }
    }

    //! set time in parameter class
    void setTime (double t)
    {
      param.setTime(t);
    }

  private:
    T& param;  // parameter class
    OP& op;    // local operator
};


#endif /* NERNST_PLANCK_POWER_OPERATOR_HH */
