/*
 * zero_potential.hh
 *
 *  Created on: Mar 30, 2012
 *      Author: jpods
 */
#ifndef DUNE_AX1_ZERO_POTENTIAL_HH
#define DUNE_AX1_ZERO_POTENTIAL_HH

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/acme1MD/common/acme1MD_parametertree.hh>

template<typename GV, typename RField, int dim>
class ZeroPotential :
  public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RField,dim>,
  ZeroPotential<GV,RField,dim> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RField,dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, ZeroPotential<GV,RField,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    typedef RField RF;
    typedef RangeType RT;

    ZeroPotential(const GV& gv, const Acme1MDParameters& params_)
    : BaseT(gv)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      y[0] = 0.;
    }

    /*
    inline const GV& getGridView () const
    {
    }
    */

    // set time for subsequent evaluation
    void setTime (double t)
    {
    }
};

#endif /* DUNE_AX1_ZERO_POTENTIAL_HH */
