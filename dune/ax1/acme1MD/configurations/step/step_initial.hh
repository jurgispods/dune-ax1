/*
 * step_initial.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_STEP_INITIAL_HH
#define DUNE_AX1_STEP_INITIAL_HH

#include <dune/ax1/acme1MD/common/acme1MD_parametertree.hh>

template<typename GV, typename RF, int dim, typename STEP_POT>
class StepInitialPot : public STEP_POT
{
  public:
    typedef typename STEP_POT::Traits Traits;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    StepInitialPot(const GV& gv_, const Acme1MDParameters& params_)
    : STEP_POT(gv_, params_)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      STEP_POT::evaluateGlobal(x, y);
    }

    // force time to be zero as this is the initial condition
    void setTime (double t)
    {
      this->time = 0.0;
    }
};

template<typename GV, typename RF, int dim, typename STEP_CON>
class StepInitialCon : public STEP_CON
{
  public:
    typedef typename STEP_CON::Traits Traits;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    StepInitialCon(const GV& gv_, const Acme1MDParameters& params_)
    : STEP_CON(gv_, params_)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      STEP_CON::evaluateGlobal(x, y);
    }

    // force time to be zero as this is the initial condition
    void setTime (double t)
    {
      this->time = 0.0;
    }
};

#endif /* DUNE_AX1_STEP_INITIAL_HH */
