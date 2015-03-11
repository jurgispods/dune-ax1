/*
 * step_initial.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_STEP_INITIAL_HH
#define DUNE_AX1_STEP_INITIAL_HH

#include <dune/ax1/acme1/common/acme1_parametertree.hh>

template<typename GV, typename RF, int dim, typename STEP_CON>
class StepInitial : public STEP_CON
{
  public:
    typedef typename STEP_CON::Traits Traits;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    StepInitial(const GV& gv_, const Acme1Parameters& params_)
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
