/*
 * hamburger_initial.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME0_HAMBURGER_INITIAL_HH
#define DUNE_AX1_ACME0_HAMBURGER_INITIAL_HH

#include <dune/ax1/acme0/common/acme0_parametertree.hh>

template<typename GV, typename RF, int dim, typename HAMBURGER_CON>
class HamburgerInitial : public HAMBURGER_CON
{
  public:
    typedef typename HAMBURGER_CON::Traits Traits;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    HamburgerInitial(const GV& gv_, const Acme0Parameters& params_)
    : HAMBURGER_CON(gv_, params_)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      HAMBURGER_CON::evaluateGlobal(x, y);
    }

    // force time to be zero as this is the initial condition
    void setTime (double t)
    {
      this->time = 0.0;
    }
};

#endif /* DUNE_AX1_ACME0_HAMBURGER_INITIAL_HH */
