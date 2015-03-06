/*
 * hamburger_initial.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_HAMBURGER_INITIAL_HH
#define DUNE_AX1_HAMBURGER_INITIAL_HH

#include <dune/ax1/acme1MD/common/acme1MD_parametertree.hh>

template<typename GV, typename RF, int dim, typename HAMBURGER_POT>
class HamburgerInitialPot : public HAMBURGER_POT
{
  public:
    typedef typename HAMBURGER_POT::Traits Traits;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    HamburgerInitialPot(const GV& gv_, const Acme1MDParameters& params_)
    : HAMBURGER_POT(gv_, params_)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      HAMBURGER_POT::evaluateGlobal(x, y);
    }

    // force time to be zero as this is the initial condition
    void setTime (double t)
    {
      this->time = 0.0;
    }
};

template<typename GV, typename RF, int dim, typename HAMBURGER_CON>
class HamburgerInitialCon : public HAMBURGER_CON
{
  public:
    typedef typename HAMBURGER_CON::Traits Traits;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    HamburgerInitialCon(const GV& gv_, const Acme1MDParameters& params_)
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

#endif /* DUNE_AX1_HAMBURGER_INITIAL_HH */
