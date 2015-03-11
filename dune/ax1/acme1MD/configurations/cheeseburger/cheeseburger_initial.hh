/*
 * cheeseburger_initial.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_CHEESEBURGER_INITIAL_HH
#define DUNE_AX1_CHEESEBURGER_INITIAL_HH

#include <dune/ax1/acme1MD/common/acme1MD_parametertree.hh>

template<typename GV, typename RF, int dim, typename CHEESEBURGER_POT>
class CheeseburgerInitialPot : public CHEESEBURGER_POT
{
  public:
    typedef typename CHEESEBURGER_POT::Traits Traits;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    CheeseburgerInitialPot(const GV& gv_, const Acme1MDParameters& params_)
    : CHEESEBURGER_POT(gv_, params_)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      CHEESEBURGER_POT::evaluateGlobal(x, y);
    }

    // force time to be zero as this is the initial condition
    void setTime (double t)
    {
      this->time = 0.0;
    }
};

template<typename GV, typename RF, int dim, typename CHEESEBURGER_CON>
class CheeseburgerInitialCon : public CHEESEBURGER_CON
{
  public:
    typedef typename CHEESEBURGER_CON::Traits Traits;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    CheeseburgerInitialCon(const GV& gv_, const Acme1MDParameters& params_)
    : CHEESEBURGER_CON(gv_, params_)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      CHEESEBURGER_CON::evaluateGlobal(x, y);
    }

    // force time to be zero as this is the initial condition
    void setTime (double t)
    {
      this->time = 0.0;
    }
};

#endif /* DUNE_AX1_CHEESEBURGER_INITIAL_HH */
