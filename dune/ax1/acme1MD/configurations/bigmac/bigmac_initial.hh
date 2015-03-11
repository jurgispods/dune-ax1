/*
 * bigmac_initial.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_BIGMAC_INITIAL_HH
#define DUNE_AX1_BIGMAC_INITIAL_HH

#include <dune/ax1/acme1MD/common/acme1MD_parametertree.hh>

template<typename GV, typename RF, int dim, typename BIGMAC_POT>
class BigmacInitialPot : public BIGMAC_POT
{
  public:
    typedef typename BIGMAC_POT::Traits Traits;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    BigmacInitialPot(const GV& gv_, const Acme1MDParameters& params_)
    : BIGMAC_POT(gv_, params_)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      BIGMAC_POT::evaluateGlobal(x, y);
    }

    // force time to be zero as this is the initial condition
    void setTime (double t)
    {
      this->time = 0.0;
    }
};

template<typename GV, typename RF, int dim, typename BIGMAC_CON>
class BigmacInitialCon : public BIGMAC_CON
{
  public:
    typedef typename BIGMAC_CON::Traits Traits;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    BigmacInitialCon(const GV& gv_, const Acme1MDParameters& params_)
    : BIGMAC_CON(gv_, params_)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      BIGMAC_CON::evaluateGlobal(x, y);
    }

    // force time to be zero as this is the initial condition
    void setTime (double t)
    {
      this->time = 0.0;
    }
};

#endif /* DUNE_AX1_BIGMAC_INITIAL_HH */
