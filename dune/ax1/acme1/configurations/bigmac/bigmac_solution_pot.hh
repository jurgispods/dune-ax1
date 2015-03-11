/*
 * bigmac_solution_pot.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_BIGMAC_SOLUTION_POT_HH
#define DUNE_AX1_BIGMAC_SOLUTION_POT_HH

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/acme1/common/acme1_parametertree.hh>

template<typename GV, typename RF, int dim>
class BigmacPot :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
  BigmacPot<GV,RF,dim> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, BigmacPot<GV,RF,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    BigmacPot(const GV& gv_, const Acme1Parameters& params_)
    : BaseT(gv_),
      gv(gv_),
      params(params_),
      time(0.0)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
    	//
    	double A = 1.0;
    	double v = 1.0;

      y[0] =   2.0 * log ( cosh( A * ( x   - v * time ) ) ) - x   * v
      			 - 2.0 * log ( cosh( A * ( 5.0 - v * time ) ) ) + 5.0 * v;
    }

    inline const GV& getGridView () const
    {
     return gv;
    }

    // set time for subsequent evaluation
    void setTime (double t)
    {
     time = t;
    }

  private:
    const GV& gv;
    const Acme1Parameters& params;
    RF time;
};


#endif /* DUNE_AX1_BIGMAC_SOLUTION_POT_HH */
