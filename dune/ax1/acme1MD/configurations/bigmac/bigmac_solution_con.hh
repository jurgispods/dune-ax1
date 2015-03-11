/*
 * bigmac_solution_con.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_BIGMAC_SOLUTION_CON_HH
#define DUNE_AX1_BIGMAC_SOLUTION_CON_HH

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/acme1MD/common/acme1MD_parametertree.hh>

template<typename GV, typename RF, int dim>
class BigmacCon :
  public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
    BigmacCon<GV,RF,dim> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, BigmacCon<GV,RF,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    BigmacCon(const GV& gv_, const Acme1MDParameters& params_)
    : BaseT(gv_),
      gv(gv_),
      params(params_),
      time(0.0)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      //
    	double A = 1.0;
    	double B = 0.0;
    	double v = 1.0;

      y[0] = - 2.0 * A * A * ( 1.0 - std::pow(tanh( A * ( x - v * time ) + B ),2) );
    }

    inline const GV& getGridView () const
    {
      return gv;
    }

    // set time for subsequent evaluation
    virtual void setTime (double t)
    {
      time = t;
    }


  private:
    const GV& gv;
    const Acme1MDParameters& params;
  protected:
    RF time;
};

#endif /* DUNE_AX1_BIGMAC_SOLUTION_CON_HH */
