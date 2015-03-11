/*
 * hamburger_solution_con.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_HAMBURGER_SOLUTION_CON_HH
#define DUNE_AX1_HAMBURGER_SOLUTION_CON_HH

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/acme1/common/acme1_parametertree.hh>

template<typename GV, typename RF, int dim>
class HamburgerCon :
  public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
    HamburgerCon<GV,RF,dim> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, HamburgerCon<GV,RF,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    HamburgerCon(const GV& gv_, const Acme1Parameters& params_)
    : BaseT(gv_),
      gv(gv_),
      params(params_),
      time(0.0)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      //
      double lambda = 0.2*con_pi;
      double A = -0.5*con_pi;
      double B = 1.2;

      y[0] = 2.0*lambda*lambda*(1.0+B*sin(lambda*x+A)*exp(lambda*lambda*time))/std::pow((B*exp(lambda*lambda*time)+sin(lambda*x+A)),2);
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
    const Acme1Parameters& params;
  protected:
    RF time;
};

#endif /* DUNE_AX1_HAMBURGER_SOLUTION_CON_HH */
