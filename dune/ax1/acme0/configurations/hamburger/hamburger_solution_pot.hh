/*
 * hamburger_solution_pot.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME0_HAMBURGER_SOLUTION_POT_HH
#define DUNE_AX1_ACME0_HAMBURGER_SOLUTION_POT_HH

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/acme0/common/acme0_parametertree.hh>

template<typename GV, typename RF, int dim>
class HamburgerPot :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
  HamburgerPot<GV,RF,dim> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, HamburgerPot<GV,RF,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    HamburgerPot(const GV& gv_, const Acme0Parameters& params_)
    : BaseT(gv_),
      gv(gv_),
      params(params_),
      time(0.0)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      double lambda = 0.2*con_pi;
      double A = -0.5*con_pi;
      double B = 1.2;
      double x0 = 5.0;

      y[0] = 2.0*(log( (B*exp(lambda*lambda*time)+sin(lambda*x +A)) /
                        (B*exp(lambda*lambda*time)+sin(lambda*x0+A)) ));
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
    const Acme0Parameters& params;
    RF time;
};


#endif /* DUNE_AX1_ACME0_HAMBURGER_SOLUTION_POT_HH */
