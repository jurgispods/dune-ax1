/*
 * cheeseburger_solution_con.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_CHEESEBURGER_SOLUTION_CON_HH
#define DUNE_AX1_CHEESEBURGER_SOLUTION_CON_HH

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/acme1/common/acme1_parametertree.hh>

template<typename GV, typename RF, int dim>
class CheeseburgerCon :
  public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
    CheeseburgerCon<GV,RF,dim> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, CheeseburgerCon<GV,RF,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    CheeseburgerCon(const GV& gv_, const Acme1Parameters& params_)
    : BaseT(gv_),
      gv(gv_),
      params(params_),
      time(0.0)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      //

      y[0] = 1.0 / ( 1.0 + time );
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

#endif /* DUNE_AX1_CHEESEBURGER_SOLUTION_CON_HH */
