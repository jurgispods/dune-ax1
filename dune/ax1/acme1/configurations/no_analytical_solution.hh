/*
 * no_analytical_solution.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_NO_ANALYTICAL_SOLUTION_HH
#define DUNE_AX1_NO_ANALYTICAL_SOLUTION_HH

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/acme1/common/acme1_parametertree.hh>

template<typename GV, typename RF, int dim>
class NoAnalyticalSolution :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
  NoAnalyticalSolution<GV,RF,dim> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, NoAnalyticalSolution<GV,RF,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    NoAnalyticalSolution(const GV & gv, const Acme1Parameters& params_)
    : BaseT(gv)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      y[0] = 0.;
    }

    /*
    inline const GV& getGridView () const
    {
    }
    */

    // set time for subsequent evaluation
    void setTime (double t)
    {}
};

#endif /* DUNE_AX1_NO_ANALYTICAL_SOLUTION_HH */
