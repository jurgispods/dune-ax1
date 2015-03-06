/*
 * no_analytical_solution.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_NO_ANALYTICAL_SOLUTION_HH
#define DUNE_AX1_NO_ANALYTICAL_SOLUTION_HH

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>

template<typename GV, typename RField, int dim>
class NoAnalyticalSolution :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RField,dim>,
  NoAnalyticalSolution<GV,RField,dim> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RField,dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, NoAnalyticalSolution<GV,RField,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    typedef RField RF;
    typedef RangeType RT;

    NoAnalyticalSolution(const GV & gv, const Acme2CylParameters& params_)
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
