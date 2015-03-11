/*
 * step_initial.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2_STEP_INITIAL_HH
#define DUNE_AX1_ACME2_STEP_INITIAL_HH

#include <dune/ax1/acme2/common/acme2_parametertree.hh>

template<typename GV, typename RF, int dim>
class StepInitialPot : public Dune::PDELab::AnalyticGridFunctionBase<
                                  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
                                  StepInitialPot<GV,RF,dim> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, StepInitialPot<GV,RF,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    StepInitialPot(const GV& gv_, const Acme2Parameters& params_)
    : BaseT(gv_),
      gv(gv_),
      params(params_),
      time(0.0)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      y = 0.0;
    }

    // force time to be zero as this is the initial condition
    void setTime (double t)
    {
      time = 0.0;
    }

    // Dummy method to comply with CustomBoundaryPotential interface
    void setPotValues(const std::vector<RF>& potValues_)
    {
    }

  private:
    const GV& gv;
    const Acme2Parameters& params;
    RF time;
};

template<typename GV, typename RF, int dim>
class StepInitialCon :  public Dune::PDELab::AnalyticGridFunctionBase<
                                  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
                                  StepInitialCon<GV,RF,dim> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, StepInitialCon<GV,RF,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    StepInitialCon(const GV& gv_, const Acme2Parameters& params_)
    : BaseT(gv_),
      gv(gv_),
      params(params_),
      time(0.0)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      double xmin = params.xMin();
      double xmax = params.xMax();

      double ymin = params.yMin();
      double ymax = params.yMax();

      DomainType center;
      center[0] = xmin + 0.5*(xmax-xmin);
      center[1] = ymin + 0.5*(ymax-ymin);

      double d = 10.;

      if(x[0] > center[0]-d && x[0] < center[0]+d
          && x[1] > center[1]-d && x[1] < center[1]+d)
      {
        y = 1.0;
      } else {
        y = 0.0;
      }
    }

    // force time to be zero as this is the initial condition
    void setTime (double t)
    {
      time = 0.0;
    }

  private:
    const GV& gv;
    const Acme2Parameters& params;
    RF time;
};

#endif /* DUNE_AX1_ACME2_STEP_INITIAL_HH */
