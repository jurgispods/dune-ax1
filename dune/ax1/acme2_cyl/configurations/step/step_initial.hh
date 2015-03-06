/*
 * step_initial.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_STEP_INITIAL_HH
#define DUNE_AX1_ACME2CYL_STEP_INITIAL_HH

#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>

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

    StepInitialPot(const GV& gv_, const Acme2CylParameters& params_)
    : BaseT(gv_),
      gv(gv_),
      params(params_),
      time(0.0)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      bool top = x[1]+1e-6 > params.yMax();
      bool bottom = x[1]-1e-6 < params.yMin();
      bool left = x[0]-1e-6 < params.xMin();
      bool right = x[0]+1e-6 > params.xMax();

      if(top || bottom || left || right)
      {
        int cdim = 0;
        if(left || right) cdim = 1;

        std::stringstream filename;
        filename << params.getOutputPrefix();
        if(USE_CYLINDER_COORDINATES)
        {
          filename << "acme2_cyl_boundary_";
        } else {
          filename << "acme2_2d_boundary_";
        }

        if(bottom) filename << "bottom";
        else if(left) filename << "left";
        else if(right) filename << "right";
        else if(top) filename << "top";
        filename << ".dat";

        std::ifstream in(filename.str());
        if(! in.good())
            DUNE_THROW(Dune::Exception, "Simulation state file " << filename.str() << " could not be found!");

        double xx = -1, yy = -1;

        while(xx < x[cdim]-1e-6)
        {
          in >> xx >> yy;
        }
        // Hopefully we found a matching x position
        //debug_jochen << "At boundary coordinate " << x << ", using value " << yy << std::endl;
        y = yy;

        return;
      }

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
    const Acme2CylParameters& params;
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

    StepInitialCon(const GV& gv_, const Acme2CylParameters& params_)
    : BaseT(gv_),
      gv(gv_),
      params(params_),
      time(0.0),
      initConcEx(NUMBER_OF_SPECIES),
      initConcIn(NUMBER_OF_SPECIES)
    {
      for(int i=0; i<NUMBER_OF_SPECIES; ++i)
      {
        initConcEx[i] = params.solution_ex.get(ION_NAMES[i], -1.);
        initConcIn[i] = params.solution_in.get(ION_NAMES[i], -1.);
      }
    }

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      double xmin = params.xMin();
      double xmax = params.xMax();

      double ymin = params.yMin();
      double ymax = params.yMax();

      DomainType pos;
      // Abuse config file parameters stimulation.position_{x,y} for position
      pos[0] = params.stimulation.get("position_x", xmin + 0.5*(xmax-xmin));
      pos[1] = params.stimulation.get("position_y", ymin + 0.5*(ymax-ymin));

      // Abuse config file parameter 'dMemb' to be used as the step width/height
      double d = params.general.get("d_memb", 10.0);
      double EPS = 1e-6;

      bool xInside = x[0] > pos[0]-d && x[0] < pos[0]+d;
      bool yInside = x[1] > pos[1]-d && x[1] < pos[1]+d;
      bool inside = xInside && yInside;

      bool xBoundary = (std::abs(x[0]-pos[0]+d) < EPS || std::abs(x[0]-pos[0]-d) < EPS);
      bool yBoundary = (std::abs(x[1]-pos[1]+d) < EPS || std::abs(x[1]-pos[1]-d) < EPS);
      bool boundary = (xBoundary && yBoundary) || (xBoundary && yInside) || (yBoundary && xInside);

      bool lowerBoundary = yBoundary && xInside && (x[1]-EPS < params.yMin());

      y = 0.0;

      // Using the following 'sophisticated' algorithm (working for a rectangle pulse on a 2D Cartesian
      // rectangle grid using linear conforming finite elements only),
      // the total source should equal rectangle area times source (density) value!
      // This should maximize the comparability with analytical solutions in Matlab
      if(inside || lowerBoundary)
      {
        // Use extracellular sodium concentration as the concentration step value
        y = 1.0 * initConcEx[Na];

        debug_jochen << "[StepInitialCon] @ " << x << ", y = " << y << std::endl;
      } else if(boundary)
      {
        // This is driving my nuts...
        double scale_factor = params.dYMin() / (std::max(params.dX(),params.dY()) + params.dYMin());

        if(xBoundary && yBoundary)
        {
          // Add corners of source area, take 1/4 of source
          y = scale_factor * scale_factor * initConcEx[Na];
        } else {
          // Add edges of source area, take 1/2 of source
          y = scale_factor * initConcEx[Na];
        }

        debug_jochen << "[StepInitialCon] @ " << x << ", y = " << y << std::endl;
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
    const Acme2CylParameters& params;
    RF time;

    std::vector<RF> initConcEx;
    std::vector<RF> initConcIn;
};

#endif /* DUNE_AX1_ACME2CYL_STEP_INITIAL_HH */
