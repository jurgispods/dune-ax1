/*
 * step_initial.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_ES_INITIAL_HH
#define DUNE_AX1_ACME2CYL_ES_INITIAL_HH

#include <tuple>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>

template<typename GV, typename RF, int dim>
class ESInitialPot : public Dune::PDELab::AnalyticGridFunctionBase<
                                  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
                                  ESInitialPot<GV,RF,dim> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, ESInitialPot<GV,RF,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    ESInitialPot(const GV& gv_, const Acme2CylParameters& params_)
    : BaseT(gv_),
      gv(gv_),
      params(params_),
      time(0.0),
      useLoadedData(params.boundary.get("useTimeDependentBoundaryValues",false)),
      boundaryLocation(params.boundary.get("loadBoundary",std::string("bottom"))),
      isLoadBoundaryDirichlet(params.isBoundaryDirichlet_Potential(boundaryLocation)),
      filename(params.boundary.get("boundaryLoadFilename","boundary_pot_bottom.dat")),
      pot_in(filename),
      nLinesPerTimeStep(0),
      boundaryValues(),
      epsTime(params.boundary.get("boundaryTimePrecision",1e-4))
    {
      if(useLoadedData)
      {
        if(boundaryLocation != "top" && boundaryLocation != "bottom" && boundaryLocation != "membrane")
          DUNE_THROW(Dune::NotImplemented, "Loading left and right boundaries currently not implemented!");

        if(! isLoadBoundaryDirichlet)
          DUNE_THROW(Dune::NotImplemented, "Loading boundary value types other than Dirichlet are not implemented"
              << " for configuration 'ES'!");

        if(! pot_in.good())
          DUNE_THROW(Dune::Exception, "Simulation state file " << filename << " could not be found!");

        debug_info << "[ESInitialPot::ESInitialPot()] Using file '" << filename
                    << "' for loading boundary data!" << std::endl;

        double t = -1;

        std::string line;
        std::getline(pot_in, line);
        int nLine = 0;

        debug_jochen << "Searching for initial time value " << time << std::endl;

        // Search for time step
        while(t < time && std::abs(t-time) > epsTime && pot_in.good())
        {
         std::getline(pot_in, line);
         nLine++;

         size_t pos = line.find("# time:");
         if(pos == std::string::npos)
           continue;
         else
           debug_jochen << "Found '# time' token in line " << (nLine+1) << std::endl;

         std::stringstream line_str(line.substr(8));
         line_str >> t;

         debug_verb << "Checking time value t=" << t << "..." << std::endl;
        }
        // We reached the end of the file while reading; assume there is a newline at the end of the file
        if(! pot_in.good())
        {
          nLinesPerTimeStep = nLine-1;
        } else {
          // We found a line starting a new time value block; substract 3 (current line + 2 newlines)
          nLinesPerTimeStep = nLine-3;
        }
        debug_jochen << "Number of lines per time step: " << nLinesPerTimeStep << std::endl;
        boundaryValues.resize(nLinesPerTimeStep);

        // Go back to the beginning
        pot_in.clear();
        pot_in.seekg(0, std::ios::beg);
      }
    }

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      bool top = x[1] + 1e-6 > params.yMax();
      bool bottom = x[1] - 1e-6 < params.yMin();
      bool left = x[0] - 1e-6 < params.xMin();
      bool right = x[0] + 1e-6 > params.xMax();
      // TODO Extend this to the case of multiple membranes / membrane element layers
      bool membrane = x[1] + 1e-6 > params.yMemb()[0] && x[1] - 1e-6 < params.yMemb()[0] + params.dMemb();

      //debug_jochen << "ESInitialPot: About to evaluate at " << x << std::endl;
      // Always FUCKING assign a value to y!! This function is not only used to evaluate boundary values, but also
      // for interpolating initial values on the whole domain. It is crucial to set y to some value in any case,
      // as otherwise bad things might happen, as y contains undefined data.
      y = 0.0;

      if (useLoadedData && isLoadBoundaryDirichlet)
      {
        if ((bottom && boundaryLocation == "bottom")
            || (membrane && boundaryLocation == "membrane")
            || (top && boundaryLocation == "top"))
        {

          double x_load = -1;
          double y_load = -1;
          double value = 0;
          // Search for x coordinate
          int i = 0;
          while (x_load < x[0] && i<boundaryValues.size())
          {
            x_load = std::get<0>(boundaryValues[i]);
            y_load = std::get<1>(boundaryValues[i]);
            value = std::get<2>(boundaryValues[i]);
            i++;
          }
          debug_jochen << "  At x=" << x[0] << ", found old coordinate x="
              << x_load << " -> value: " << value << std::endl;

          // Make sure we have found a matching x coordinate
          if (std::abs(x_load - x[0]) > 1e-6)
            DUNE_THROW(Dune::Exception,
                "Could not find boundary value for x coordinate " << x[0] << ", found x=" << x_load << " instead (difference: " << std::abs(x_load-x[0]) << ")!");

          y = value;

          // convert [mV] -> [1] (not necessary when loading from file)
          //y *= con_e / (1.0e3 * con_k * 279.45);

//        } else {
//          debug_jochen << "No loaded data => no evaluation! y=" << y << std::endl;
//        }

          return;
        }
      } else {
        // Not loading data:

        // Deactivated
        if(false)
        {
          // Use a smooth Gaussian function
          const double baseline = -65;
          const double a = 110; // amplitude
          const double x0 = -2e6; // negative start coordinate of traveling signal
          const double v = 1e3; // velocity in nm/µs, here equal to 1m/s
          const double b = x0 + time * v; // center of the peak
          const double c = 1e6; // standard deviation (in nm)

          // evaluate Gaussian function
          double exponent = -(x[0] - b) * (x[0] - b) / (2 * c * c);
          y = baseline + a * std::exp(exponent);

          // convert [mV] -> [1]
          y = convertFrom_mV(y);
          return;
        }

        // Default: return parameter file value
        // TODO Put the whole functionality from config classes into one single file, see LaplaceParameters!
        // ES
        if(! params.useMembrane() || x[1] + 1e-6 > params.yMemb()[0] + params.dMemb())
        {
          y = convertFrom_mV(params.solution_ex.get("pot", 0.0));
          return;
        }
        // CY
        if(x[1] - 1e-6 < params.yMemb()[0] )
        {
          y = convertFrom_mV(params.solution_in.get("pot", -65.0));
          //debug_jochen << "g= " << y << std::endl;
          return;
        }

        DUNE_THROW(Dune::Exception, "Node is neither ES nor CYTOSOL!");
      }
    }

    // force time to be zero as this is the initial condition
    void setTime (double t)
    {
      time = t;

      // New method: Search for next time step in provided file
      if(useLoadedData)
      {
        assert(pot_in.good());
        //debug_jochen << "Trying to read time=" << time << "@ x=" << x[0] << ", skipping first "
        //    << startLine << " lines!" << std::endl;

        debug_info << "= Loading '" << boundaryLocation << "' boundary data from file for requested time "
            << time << "..." << std::endl;

        double t = -1;

        std::string line;

        // Search for time step
        while(t < time && std::abs(t-time) > epsTime && pot_in.good())
        {
         std::getline(pot_in, line);

         size_t pos = line.find("# time:");
         if(pos == std::string::npos)
           continue;

         std::stringstream line_str(line.substr(8));
         line_str >> t;

         debug_verb << "Checking time value t=" << t << "..." << std::endl;
        }
        debug_jochen << "At requested time " << time << ", found old timestep t=" << t << std::endl;

        // Make sure we have found a matching timestep
        if(std::abs(t-time) > epsTime)
         DUNE_THROW(Dune::Exception, "Could not find boundary value for time "
             << time  << ", found t=" << t << " instead!");

        // Cache all boundary entries belonging to the current timestep for quicker loading
        double x_load = -1;
        double y_load = -1;
        double value = 0;
        for(int i=0; i<nLinesPerTimeStep && pot_in.good(); i++)
        {
          std::getline(pot_in, line);

          std::stringstream line_str(line);
          line_str >> x_load >> y_load >> value;
          //debug_jochen << " Processing line " << line << std::endl;
          debug_verb << " Found old coordinate x=" << x_load << std::endl;

          // TODO Optimization: Only push those coordinates that are actually present on this processor!
          boundaryValues[i] = std::tuple<double,double,double>(x_load,y_load,value);
        }
      }
    }

    RF convertFrom_mV(RF pot) const
    {
      RF y = pot * (con_e / (1.0e3 * con_k * 279.45));
      return y;
    }

    // Dummy method to comply with CustomBoundaryPotential interface
    void setPotValues(const std::vector<RF>& potValues_)
    {}

  private:
    const GV& gv;
    const Acme2CylParameters& params;
    RF time;
    const bool useLoadedData;
    const std::string boundaryLocation;
    const bool isLoadBoundaryDirichlet;
    std::string filename;
    std::ifstream pot_in;
    int nLinesPerTimeStep;
    std::vector<std::tuple<double,double,double> > boundaryValues;
    double epsTime;
};

#endif /* DUNE_AX1_ACME2CYL_STEP_INITIAL_HH */
