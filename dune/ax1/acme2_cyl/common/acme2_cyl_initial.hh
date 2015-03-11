#ifndef DUNE_AX1_ACME2CYL_INITIAL_HH
#define DUNE_AX1_ACME2CYL_INITIAL_HH

#include <valarray>

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>
#include <dune/ax1/acme2_cyl/configurations/zero_potential.hh>

//! This is a generic wrapper class that internally calls the INITIAL_GF object to
//! calculate initial values
template<typename GV, typename RF, int numSpecies, typename PHYSICS, typename INITIAL_GF>
class InitialCon
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<GV,RF,numSpecies,Dune::FieldVector<RF,numSpecies> >,
    InitialCon<GV,RF,numSpecies,PHYSICS,INITIAL_GF> >
{

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,numSpecies,Dune::FieldVector<RF,numSpecies> > Traits;

  //! construct from grid view
  InitialCon (const GV& gv_, PHYSICS& physics_, INITIAL_GF& initialGF_)
  : gv(gv_),
    physics(physics_),
    initialGF(initialGF_),
    params(physics_.getParams()),
    time(0.0),
    initConcEx(numSpecies),
    initConcIn(numSpecies)
  {
    for(int i=0; i<numSpecies; ++i)
    {
      initConcEx[i] = params.solution_ex.get(ION_NAMES[i], -1.);
      initConcIn[i] = params.solution_in.get(ION_NAMES[i], -1.);
    }
    std::string initConcStr = params.initConc();
  }

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                              typename Traits::RangeType& y) const
  {
    int elemIndex = physics.getElementIndex(e);
    int subdomainIndex = physics.getSubdomainIndex(e);

    for(int i=0; i<numSpecies; ++i)
    {
      switch(subdomainIndex)
      {
        case CYTOSOL:
        {
          y[i] = initConcIn[i];
          break;
        }
        case ES:
        {
          y[i] = initConcEx[i];
          break;
        }
        case MEMBRANE:
        {
          // No charge carriers inside the membrane!
          y[i] = 0.0;
          return;
        }
        default:
          DUNE_THROW(Dune::Exception, "Element has an unknown subdomain index!");
      }
    }

    // Override default behaviour by evaluating custom initial gridfunction
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    const int dim = Traits::GridViewType::Grid::dimension;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    initialGF.evaluateGlobal(x,y);
  }
  
  inline const GV& getGridView () const
  {
    return gv;
  }

  // set time for subsequent evaluation
  void setTime (double t)
  {
    time = t;
    initialGF.setTime(t);
  }

private:
  PHYSICS physics;
  INITIAL_GF& initialGF;
  const GV& gv;
  const Acme2CylParameters& params;
  RF time;
  std::vector<RF> initConcEx;
  std::vector<RF> initConcIn;
};


// =============== Different functions for initial concentrations ===============================

/**
 * Trivial class that does nothing; can be used to keep default 'bulk' initial values
 * (constant values for each electrolyte subdomain)
 */
template<typename GV, typename RF, int numSpecies>
class InitialVoid :
  public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,numSpecies>,
  InitialVoid<GV,RF,numSpecies> >
{
  public:
      typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,numSpecies> Traits;
      typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InitialVoid<GV,RF,numSpecies> > BaseT;

      typedef typename Traits::DomainType DomainType;
      typedef typename Traits::RangeType RangeType;

      InitialVoid(const GV& gv, const Acme2CylParameters& params_)
      : BaseT(gv)
      {}

      // Does nothing
      inline void evaluateGlobal(const DomainType & x, RangeType & y) const
      {
      }

      // Does nothing
      void setTime (double t)
      {
      }

  private:
};


/**
 * Step initial concentrations: constant concentrations for each domain (intra-/extracellular)
 */
template<typename GV, typename RF, int numSpecies>
class InitialStep :
  public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,numSpecies>,
  InitialStep<GV,RF,numSpecies> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,numSpecies> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InitialStep<GV,RF,numSpecies> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    InitialStep(const GV& gv, const Acme2CylParameters& params_)
    : BaseT(gv),
      params(params_),
      yMemb(params.yMemb()),
      initConcEx(numSpecies),
      initConcIn(numSpecies),
      smoothStep(false)
    {
      for(int i=0; i<numSpecies; ++i)
      {
        initConcEx[i] = params.solution_ex.get(ION_NAMES[i], -1.);
        initConcIn[i] = params.solution_in.get(ION_NAMES[i], -1.);
      }
      std::string initConcStr = params.initConc();
      if(initConcStr == "smoothStep")
      {
        smoothStep = true;
      }
    }

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      if(params.nMembranes() > 1)
        DUNE_THROW(Dune::NotImplemented, "Step configuration works for one membrane only!");


      for ( int i=0; i<numSpecies; ++i )
      {
        // Intracellular domain
        if (x[1] < yMemb[0] || std::abs(x[1] - yMemb[0]) < 1e-12)
        {
          y[i] = initConcIn[i];
        }
        // Extracellular domain
        else if(x[1] > yMemb[0] || std::abs(x[1] - yMemb[0]) < 1e-12)
        {
          y[i] = initConcEx[i];
        // Membrane
        } else {
          // Smooth discontinuity in initial concentrations
          if(smoothStep)
          {
            typename Traits::DomainFieldType y_scale = (x[1] - yMemb[0]) / (params.dMemb());
            y[i] = initConcIn[i] + (initConcEx[i]-initConcIn[i]) * (y_scale*y_scale*(3-2*y_scale));
          // No charge carriers inside membrane
          } else {
            if (x[1] < yMemb[0] + 0.5 * params.dMemb()) y[i] = initConcIn[i];
            else                                              y[i] = initConcEx[i];
          }
        }
      }
    }

    // Does nothing
    void setTime (double t)
    {
    }


  private:
    const Acme2CylParameters& params;
    const std::vector<double> yMemb;
    std::vector<RF> initConcEx;
    std::vector<RF> initConcIn;
    bool smoothStep;
};

/**
 * SquarePulse initial concentrations: square pulse with constant concentrations for specified x-range, 0 elsewhere
 */
template<typename GV, typename RF, int numSpecies>
class InitialSquarePulse :
  public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,numSpecies>,
  InitialSquarePulse<GV,RF,numSpecies> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,numSpecies> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InitialSquarePulse<GV,RF,numSpecies> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    InitialSquarePulse(const GV& gv, const Acme2CylParameters& params_)
    : BaseT(gv),
      params(params_)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      if(params.nMembranes() > 1)
        DUNE_THROW(Dune::NotImplemented, "Step configuration works for one membrane only!");

      double dpulse = 5.0;
      double xCenter = 0.5 * (params.xMax() - params.xMin());
      double yCenter = 0.5 * (params.yMemb()[0] - params.yMin());
      double pulse_amplitude = 1.0;

      y = 0.0;

      // square pulse
      if ( x[0] > xCenter-dpulse and x[0] < xCenter+dpulse
          && x[1] > yCenter-dpulse and x[1] < yCenter+dpulse)
      {
        y[0] = pulse_amplitude;
      } else {
        y[0] = 1.0;
      }

      if(numSpecies > 2)
      {
        // potassium
        y[1] = 1.0;
      }
      if(numSpecies > 2)
      {
        // chloride
        y[2] = 2.0;
      }
    }

    // Does nothing
    void setTime (double t)
    {
    }


  private:
    const Acme2CylParameters& params;
};

/**
 * TrianglePulse initial concentrations: triangle pulse for specified x-range, 0 elsewhere
 */
template<typename GV, typename RF, int numSpecies>
class InitialTrianglePulse :
  public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,numSpecies>,
  InitialTrianglePulse<GV,RF,numSpecies> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,numSpecies> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InitialTrianglePulse<GV,RF,numSpecies> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    InitialTrianglePulse(const GV & gv, const Acme2CylParameters& params_)
    : BaseT(gv),
      params(params_)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      double xpulse = 0.5 * (params.xMax() - params.dMemb());
      double pulse_amplitude = 1.0;

      // triangle pulse
      if ( x[0] == 0.0 )
      {
       y[0] = pulse_amplitude;
       if(numSpecies > 1) y[1] = y[0];
       if(numSpecies > 2) y[2] = y[0] * 2.0;
      }
      else if ( x[0] < 0.0 and x[0] > -xpulse )
      {
       y[0] = pulse_amplitude * ( 1.0 + x[0] / xpulse );
       if(numSpecies > 1) y[1] = y[0];
       if(numSpecies > 2) y[2] = y[0] * 2.0;
      }
      else if ( x[0] > 0.0 and x[0] < xpulse )
      {
       y[0] = pulse_amplitude * ( 1.0 - x[0] / xpulse );
       if(numSpecies > 1) y[1] = y[0];
       if(numSpecies > 2) y[2] = y[0] * 2.0;
      }
      else
      {
       y[0] = 0.0;
       if(numSpecies > 1) y[1] = 0.0;
       if(numSpecies > 2) y[2] = 0.0;
      }
    }

    // Does nothing
    void setTime (double t)
    {
    }

  private:
    const Acme2CylParameters& params;
};

// initial potential ###############################################################################

template<typename GV, typename RF, typename Physics, typename INITIAL_GF>
class InitialPot
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
    InitialPot<GV,RF,Physics,INITIAL_GF> >
{

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  InitialPot (const GV& gv_, Physics& physics_, INITIAL_GF& initialGF_, int intorder_=2) :
    gv(gv_),
    physics(physics_),
    initialGF(initialGF_)
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    //const int dim = Traits::GridViewType::Grid::dimension;
    //typedef typename Traits::GridViewType::Grid::ctype ctype;
    //Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);
    
    initialGF.evaluate(e,xlocal,y);
  }

  
  inline const GV& getGridView() const
  {
    return gv;
  }

  //! set time for subsequent evaluation
  void setTime (double t)
  {
    time = t;
    initialGF.setTime(t);
  }


private:
  const GV& gv;
  Physics & physics;
  INITIAL_GF& initialGF;
  RF time;

};
#endif /* DUNE_AX1_ACME2CYL_INITIAL_HH */
