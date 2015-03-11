#ifndef DUNE_AX1_ACME0_INITIAL_HH
#define DUNE_AX1_ACME0_INITIAL_HH

#include <valarray>

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/acme0/common/acme0_parametertree.hh>
#include <dune/ax1/acme0/configurations/hamburger/hamburger_config.hh>

template<typename GV, typename RF, int numSpecies>
class InitialCon
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<
    GV,RF,numSpecies,Dune::FieldVector<RF,numSpecies> >, InitialCon<GV,RF,numSpecies> >
{

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,numSpecies,Dune::FieldVector<RF,numSpecies> > Traits;

  //! construct from grid view
  InitialCon (const GV& gv_, const Acme0Parameters& params_)
  : gv(gv_),
    params(params_),
    step(gv,params),
    squarePulse(gv,params),
    trianglePulse(gv,params)
  {
    std::string initConcStr = params.initConc();
    if(initConcStr == "smoothStep")
    {
      step.setSmoothStep(true);
    }
  }

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                              typename Traits::RangeType& y) const
  {

    typedef typename Traits::GridViewType::Grid::ctype ctype;
    const int dim = Traits::GridViewType::Grid::dimension;

    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    std::string initConcStr = params.initConc();

    //std::cout << "##### " << elementMapper.map(e) << std::endl;
    
    // Switch for initial condition
    if(initConcStr == "step" || initConcStr == "smoothStep")
    {
      step.evaluateGlobal(x,y);
    } else if (initConcStr == "squarePulse") {
      squarePulse.evaluateGlobal(x,y);
    } else if(initConcStr == "trianglePulse") {
      trianglePulse.evaluateGlobal(x,y);
    } else {
      DUNE_THROW(Dune::Exception,
          "No initial concentration configuration could be found for input 'initConc=" << initConcStr << "'!");
    }
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


  // =============== Different functions for initial concentrations ===============================

  /**
   * Step initial concentrations: constant concentrations for each domain (intra-/extracellular)
   */
  template<int dim>
  class Step :
    public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
    Step<dim> >
  {
    public:
      typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
      typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, Step<dim> > BaseT;

      typedef typename Traits::DomainType DomainType;
      typedef typename Traits::RangeType RangeType;

      Step(const GV & gv, const Acme0Parameters& params_)
      : BaseT(gv),
        params(params_),
        initConcEx(NUMBER_OF_SPECIES),
        initConcIn(NUMBER_OF_SPECIES),
        smoothStep(false)
      {
        for(int i=0; i<NUMBER_OF_SPECIES; ++i)
        {
          initConcEx[i] = params.solution_ex.get(ION_NAMES[i], -1.);
          initConcIn[i] = params.solution_in.get(ION_NAMES[i], -1.);
        }
      }

      inline void evaluateGlobal(const DomainType & x, RangeType & y) const
      {
        bool useLogScaling = params.useLogScaling();
        // Left and right boundary of membrane
        DomainType x_l(-0.5 * params.general.get("d_memb", 50.));
        DomainType x_r( 0.5 * params.general.get("d_memb", 50.));

        for ( int i=0; i<numSpecies; ++i )
        {
          // Intracellular domain
          if (x[0] < x_l)
          {
            y[i] = initConcIn[i];
          }
          // Extracellular domain
          else if(x[0] > x_r)
          {
            y[i] = initConcEx[i];
          // Membrane
          } else {
            // Smooth discontinuity in initial concentrations
            if(smoothStep)
            {
              DomainType x_scale = (x-x_l)/(x_r-x_l);
              y[i] = initConcIn[i] + (initConcEx[i]-initConcIn[i]) * (x_scale*x_scale*(3-2*x_scale));

            } else {
              // No charge carriers inside membrane
              if(params.useMembrane())
              {
                if (x[0] < 0.0) y[i] = 0.0;
                else            y[i] = 0.0;
              } else {
                if (x[0] < 0.0) y[i] = initConcIn[i];
                else            y[i] = initConcEx[i];
              }
            }
          }
          if (useLogScaling) y[i] = std::log(y[i]);
        }
      }

      void setSmoothStep(bool smoothStep_)
      {
        smoothStep = smoothStep_;
      }


    private:
      const Acme0Parameters& params;
      std::vector<RF> initConcEx;
      std::vector<RF> initConcIn;
      bool smoothStep;
  };

  /**
   * SquarePulse initial concentrations: square pulse with constant concentrations for specified x-range, 0 elsewhere
   */
  template<int dim>
  class SquarePulse :
    public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
    SquarePulse<dim> >
  {
    public:
      typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
      typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, SquarePulse<dim> > BaseT;

      typedef typename Traits::DomainType DomainType;
      typedef typename Traits::RangeType RangeType;

      SquarePulse(const GV & gv, const Acme0Parameters& params_)
      : BaseT(gv),
        params(params_)
      {}

      inline void evaluateGlobal(const DomainType & x, RangeType & y) const
      {
        double xpulse = 0.5 * (params.xMax() - params.dMemb());
        double pulse_amplitude = 1.0;

        // square pulse

        if ( x[0] > -xpulse and x[0] < xpulse )
        {
          y[0] = pulse_amplitude;
          if(NUMBER_OF_SPECIES > 1) y[1] = 0.0;
          if(NUMBER_OF_SPECIES > 2) y[2] = 0.0;
        }
        else
        {
          y[0] = 0.0;
          if(NUMBER_OF_SPECIES > 1) y[1] = 0.0;
          if(NUMBER_OF_SPECIES > 2) y[2] = 0.0;
        }
      }


    private:
      const Acme0Parameters& params;
  };

  /**
   * TrianglePulse initial concentrations: triangle pulse for specified x-range, 0 elsewhere
   */
  template<int dim>
  class TrianglePulse :
    public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
    TrianglePulse<dim> >
  {
    public:
      typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
      typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, TrianglePulse<dim> > BaseT;

      typedef typename Traits::DomainType DomainType;
      typedef typename Traits::RangeType RangeType;

      TrianglePulse(const GV & gv, const Acme0Parameters& params_)
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
         if(NUMBER_OF_SPECIES > 1) y[1] = y[0];
         if(NUMBER_OF_SPECIES > 2) y[2] = y[0] * 2.0;
        }
        else if ( x[0] < 0.0 and x[0] > -xpulse )
        {
         y[0] = pulse_amplitude * ( 1.0 + x[0] / xpulse );
         if(NUMBER_OF_SPECIES > 1) y[1] = y[0];
         if(NUMBER_OF_SPECIES > 2) y[2] = y[0] * 2.0;
        }
        else if ( x[0] > 0.0 and x[0] < xpulse )
        {
         y[0] = pulse_amplitude * ( 1.0 - x[0] / xpulse );
         if(NUMBER_OF_SPECIES > 1) y[1] = y[0];
         if(NUMBER_OF_SPECIES > 2) y[2] = y[0] * 2.0;
        }
        else
        {
         y[0] = 0.0;
         if(NUMBER_OF_SPECIES > 1) y[1] = 0.0;
         if(NUMBER_OF_SPECIES > 2) y[2] = 0.0;
        }
      }

    private:
      const Acme0Parameters& params;
  };

private:
  const GV& gv;
  const Acme0Parameters& params;
  RF time;
  Step<NUMBER_OF_SPECIES> step;
  SquarePulse<NUMBER_OF_SPECIES> squarePulse;
  TrianglePulse<NUMBER_OF_SPECIES> trianglePulse;
};

// initial potential ###############################################################################

template<typename GV, typename RF, typename Physics, typename GF_CD>
class InitialPot
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >
  , InitialPot<GV,RF,Physics,GF_CD> >
{

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

  //! construct from grid view
  InitialPot (const GV& gv_, Physics & physics_, GF_CD& gfChargeDensity_) :
    gv(gv_), physics(physics_), gfChargeDensity(gfChargeDensity_), elementMapper(gv),
    xmax(physics.getParams().xMax()), xmin(physics.getParams().xMin()),
    calculateAll(false)
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);
    
    int elemIndex = elementMapper.map(e);

    /*
    // We are on the boundary
    if (x[0] < xmin+1.0e-6 || x[0] >  xmax-1.0e-6)
    {
      //std::valarray<RF> source;
      //std::valarray<RF> pos;
      //Tools::getSolutionVectorDG(gfChargeDensity, 0, pos, source);
      //source = source*physics.getPoissonConstant();

      y = calculateExactPotential(x);

      // Left boundary: Cytosol
      if (x[0] < xmin+1.0e-6 )
      {
        //y=Tools::multipoleExpansion1D(physics.position, source, physics.permittivity, -4.0, 3);
        //y=Tools::potentialExact(physics.getPosition(), source, physics.getPermittivity(), xmin);
        //debug_verb << "### [left] boundary potential values: " << y << std::endl;
        //y = 0.0;
      }

      // Right boundary: Extracellular fluid
      if ( x[0] >  xmax-1.0e-6 )
      {
        //y=Tools::multipoleExpansion1D(physics.position, source, physics.permittivity, 4.0, 3);
        //y=Tools::potentialExact(physics.getPosition(), source, physics.getPermittivity(), xmax);
        //debug_verb << "### [right] boundary potential values: " << y << std::endl;
        y = 0.0;
      }

      // for testing: set potential to zero
      //y = 0.0;
    }
    // We are not on the boundary
    else
    {
      if(not calculateAll)
      {
        y = 0.0;
      } else {
        // Calculate potential values also in interior!
        y = calculateExactPotential(x);
      }

    }
    //y = ( 500.0 + x[0] ) / 100.0;

     */
    y = 0.0;

  }


  inline RF calculateExactPotential(const typename Traits::DomainType& x) const
  {
    const int dim = GV::dimension;
    const int intorder = 2;

    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef typename Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;
    typedef typename GF_CD::Traits::DomainType DT;
    typedef typename GF_CD::Traits::DomainFieldType DF;
    typedef typename GF_CD::Traits::RangeType RT;

    ElementMapper elementMapper(gv);
    RF potExact = 0.0;
    for (ElementIterator eit=gv.template begin<0>(); eit!=gv.template end<0>(); ++eit)
    {
      int elemIndex = elementMapper.map(*eit);

      RF permittivity = physics.getPermittivity(physics.getSubdomainIndex(elemIndex));

      // select quadrature rule
      Dune::GeometryType gt = eit->geometry().type();
      const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

      for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        RF weight = it->weight();
        RF integrationElement = eit->geometry().integrationElement(it->position());

        RT chargeDensity(0.0);
        gfChargeDensity.evaluate(*eit, it->position(), chargeDensity);
        //debug_verb << "cd = " << chargeDensity << std::endl;

        RF source = 0.0;
        source = chargeDensity[0]*physics.getPoissonConstant();

        source *= (-0.5/permittivity) * std::abs(x - eit->geometry().global(it->position()));

        potExact += (source * weight * integrationElement);
      }
    }
    return potExact;
  }
  
  inline const GV& getGridView() const
  {
    return gv;
  }

  //! set time for subsequent evaluation
  void setTime (double t)
  {
    time = t;
  }

  void setCalculateAll(bool calculateAll_)
  {
    calculateAll = calculateAll_;
  }



private:
  const GV& gv;
  Physics & physics;
  GF_CD& gfChargeDensity;
  RF time;
  ElementMapper elementMapper;
  RF xmax;
  RF xmin;
  bool calculateAll;
};
#endif /* DUNE_AX1_ACME0_INITIAL_HH */
