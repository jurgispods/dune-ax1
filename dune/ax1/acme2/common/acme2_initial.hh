#ifndef DUNE_AX1_ACME2_INITIAL_HH
#define DUNE_AX1_ACME2_INITIAL_HH

#include <valarray>

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/acme2/common/acme2_parametertree.hh>
#include <dune/ax1/acme2/configurations/zero_potential.hh>

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
    /*
    bool isMembraneBoundary = (std::abs(x[1]-params.yMemb()) < 1e-12
        || std::abs(x[1] - params.yMemb() - params.dMemb()) < 1e-12);

    // No charge carriers inside the membrane!
    if(physics.isMembrane(e) && not isMembraneBoundary) // Insane in the membrane!
    {
      for(int i=0; i<numSpecies; ++i)
      {
        y[i] = 0.0;
      }
      return;
    }
    */

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
          DUNE_THROW(Dune::Exception, "Element has an unknown group index!");
      }
    }

    // Override default behaviour by evaluating custom intitial gridfunction
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
  }

private:
  PHYSICS physics;
  INITIAL_GF& initialGF;
  const GV& gv;
  const Acme2Parameters& params;
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

      InitialVoid(const GV& gv, const Acme2Parameters& params_)
      : BaseT(gv)
      {}

      // Does nothing
      inline void evaluateGlobal(const DomainType & x, RangeType & y) const
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

    InitialStep(const GV& gv, const Acme2Parameters& params_)
    : BaseT(gv),
      params(params_),
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
      for ( int i=0; i<numSpecies; ++i )
      {
        // Intracellular domain
        if (x[1] < params.yMemb() || std::abs(x[1] - params.yMemb()) < 1e-12)
        {
          y[i] = initConcIn[i];
        }
        // Extracellular domain
        else if(x[1] > params.yMemb() || std::abs(x[1] - params.yMemb()) < 1e-12)
        {
          y[i] = initConcEx[i];
        // Membrane
        } else {
          // Smooth discontinuity in initial concentrations
          if(smoothStep)
          {
            typename Traits::DomainFieldType y_scale = (x[1] - params.yMemb()) / (params.dMemb());
            y[i] = initConcIn[i] + (initConcEx[i]-initConcIn[i]) * (y_scale*y_scale*(3-2*y_scale));
          // No charge carriers inside membrane
          } else {
            if (x[1] < params.yMemb() + 0.5 * params.dMemb()) y[i] = initConcIn[i];
            else                                              y[i] = initConcEx[i];
          }
        }
      }
    }


  private:
    const Acme2Parameters& params;
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

    InitialSquarePulse(const GV& gv, const Acme2Parameters& params_)
    : BaseT(gv),
      params(params_)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      double dpulse = 5.0;
      double xCenter = 0.5 * (params.xMax() - params.xMin());
      double yCenter = 0.5 * (params.yMemb() - params.yMin());
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


  private:
    const Acme2Parameters& params;
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

    InitialTrianglePulse(const GV & gv, const Acme2Parameters& params_)
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

  private:
    const Acme2Parameters& params;
};

// initial potential ###############################################################################

template<typename GV, typename RF, typename Physics, typename GF_CD, typename INITIAL_GF>
class InitialPot
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
    InitialPot<GV,RF,Physics,GF_CD,INITIAL_GF> >
{

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  InitialPot (const GV& gv_, Physics& physics_, GF_CD& gfChargeDensity_, INITIAL_GF& initialGF_, int intorder_=2) :
    gv(gv_), physics(physics_), gfChargeDensity(gfChargeDensity_), initialGF(initialGF_), intorder(intorder_),
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
    
    //if(physics.getSubdomainIndex(e) == CYTOSOL)
    //{
    //  y = physics.convertFrom_mV(-65.);
    //} else {
      initialGF.evaluate(e,xlocal,y);
    //}

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



    y = 0.0;
    */

  }


  inline RF calculateExactPotential(const typename Traits::DomainType& x) const
  {
    const int dim = GV::dimension;
    assert(dim == 1);

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
  INITIAL_GF& initialGF;
  const int intorder;
  RF time;
  RF xmax;
  RF xmin;
  bool calculateAll;

};
#endif /* DUNE_AX1_ACME2_INITIAL_HH */
