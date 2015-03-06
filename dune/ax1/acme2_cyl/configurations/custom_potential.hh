#ifndef DUNE_AX1_CUSTOM_POTENTIAL_HH
#define DUNE_AX1_CUSTOM_POTENTIAL_HH

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>

template<typename GV, typename RField, int dim>
class CustomBoundaryPotential :
  public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RField,dim>,
    CustomBoundaryPotential<GV,RField,dim> >
{
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RField,dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, CustomBoundaryPotential<GV,RField,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    typedef RField RF;
    typedef RangeType RT;

    CustomBoundaryPotential(const GV& gv, const Acme2CylParameters& params_)
    : BaseT(gv),
      params(params_),
      potValues(4, 0.0),
      time(0.0)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      assert(params.nMembranes() <= 1);
      // ES
      if(! params.useMembrane() || x[1] + 1e-6 > params.yMemb()[0] + params.dMemb())
      {
        y[0] = convertFrom_mV(params.solution_ex.get("pot", 0.0));
        //debug_jochen << ", g= " << y << std::endl;
        return;
      }
      // CY
      if(x[1] - 1e-6 < params.yMemb()[0] )
      {
        y[0] = convertFrom_mV(params.solution_in.get("pot", 0.0));
        //debug_jochen << ", g= " << y << std::endl;
        return;
      }

      // Case membrane element at time other than initial condition
      if(time > 0.0)
      {
        DUNE_THROW(Dune::Exception, "Element is neither ES nor CYTOSOL!");
      }

      // Membrane at initial time=0
      y[0] = convertFrom_mV(params.solution_ex.get("pot", 0.0));
      //debug_jochen << ", g= " << y << std::endl;
      return;

      /*
       * Deactivate this crap for now. The whole initialization structure will be simplified when
       * incorporating everything into one single parameter class. The following initialization was
       * never really used at all, and it wasn't very useful anyway.
      // Lower boundary
      if(x[1]-1e-6 < params.yMin())
      {
        y[0] = potValues[0];
      }
      // Left boundary
      if(x[0]-1e-6 < params.xMin())
      {
        y[0] = potValues[1];
      }
      // Right boundary
      if(x[0]+1e-6 > params.xMax())
      {
        y[0] = potValues[2];
      }
      // Upper boundary
      if(x[1]+1e-6 > params.yMax())
      {
        y[0] = potValues[3];
      }
      */
    }

    /*
    inline const GV& getGridView () const
    {
    }
    */

    // set time for subsequent evaluation
    void setTime (double t)
    {
      time = t;
    }

    void setPotValues(const std::vector<RField>& potValues_)
    {
      potValues = potValues_;
    }

  private:
    RF convertFrom_mV(RF pot) const
    {
      RF y = pot * (con_e / (1.0e3 * con_k * 279.45));
      return y;
    }

    const Acme2CylParameters& params;
    std::vector<RField> potValues;
    double time;
};

#endif /* DUNE_AX1_CUSTOM_POTENTIAL_HH */
