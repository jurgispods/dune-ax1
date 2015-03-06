#ifndef DUNE_AX1_CUSTOM_POTENTIAL_HH
#define DUNE_AX1_CUSTOM_POTENTIAL_HH

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/acme2/common/acme2_parametertree.hh>

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

    CustomBoundaryPotential(const GV& gv, const Acme2Parameters& params_)
    : BaseT(gv),
      params(params_),
      potValues(4, 0.0)
    {}

    inline void evaluateGlobal(const DomainType & x, RangeType & y) const
    {
      // Zero on the interior
      y[0] = 0.;

      // Upper boundary
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
    }

    /*
    inline const GV& getGridView () const
    {
    }
    */

    // set time for subsequent evaluation
    void setTime (double t)
    {
    }

    void setPotValues(const std::vector<RField>& potValues_)
    {
      potValues = potValues_;
    }

  private:
    const Acme2Parameters& params;
    std::vector<RField> potValues;
};

#endif /* DUNE_AX1_CUSTOM_POTENTIAL_HH */
