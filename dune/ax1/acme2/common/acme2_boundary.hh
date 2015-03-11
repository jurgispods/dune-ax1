#ifndef DUNE_AX1_ACME2_BOUNDARY_HH
#define DUNE_AX1_ACME2_BOUNDARY_HH

#include <dune/common/fvector.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/constraints/constraintsparameters.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>

/*! Adapter that extracts boundary condition type function from parameter class
 *  \tparam PARAMS  model of ConvectionDiffusionParameterInterface
 */
template<typename PARAMS, typename PHYSICS>
class BCTypeSingleCon : public Dune::PDELab::DirichletConstraintsParameters
{
  public:
    BCTypeSingleCon(const typename PARAMS::Traits::GridViewType& gv_, PARAMS& params_, PHYSICS& physics_, const int ionSpecies_ )
      : gv( gv_ ), params( params_ ), physics(physics_), ionSpecies(ionSpecies_)
    {}

    template<typename I>
    bool isDirichlet(const I& ig, const Dune::FieldVector<typename I::ctype, I::dimension-1>& coord) const
    {
      //t.setIonSpecies(ionSpecies);
      return(params.bctype(ig.intersection(), coord, ionSpecies) == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet);
    }

    template<typename I>
    bool isNeumann(const I& ig, const Dune::FieldVector<typename I::ctype, I::dimension-1>& coord) const
    {
      //t.setIonSpecies(ionSpecies);
      return(params.bctype(ig.intersection(), coord, ionSpecies) == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
    }

    int getIonSpecies()
    {
      return ionSpecies;
    }
  
private:
  const typename PARAMS::Traits::GridViewType& gv;
  PARAMS& params;
  PHYSICS& physics;
  const int ionSpecies;
};


/*! Adapter that extracts boundary condition type function from parameter class
 *  \tparam PARAMS  model of ConvectionDiffusionParameterInterface
 */
template<typename PARAMS, typename PHYSICS>
class BCTypePot : public Dune::PDELab::DirichletConstraintsParameters   /*@\label{bcp:base}@*/
{

public:

  BCTypePot(const typename PARAMS::Traits::GridViewType& gv_, const PARAMS& params_, PHYSICS& physics_ )
    : gv( gv_ ), params( params_ ), physics(physics_)
  {
  }

  template<typename I>
  bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1>& coord) const
  {
    // No boundary condition on the membrane
    //if(physics.isMembraneInterface(ig.intersection())) return false;

    return( params.bctype( ig.intersection(), coord )
            == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet );
  }

  template<typename I>
  bool isNeumann(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord) const
  {
    // No boundary condition on the membrane
    //if(physics.isMembraneInterface(ig.intersection())) return false;

    return( params.bctype( ig.intersection(), coord )
        == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann );
  }

private:
  const typename PARAMS::Traits::GridViewType& gv;
  const PARAMS& params;
  PHYSICS& physics;

};


/*! Adapter that extracts Dirichlet boundary conditions from parameter class
 *  \tparam PARAMS  model of ConvectionDiffusionParameterInterface
 */
template<typename PARAMS>
class DirichletValuesSingleCon
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits
      <typename PARAMS::Traits::GridViewType,
       typename PARAMS::Traits::RangeFieldType,
       1,
       Dune::FieldVector<typename PARAMS::Traits::RangeFieldType,1>
      >,
    DirichletValuesSingleCon<PARAMS> >
{
  public:
    typedef Dune::PDELab::GridFunctionTraits<typename PARAMS::Traits::GridViewType,
                                             typename PARAMS::Traits::RangeFieldType,
                                             1,Dune::FieldVector<typename PARAMS::Traits::RangeFieldType,1> > Traits;

    //! constructor
    DirichletValuesSingleCon(const typename Traits::GridViewType& g_, PARAMS& t_, const int ionSpecies_)
    : g(g_), params(t_), ionSpecies(ionSpecies_)
    {}

    //! \copydoc GridFunctionBase::evaluate()
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                          typename Traits::RangeType& y) const
    {
      y = params.g(e,x,ionSpecies);
    }

    inline const typename Traits::GridViewType& getGridView () const
    {
      return g;
    }

    inline void setTime(double time_)
    {
      params.setTime(time_);
    }

  private:
    const typename Traits::GridViewType& g;
    PARAMS& params;
    const int ionSpecies;
};

#endif /* DUNE_AX1_ACME2_BOUNDARY_HH */
