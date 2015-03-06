#ifndef DUNE_AX1_ACME1_BOUNDARY_HH
#define DUNE_AX1_ACME1_BOUNDARY_HH

#include <dune/common/fvector.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/constraints/constraintsparameters.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>

/*! Adapter that extracts boundary condition type function from parameter class
 *  \tparam T  model of ConvectionDiffusionParameterInterface
 */
template<typename T>
class BCTypeSingleCon : public Dune::PDELab::DirichletConstraintsParameters
{
  public:
    BCTypeSingleCon(const typename T::Traits::GridViewType& gv_, T& t_, const int ionSpecies_ )
      : gv( gv_ ), t( t_ ), ionSpecies(ionSpecies_)
    {}

    template<typename I>
    bool isDirichlet(const I& ig, const Dune::FieldVector<typename I::ctype, I::dimension-1>& coord) const
    {
      t.setIonSpecies(ionSpecies);
      return(t.bctype(ig.intersection(), coord) == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet);
    }

    int getIonSpecies()
    {
      return ionSpecies;
    }
  
private:
  const typename T::Traits::GridViewType& gv;
  T& t;
  const int ionSpecies;
};



/*! Adapter that extracts Dirichlet boundary conditions from parameter class
 *  \tparam T  model of ConvectionDiffusionParameterInterface
 */
template<typename T>
class DirichletValuesSingleCon
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits
      <typename T::Traits::GridViewType,
       typename T::Traits::RangeFieldType,
       1,
       Dune::FieldVector<typename T::Traits::RangeFieldType,1>
      >,
    DirichletValuesSingleCon<T> >
{
  public:
    typedef Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                             typename T::Traits::RangeFieldType,
                                             1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> > Traits;

    //! constructor
    DirichletValuesSingleCon(const typename Traits::GridViewType& g_, T& t_, const int ionSpecies_)
    : g(g_), t(t_), ionSpecies(ionSpecies_)
    {}

    //! \copydoc GridFunctionBase::evaluate()
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                          typename Traits::RangeType& y) const
    {
      t.setIonSpecies(ionSpecies);
      y = t.g(e,x);
    }

    inline const typename Traits::GridViewType& getGridView () const
    {
      return g;
    }

    inline void setTime(double time_)
    {
      t.setTime(time_);
    }

  private:
    const typename Traits::GridViewType& g;
    T& t;
    const int ionSpecies;
};

#endif /* DUNE_AX1_ACME1_BOUNDARY_HH */
