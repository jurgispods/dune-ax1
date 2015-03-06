/*
 * ax1_multidomaingridfunction.hh
 *
 *  Created on: Mar 23, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_MULTIDOMAINGRIDFUNCTIONRESTRICTION_HH
#define DUNE_AX1_MULTIDOMAINGRIDFUNCTIONRESTRICTION_HH

#include <dune/ax1/common/constants.hh>


/**
 * This class evaluates a given grid function living on the multidomain grid only of the specified subdomain GV
 */
template<typename GV, typename DGF, typename PHYSICS>
class Ax1MultiDomainGridFunctionRestriction
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,
               typename DGF::Traits::RangeFieldType,
               DGF::Traits::RangeType::dimension,
               typename DGF::Traits::RangeType>,
               Ax1MultiDomainGridFunctionRestriction<GV, DGF, PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,                           // grid view type
               typename DGF::Traits::RangeFieldType,                     // range field type (double)
               DGF::Traits::RangeType::dimension,                        // number of components of image (1)
               typename DGF::Traits::RangeType                           // image type (Dune::FieldVector<double, 1>)
               > Traits;

  //typedef typename DGF::Traits Traits;

  typedef typename DGF::Traits::RangeFieldType RF;
  typedef typename DGF::Traits::RangeType RT;


  //! constructor
  Ax1MultiDomainGridFunctionRestriction (const GV& gv_, DGF& dgf_, PHYSICS& physics_)
    : gridView(gv_), dgf(dgf_), physics(physics_)
  {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    assert(dgf.getGridView().contains(gridView.grid().multiDomainEntity(e)));

    dgf.evaluate(gridView.grid().multiDomainEntity(e),x,y);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  const GV& gridView;
  DGF& dgf;
  PHYSICS& physics;
};


#endif /* DUNE_AX1_MULTIDOMAINGRIDFUNCTIONRESTRICTION_HH */
