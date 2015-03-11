/*
 * ax1_multidomaingridfunction.hh
 *
 *  Created on: Mar 23, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_MULTIDOMAINGRIDFUNCTIONEXTENSION_HH
#define DUNE_AX1_MULTIDOMAINGRIDFUNCTIONEXTENSION_HH

#include <dune/ax1/common/constants.hh>

/**
 * This class puzzles together the grid function on two subdomains to yield a grid function on the
 * multidomain grid
 */
template<typename GV, typename DGF, typename PHYSICS>
class Ax1MultiDomainGridFunctionExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,
               typename DGF::Traits::RangeFieldType,
               DGF::Traits::RangeType::dimension,
               typename DGF::Traits::RangeType
               >,
               Ax1MultiDomainGridFunctionExtension<GV, DGF, PHYSICS> >
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
  typedef DGF BASE_DGF;


  //! constructor
  Ax1MultiDomainGridFunctionExtension (const GV& gv_, const DGF& dgf_, const PHYSICS& physics_)
    : gridView(gv_), dgf(dgf_), physics(physics_)
  {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    if(not physics.isMembrane(e))
    {
      typename PHYSICS::SubDomainElementPointer sdep = dgf.getGridView().grid().subDomainEntityPointer(e);
      assert(dgf.getGridView().contains(*sdep));

      dgf.evaluate(*sdep,x,y);
    }
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

  void setTime(double time)
  {
    dgf.setTime(time);
  }

private:
  const GV& gridView;
  const DGF& dgf;
  const PHYSICS& physics;

};


#endif /* DUNE_AX1_MULTIDOMAINGRIDFUNCTIONEXTENSION_HH */
