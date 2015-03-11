/*
 * membranepotentialgridfunction.hh
 *
 *  Created on: Feb 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_MEMBRANEPOTENTIALGRIDFUNCTION_HH
#define DUNE_AX1_MEMBRANEPOTENTIALGRIDFUNCTION_HH

#include <dune/ax1/common/constants.hh>


template<typename DGF_POT, typename PHYSICS>
class MembranePotentialGridFunction
  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::IntersectionGridFunctionTraits<typename DGF_POT::Traits::GridViewType,
                                       typename DGF_POT::Traits::RangeFieldType,
                                       1,
                                       Dune::FieldVector<typename DGF_POT::Traits::RangeFieldType,1> >,
      MembranePotentialGridFunction<DGF_POT, PHYSICS> >
{
public:
  typedef Dune::PDELab::IntersectionGridFunctionTraits<typename DGF_POT::Traits::GridViewType,   // grid view type
                                           typename DGF_POT::Traits::RangeFieldType, // image field type (double)
                                           1,                                        // number of components of image (1)
                                           Dune::FieldVector<typename DGF_POT::Traits::RangeFieldType,1> // image type (Dune::FieldVector<double, 1>)
                                           > Traits;

  //! constructor
  MembranePotentialGridFunction (DGF_POT& dgfPot_, PHYSICS& physics_)
    : dgfPot(dgfPot_), gridView(dgfPot_.getGridView()), physics(physics_)
  {
    /*
    typename Traits::RangeType rt;
    typename Traits::RangeFieldType rf;
    typename Traits::DomainType dt;
    typename Traits::DomainFieldType df;
    typename DGF_CON::Traits::RangeType dgfConRT;

    debug_verb << "===== ChargeDensityGridFunction TYPEINFO ====" << std::endl;
    debug_verb << "Traits::RangeType = " << Tools::getTypeName(rt) << std::endl;
    debug_verb << "Traits::RangeFieldType = " << Tools::getTypeName(rf) << std::endl;
    debug_verb << "Traits::DomainType = " << Tools::getTypeName(dt) << std::endl;
    debug_verb << "Traits::DomainFieldType = " << Tools::getTypeName(df) << std::endl;
    debug_verb << "DGF_CON::Traits::RangeType = " << Tools::getTypeName(dgfConRT) << std::endl<< std::endl;
    */
  }

  template<typename I>
  inline void evaluate (const I &is,
                          const typename Traits::DomainType &x,
                          typename Traits::RangeType &y) const
  {
    y = 0.0;

    if(physics.isMembraneInterface(is))
    {
      physics.getMembranePotential(is, dgfPot, y);
      //debug_jochen << "Evaluate memb_pot @" << is.geometry().center() << ": " << y << std::endl;
    }
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  DGF_POT& dgfPot;
  const typename Traits::GridViewType& gridView;
  PHYSICS& physics;
};


template<typename SubGV, typename DGF_POT, typename PHYSICS>
class MembranePotentialGridFunction_MembGV
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<SubGV,
                                       typename DGF_POT::Traits::RangeFieldType,
                                       1,
                                       Dune::FieldVector<typename DGF_POT::Traits::RangeFieldType,1> >,
      MembranePotentialGridFunction_MembGV<SubGV,DGF_POT,PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<SubGV,   // grid view type
                                           typename DGF_POT::Traits::RangeFieldType, // image field type (double)
                                           1,                                        // number of components of image (1)
                                           Dune::FieldVector<typename DGF_POT::Traits::RangeFieldType,1> // image type (Dune::FieldVector<double, 1>)
                                           > Traits;

  //! constructor
  MembranePotentialGridFunction_MembGV (const SubGV& subGV, const DGF_POT& dgfPot_, const PHYSICS& physics_)
    : gridView(subGV),
      dgfPot(dgfPot_),
      physics(physics_)
  {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    assert(physics.isMembrane(e));

    physics.getMembranePotential(gridView.grid().multiDomainEntity(e), dgfPot, y);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  const typename Traits::GridViewType& gridView;
  const DGF_POT& dgfPot;
  const PHYSICS& physics;
};

#endif /* DUNE_AX1_MEMBRANEPOTENTIALGRIDFUNCTION_HH */
