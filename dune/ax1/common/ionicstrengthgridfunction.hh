/*
 * ionicstrengthgridfunction.hh
 *
 *  Created on: 13.12.2011
 *      Author: hannes
 */

#ifndef DUNE_AX1_IONICSTRENGTHGRIDFUNCTION_HH_
#define DUNE_AX1_IONICSTRENGTHGRIDFUNCTION_HH_

#include <dune/ax1/common/constants.hh>


template<typename DGF_CON, typename PHYSICS>
class IonicStrengthGridFunction
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename DGF_CON::Traits::GridViewType,
                                                                           typename DGF_CON::RF,
                                                                           1,
                                                                           Dune::FieldVector<typename DGF_CON::RF,1> >
                                          ,IonicStrengthGridFunction<DGF_CON, PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename DGF_CON::Traits::GridViewType,   // grid view type
                                           typename DGF_CON::RF,                     // range field type (double)
                                           1,                                        // number of components of image (1)
                                           Dune::FieldVector<typename DGF_CON::RF,1> // image type (Dune::FieldVector<double, 1>)
                                           > Traits;

  //! constructor
  IonicStrengthGridFunction (DGF_CON& dgfCon_, PHYSICS& physics_)
    : dgfCon(dgfCon_), gridView(dgfCon_.getGridView()), physics(physics_)
  {
    /*
    typename Traits::RangeType rt;
    typename Traits::RangeFieldType rf;
    typename Traits::DomainType dt;
    typename Traits::DomainFieldType df;
    typename DGF_CON::Traits::RangeType dgfConRT;

    debug_verb << "===== IonicStrengthGridFunction TYPEINFO ====" << std::endl;
    debug_verb << "Traits::RangeType = " << Tools::getTypeName(rt) << std::endl;
    debug_verb << "Traits::RangeFieldType = " << Tools::getTypeName(rf) << std::endl;
    debug_verb << "Traits::DomainType = " << Tools::getTypeName(dt) << std::endl;
    debug_verb << "Traits::DomainFieldType = " << Tools::getTypeName(df) << std::endl;
    debug_verb << "DGF_CON::Traits::RangeType = " << Tools::getTypeName(dgfConRT) << std::endl<< std::endl;
    */
  }

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    typename DGF_CON::Traits::RangeType conc;
    dgfCon.evaluate(e,x,conc);
    for(int j=0; j<NUMBER_OF_SPECIES; ++j)
    {
      y += 0.5 * physics.getValence(j) * physics.getValence(j) * conc[j];
    }
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  DGF_CON& dgfCon;
  const typename Traits::GridViewType& gridView;
  PHYSICS& physics;
};

#endif /* DUNE_AX1_IONICSTRENGTHGRIDFUNCTION_HH_ */
