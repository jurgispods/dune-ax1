/*
 * chargedensitygridfunction.hh
 *
 *  Created on: Sep 7, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_IONFLUXGRIDFUNCTION_HH
#define DUNE_AX1_IONFLUXGRIDFUNCTION_HH

#include <dune/ax1/common/constants.hh>


template<typename DGF_CON, typename DGF_CON_GRAD, typename DGF_POT_GRAD, typename PHYSICS>
class IonFluxGridFunction
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename DGF_CON_GRAD::Traits::GridViewType,
             //typename DGF_CON_GRAD::RF,
             typename DGF_CON_GRAD::Traits::RangeFieldType,
             DGF_CON_GRAD::Traits::RangeType::dimension,
             typename DGF_CON_GRAD::Traits::RangeType >,
             IonFluxGridFunction<DGF_CON, DGF_CON_GRAD, DGF_POT_GRAD, PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename DGF_CON_GRAD::Traits::GridViewType,   // grid view type
             //typename DGF_CON_GRAD::RF,
             typename DGF_CON_GRAD::Traits::RangeFieldType,  // range field type (double)
             DGF_CON_GRAD::Traits::RangeType::dimension,      // dimension of image
             typename DGF_CON_GRAD::Traits::RangeType > Traits;

  //! constructor
  IonFluxGridFunction (DGF_CON& dgfCon_, DGF_CON_GRAD& dgfConGrad_, DGF_POT_GRAD& dgfPotGrad_, PHYSICS& physics_)
    : dgfCon(dgfCon_),
      dgfConGrad(dgfConGrad_),
      dgfPotGrad(dgfPotGrad_),
      gridView(dgfCon_.getGridView()),
      physics(physics_)
  {
    /*
    typename Traits::RangeType rt;
    typename Traits::RangeFieldType rf;
    typename Traits::DomainType dt;
    typename Traits::DomainFieldType df;
    typename DGF_POT_GRAD::Traits::RangeType dgfPotGradRT;
    typename DGF_POT_GRAD::Traits::RangeFieldType dgfPotGradRFT;
    //typename DGF_POT_GRAD::RF dgfPotGradRF;
    typename DGF_CON::Traits::RangeType dgfConRT;
    typename DGF_CON::Traits::RangeFieldType dgfConRFT;
    typename DGF_CON::RF dgfConRF;
    typename DGF_CON_GRAD::Traits::RangeType dgfConGradRT;
    typename DGF_CON_GRAD::Traits::RangeFieldType dgfConGradRFT;
    typename DGF_CON_GRAD::RF dgfConGradRF;

    debug_verb << "===== ChargeDensityGridFunction TYPEINFO ====" << std::endl;
    debug_verb << "Traits::RangeType = " << Tools::getTypeName(rt) << std::endl;
    debug_verb << "Traits::RangeFieldType = " << Tools::getTypeName(rf) << std::endl;
    debug_verb << "Traits::DomainType = " << Tools::getTypeName(dt) << std::endl;
    debug_verb << "Traits::DomainFieldType = " << Tools::getTypeName(df) << std::endl;
    debug_verb << "DGF_POT_GRAD::Traits::RangeType = " << Tools::getTypeName(dgfPotGradRT) << std::endl;
    debug_verb << "DGF_POT_GRAD::Traits::RangeFieldType = " << Tools::getTypeName(dgfPotGradRFT) << std::endl;
    //debug_verb << "DGF_POT_GRAD::RF = " << Tools::getTypeName(dgfPotGradRF) << std::endl<< std::endl;
    debug_verb << "DGF_CON::Traits::RangeType = " << Tools::getTypeName(dgfConRT) << std::endl;
    debug_verb << "DGF_CON::Traits::RangeFieldType = " << Tools::getTypeName(dgfConRFT) << std::endl;
    debug_verb << "DGF_CON::RF = " << Tools::getTypeName(dgfConRF) << std::endl<< std::endl;
    debug_verb << "DGF_CON_GRAD::Traits::RangeType = " << Tools::getTypeName(dgfConGradRT) << std::endl;
    debug_verb << "DGF_CON_GRAD::Traits::RangeFieldType = " << Tools::getTypeName(dgfConGradRFT) << std::endl;
    debug_verb << "DGF_CON_GRAD::RF = " << Tools::getTypeName(dgfConGradRF) << std::endl<< std::endl;
  */
  }

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    int elemIndex = physics.getElementIndex(e);
    y = 0.0;

    typename DGF_CON::Traits::RangeType conc;
    typename DGF_CON_GRAD::Traits::RangeType concGrad;
    typename DGF_POT_GRAD::Traits::RangeType potGrad;

    dgfCon.evaluate(e,x,conc);
    dgfConGrad.evaluate(e,x,concGrad);
    dgfPotGrad.evaluate(e,x,potGrad);

    //debug_jochen << "°°° element @ " << e.geometry().corner(0) << std::endl;

    const int dimDomain = DGF_CON_GRAD::ChildType::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimDomain;
    const int nComponents = y.size();

    for(int i=0; i<nComponents; ++i)
    {
      int j = i / dimDomain;

      y[i] = concGrad[i];
      y[i] += conc[j] * physics.getValence(j) * potGrad[i % dimDomain];

      // Minus!!
      y[i] *= - physics.getDiffCoeff(j, elemIndex);

      std::string direction = (i % dimDomain == 0) ? "x" : "y";

      //debug_jochen << ION_NAMES[j] << " flux (" << direction << ") = - "
      //    << physics.getDiffCoeff(j, elemIndex) << " * ("
      //    << concGrad[i] << " + " << conc[j] << " * " << physics.getValence(j) << " * "
      //    << potGrad[i % dimDomain]
      //    << ") = " << y[i] << std::endl;
    }
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  DGF_CON& dgfCon;
  DGF_CON_GRAD& dgfConGrad;
  DGF_POT_GRAD& dgfPotGrad;
  const typename Traits::GridViewType& gridView;
  PHYSICS& physics;
};

#endif /* DUNE_AX1_IONFLUXGRIDFUNCTION_HH */
