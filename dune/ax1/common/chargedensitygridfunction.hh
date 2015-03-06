/*
 * chargedensitygridfunction.hh
 *
 *  Created on: Sep 5, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_CHARGEDENSITYGRIDFUNCTION_HH
#define DUNE_AX1_CHARGEDENSITYGRIDFUNCTION_HH

#include <dune/ax1/common/constants.hh>


template<typename DGF_CON, typename PHYSICS>
class ChargeDensityGridFunction
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename DGF_CON::Traits::GridViewType,
             typename DGF_CON::Traits::RangeFieldType,
             1,
             Dune::FieldVector<typename DGF_CON::Traits::RangeFieldType,1> >,
             ChargeDensityGridFunction<DGF_CON, PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename DGF_CON::Traits::GridViewType,   // grid view type
             typename DGF_CON::Traits::RangeFieldType,                     // range field type (double)
             1,                                                            // number of components of image (1)
             Dune::FieldVector<typename DGF_CON::Traits::RangeFieldType,1> // image type (Dune::FieldVector<double, 1>)
             > Traits;

  //! \brief Flag to decide if an extrapolation of the concentrations should be done
  const static bool doExtrapolation = false;

  //! constructor
  ChargeDensityGridFunction (DGF_CON& dgfCon_, DGF_CON& dgfConPrevious_, PHYSICS& physics_)
    : dgfCon(dgfCon_), dgfConPrevious(dgfConPrevious_), gridView(dgfCon_.getGridView()), physics(physics_)
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

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    typename DGF_CON::Traits::RangeType conc(0.0);

    //Dune::ios_base_all_saver asd(std::cout);
    //if(physics.getElementIndex(e) == 0)
    //{
    //  dgfCon.evaluate(e,0,conc);
    //  debug_jochen << "[" << ION_NAMES[0] << "] @ " << e.geometry().global(0) << " = " << std::scientific << std::setprecision(16)
    //    << conc[0] << std::endl;
    //}

    dgfCon.evaluate(e,x,conc);

//    typename DGF_CON::Traits::RangeType concExtrapolated(conc);
//    // Extrapolate new concentrations; assume one step time stepper with constant dt here!
//    if(doExtrapolation)
//    {
//      typename DGF_CON::Traits::RangeType concPrevious;
//      dgfConPrevious.evaluate(e,x,concPrevious);
//      concExtrapolated += conc;
//      concExtrapolated -= concPrevious;
//    }
//    conc = concExtrapolated;

    //debug_jochen << " - ";
    for(int j=0; j<NUMBER_OF_SPECIES; ++j)
    {
      y += physics.getValence(j) * conc[j];
      //debug_jochen << std::scientific << std::setprecision(30) << conc[j] << " + ";
    }
    //debug_jochen << " -> " << y << std::endl;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  DGF_CON& dgfCon;
  DGF_CON& dgfConPrevious;
  const typename Traits::GridViewType& gridView;
  PHYSICS& physics;
};

#endif /* DUNE_AX1_CHARGEDENSITYGRIDFUNCTION_HH */
