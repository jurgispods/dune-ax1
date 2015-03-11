/*
 * ax1_conductivity_gridfunction.hh
 *
 *  Created on: Oct 16, 2014
 *      Author: jpods
 */

#ifndef DUNE_AX1_CONDUCTIVITY_GRIDFUNCTION_HH
#define DUNE_AX1_CONDUCTIVITY_GRIDFUNCTION_HH


#include <dune/ax1/common/constants.hh>


template<typename PARAM_POT, typename DGF_CON, typename PHYSICS>
class Ax1ConductivityGridFunction
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename DGF_CON::Traits::GridViewType,
             typename DGF_CON::Traits::RangeFieldType,
             1,
             Dune::FieldVector<typename DGF_CON::Traits::RangeFieldType,1> >,
             Ax1ConductivityGridFunction<PARAM_POT, DGF_CON, PHYSICS> >
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
  Ax1ConductivityGridFunction (PARAM_POT& param_, DGF_CON& dgfCon_, PHYSICS& physics_)
    : param(param_),
      dgfCon(dgfCon_),
      gridView(dgfCon_.getGridView()),
      physics(physics_),
      lengthScale(physics.getLengthScale()),
      timeScale(physics.getTimeScale()),
      temp(physics.getElectrolyte().getTemperature())
  {
  }

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    typename DGF_CON::Traits::RangeType conc(0.0);
    dgfCon.evaluate(e,x,conc);

    typename PARAM_POT::Traits::PermTensorType A = param.A(e, x, conc);

    // Assume isotropic conductivity here
    y = A[0][0];

    // MoriPoissonParameters::A returns a value with units (TS/LS) * (A/m)
    // => bring to S/m by multiplying by (LS/TS) and (e/(kT)) (units 1/V)
    y *= (lengthScale / timeScale) * con_e / (con_k * temp);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  PARAM_POT& param;
  DGF_CON& dgfCon;
  const typename Traits::GridViewType& gridView;
  PHYSICS& physics;
  typename Traits::RangeType lengthScale;
  typename Traits::RangeType timeScale;
  typename Traits::RangeType temp;
};



#endif /* DUNE_AX1_CONDUCTIVITY_GRIDFUNCTION_HH */
