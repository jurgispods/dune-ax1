/*
 * ax1_parameters_to_membraneflux_adapter.hh
 *
 *  Created on: Oct 7, 2014
 *      Author: jpods
 */

#ifndef DUNE_AX1_PARAMETERS_TO_MEMBRANEFLUX_ADAPTER_HH
#define DUNE_AX1_PARAMETERS_TO_MEMBRANEFLUX_ADAPTER_HH

template<typename PARAM_CON, int numSpecies = NUMBER_OF_SPECIES>
class Ax1ParametersToMembraneFluxAdapter
  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::IntersectionGridFunctionTraits<typename PARAM_CON::Traits::GridViewType,
                                       typename PARAM_CON::Traits::RangeFieldType,
                                       numSpecies,
                                       Dune::FieldVector<typename PARAM_CON::Traits::RangeFieldType,numSpecies> >,
     Ax1ParametersToMembraneFluxAdapter<PARAM_CON> >
{
public:
  typedef Dune::PDELab::IntersectionGridFunctionTraits<typename PARAM_CON::Traits::GridViewType,   // grid view type
                                           typename PARAM_CON::Traits::RangeFieldType, // image field type (double)
                                           numSpecies,  // number of components of image
                                           Dune::FieldVector<typename PARAM_CON::Traits::RangeFieldType,numSpecies> // image type (Dune::FieldVector<double, 1>)
                                           > Traits;

  Ax1ParametersToMembraneFluxAdapter(const PARAM_CON& paramCon_, const bool withoutMoriFlux_ = false)
  : paramCon(paramCon_),
    withoutMoriFlux(withoutMoriFlux_)
  {
    assert(numSpecies <= NUMBER_OF_SPECIES);
  }



  template<typename I>
  inline void evaluate (const I &is,
                          const typename Traits::DomainType &x,
                          typename Traits::RangeType &y) const
  {
    // TODO Multiply by physical constants to get the real current in Amps!
    for(int i=0; i<numSpecies; i++)
    {
      //debug_jochen << "Evaluating NernstPlanckParameters::memb_flux on intersection @ "
      //    << is.geometry().center() << std::endl;

      //If flag 'withoutMoriFlux' is explicitly set to true, call j with additional bool parameter
      // (only available in Mori parameter classes)
      if(withoutMoriFlux)
      {
        y[i] = paramCon[i]->j(is, x, false);
      } else {
        y[i] = paramCon[i]->j(is, x);
      }
    }
  }


  inline const typename Traits::GridViewType& getGridView () const
  {
    return paramCon[0]->getGridView();
  }

private:
  const PARAM_CON& paramCon;
  const bool withoutMoriFlux;
};



#endif /* DUNE_AX1_PARAMETERS_TO_MEMBRANEFLUX_ADAPTER_HH */
