/*
 * ax1_membranecurrent_gridfunction.hh
 *
 *  Created on: Feb 18, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_MEMBRANECURRENT_GRIDFUNCTION_HH
#define DUNE_AX1_MEMBRANECURRENT_GRIDFUNCTION_HH

template<typename PARAM_CON, int numSpecies = NUMBER_OF_SPECIES>
class Ax1MembraneCurrentGridFunction
  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::IntersectionGridFunctionTraits<typename PARAM_CON::Traits::GridViewType,
                                       typename PARAM_CON::Traits::RangeFieldType,
                                       numSpecies,
                                       Dune::FieldVector<typename PARAM_CON::Traits::RangeFieldType,numSpecies> >,
     Ax1MembraneCurrentGridFunction<PARAM_CON> >
{
public:
  typedef Dune::PDELab::IntersectionGridFunctionTraits<typename PARAM_CON::Traits::GridViewType,   // grid view type
                                           typename PARAM_CON::Traits::RangeFieldType, // image field type (double)
                                           numSpecies,  // number of components of image
                                           Dune::FieldVector<typename PARAM_CON::Traits::RangeFieldType,numSpecies> // image type (Dune::FieldVector<double, 1>)
                                           > Traits;

  Ax1MembraneCurrentGridFunction(const PARAM_CON& paramCon_)
  : paramCon(paramCon_)
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

      y[i] = paramCon[i]->memb_flux(is, x);
    }
  }


  inline const typename Traits::GridViewType& getGridView () const
  {
    return paramCon[0]->getGridView();
  }

private:
  const PARAM_CON& paramCon;
};

#endif /* DUNE_AX1_MEMBRANECURRENT_GRIDFUNCTION_HH */
