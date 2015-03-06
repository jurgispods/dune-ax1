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
  : paramCon(paramCon_),
    lengthScale(paramCon_[0]->getPhysics().getLengthScale()),
    timeScale(paramCon_[0]->getPhysics().getLengthScale())
  {
    assert(numSpecies <= NUMBER_OF_SPECIES);
  }

  template<typename I>
  inline void evaluate (const I &is,
                          const typename Traits::DomainType &x,
                          typename Traits::RangeType &y) const
  {
    for(int i=0; i<numSpecies; i++)
    {
      //debug_jochen << "Evaluating NernstPlanckParameters::memb_flux on intersection @ "
      //    << is.geometry().center() << std::endl;

      y[i] = paramCon[i]->j(is, x) * paramCon[0]->getPhysics().getValence(i);

    }

    // Get intersection geometry
    typedef typename Acme2CylGeometrySwitch::GeometrySwitch<typename I::Geometry>::type GEO;
    GEO geo(is.geometry());

    // Intersection volume is in units of LS^2, bring to SI units m^2
    typename Traits::RangeFieldType scale = geo.volume() * lengthScale * lengthScale;

    // Now add factors to bring membrane flux j (units (TS/LS) * mol/(m^2 * s)) to C/(m^2 * s) = A/m^2
    scale *= (lengthScale / timeScale) * con_e * con_mol;

    // Resulting membrane current in Amps
    y *= scale;
  }


  inline const typename Traits::GridViewType& getGridView () const
  {
    return paramCon[0]->getGridView();
  }

private:
  const PARAM_CON& paramCon;
  const double lengthScale;
  const double timeScale;
};

#endif /* DUNE_AX1_MEMBRANECURRENT_GRIDFUNCTION_HH */
