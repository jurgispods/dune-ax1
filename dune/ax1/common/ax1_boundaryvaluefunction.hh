/*
 * ax1_boundaryvaluefunction.hh
 *
 *  Created on: May 28, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_BOUNDARYVALUEFUNCTION_HH
#define DUNE_AX1_BOUNDARYVALUEFUNCTION_HH



#include <dune/ax1/common/constants.hh>

/**
 * This is one of the central classes of the dune-ax1 simulator. It is by design a boundary gridfunction,
 * i.e. it is defined on the electrolyte-membrane interfaces and calculates trans-membrane flux
 * by iterating over all channels present on the corresponding membrane element. For each ion species,
 * the conductance with respect to the current potential jump at the membrane as well as the inner and
 * outer concentrations is calculated and cached. The calculation is done in the updateFlux method,
 * while the channel states are advanced one time step in the updateChannels method.
 * This cached raw flux (which is not allowed to change within one time step) is multiplied by the outer
 * normal and scaled by the area of the interface in the evaluate method to yield the total flux for
 * an ion species.
 */
template<typename GF, typename PHYSICS>
class Ax1BoundaryValueFunction
  : public Dune::PDELab::BoundaryGridFunctionBase<
      Dune::PDELab::BoundaryGridFunctionTraits<typename GF::Traits::GridViewType,
                                       typename GF::Traits::RangeFieldType,
                                       GF::Traits::RangeType::dimension,
                                       typename GF::Traits::RangeType>,
     Ax1BoundaryValueFunction<GF,PHYSICS> >
{
public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<typename GF::Traits::GridViewType,   // grid view type
                                           typename GF::Traits::RangeFieldType, // image field type (double)
                                           GF::Traits::RangeType::dimension,  // number of components of image
                                           typename GF::Traits::RangeType // image type (Dune::FieldVector<double, 1>)
                                           > Traits;

  //! constructor
  Ax1BoundaryValueFunction (GF& gf_, PHYSICS& physics_)
     : gf(gf_),
       gridView(gf_.getGridView()),
       physics(physics_)
   {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I> &ig,
                            const typename Traits::DomainType &x,
                            typename Traits::RangeType &y) const
  {
    // Map x (local coordinate of intersection) to x_entity (local coordinate of inside entity)
    typename GF::Traits::DomainType x_entity = ig.geometryInInside().global(x);

    gf.evaluate(*ig.inside(), ig.geometryInInside().lo, y);

    debug_jochen << "Evaluating boundary GF @ coordinate " << (*ig.inside()).geometry().global(x_entity)
      << ": y=" << y << std::endl;

  }


private:
  GF& gf;
  typename GF::Traits::GridViewType gridView;
  PHYSICS& physics;

};

#endif /* DUNE_AX1_BOUNDARYVALUEFUNCTION_HH */
