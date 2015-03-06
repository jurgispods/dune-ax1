/*
 * ax1_permittivity_gridfunction.hh
 *
 *  Created on: Jul 8, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_PERMITTIVITY_GRIDFUNCTION_HH
#define DUNE_AX1_PERMITTIVITY_GRIDFUNCTION_HH

template<typename GV, typename PHYSICS>
class Ax1PermittivityGridFunction
  :  public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV, typename PHYSICS::FieldType, 1,
                                                                            Dune::FieldVector<typename PHYSICS::FieldType,1> >,
     Ax1PermittivityGridFunction<GV,PHYSICS> >
{
  typedef Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV, typename PHYSICS::FieldType, 1,
                                         Dune::FieldVector<typename PHYSICS::FieldType,1> >,
      Ax1PermittivityGridFunction<GV,PHYSICS> > BaseT;

public:
    typedef typename BaseT::Traits Traits;

    Ax1PermittivityGridFunction(const GV& gv_, const PHYSICS& physics_)
    : gv(gv_),
      physics(physics_)
    {}


    // Evalute grid function; return 0 if the element is part of the membrane!
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                                typename Traits::RangeType& y) const
    {
      y = physics.getPermittivity(e);
    }

    const typename Traits::GridViewType& getGridView() const
    {
      return gv;
    }

private:
    const GV& gv;
    const PHYSICS& physics;
};


#endif /* DUNE_AX1_PERMITTIVITY_GRIDFUNCTION_HH */
