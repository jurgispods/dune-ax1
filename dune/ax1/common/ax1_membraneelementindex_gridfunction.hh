/*
 * ax1_subdomainelementindex_gridfunction.hh
 *
 *  Created on: Jul 16, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_SUBDOMAINELEMENTINDEX_GRIDFUNCTION_HH
#define DUNE_AX1_SUBDOMAINELEMENTINDEX_GRIDFUNCTION_HH

template<typename GV, typename PHYSICS>
class Ax1MembraneElementIndexGridFunction

  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::IntersectionGridFunctionTraits<GV, int, 4,
                                       std::tuple<int, int, typename PHYSICS::FieldType, int> >,
     Ax1MembraneElementIndexGridFunction<GV, PHYSICS> >
{
  public:

    typedef Dune::PDELab::IntersectionGridFunctionTraits<GV, int, 4,
        std::tuple<int, int, typename PHYSICS::FieldType, int> > Traits;


    Ax1MembraneElementIndexGridFunction(const GV& gv_, const PHYSICS& physics_)
    : gv(gv_),
      physics(physics_)
    {}


    template<typename I>
    inline void evaluate (const I& is,
                          const typename Traits::DomainType& x,
                                typename Traits::RangeType& y) const
    {
      typename PHYSICS::ElementPointer ep = is.inside();

      if(not physics.isMembrane(*ep))
        DUNE_THROW(Dune::Exception, "Ax1MembraneElementIndexGridFunction must be called for an interface with an inside membrane element!");

      // subdomain element index on this processor
      int localSubdomainElementIndex = physics.getSubDomainElementIndex(*ep);
      int globalSubdomainElementIndex = localSubdomainElementIndex + physics.nOffset(0);

      // Hardcoded correction for overlap size (assume partitioning in x-direction only)
      if(physics.nOffset(0) > 0)
        globalSubdomainElementIndex -= gv.overlapSize(0);

      y = typename Traits::RangeType(globalSubdomainElementIndex,
            physics.getGroupIndex(*ep),
            physics.getPermittivity(*ep),
            gv.comm().rank());

    }

    const typename Traits::GridViewType& getGridView() const
    {
      return gv;
    }

private:
    const GV& gv;
    const PHYSICS& physics;
};

#endif /* DUNE_AX1_SUBDOMAINELEMENTINDEX_GRIDFUNCTION_HH */
