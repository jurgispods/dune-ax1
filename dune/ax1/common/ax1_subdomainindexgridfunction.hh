/*
 * ax1_subdomainindexgridfunction.hh
 *
 *  Created on: Aug 29, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_SUBDOMAININDEXGRIDFUNCTION_HH
#define DUNE_AX1_SUBDOMAININDEXGRIDFUNCTION_HH

//! Simple as that: Return the subdomainIndex of an element (one of CYTOSOL, ES, MEMBRANE)
template<typename GV, typename RF, typename PHYSICS>
class Ax1SubdomainIndexGridFunction
  : public Dune::PDELab::GridFunctionBase<
          Dune::PDELab::GridFunctionTraits<GV,
                                           RF,
                                           1, // This should be enough for every channel (gbar, g and up to 3 gating particles)
                                           Dune::FieldVector<RF,1> >,
         Ax1SubdomainIndexGridFunction<GV,RF,PHYSICS> >
{
  public:
    typedef Dune::PDELab::GridFunctionTraits<GV,   // grid view type
                                           RF, // image field type (double)
                                           1,                                        // number of components of image (k)
                                           Dune::FieldVector<RF,1> // image type (Dune::FieldVector<double, k>)
                                           > Traits;

    //! constructor
    Ax1SubdomainIndexGridFunction (const GV& gv_, const PHYSICS& physics_)
      : gridView(gv_), physics(physics_)
    {
    }

    //! \copydoc GridFunctionBase::evaluate()
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                          typename Traits::RangeType& y) const
    {
      y = physics.getSubdomainIndex(e);
    }

    inline const typename Traits::GridViewType& getGridView () const
    {
      return gridView;
    }

  private:
    const typename Traits::GridViewType& gridView;
    const PHYSICS& physics;
};

#endif /* DUNE_AX1_SUBDOMAININDEXGRIDFUNCTION_HH */
