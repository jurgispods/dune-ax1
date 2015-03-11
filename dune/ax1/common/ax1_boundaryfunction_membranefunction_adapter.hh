/*
 * ax1_multidomaingridfunction.hh
 *
 *  Created on: Mar 23, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_BOUNDARYFUNCTION_MEMBRANEFUNCTION_ADAPTER_HH
#define DUNE_AX1_BOUNDARYFUNCTION_MEMBRANEFUNCTION_ADAPTER_HH

#include <dune/ax1/common/constants.hh>


/**
 * This class converts a given grid function living on the interfaces of a subdomain SubGV
 * to a grid function living on the elements of the subdomain SubGV
 */
template<typename SubGV, typename DGF, typename PHYSICS>
class Ax1BoundaryFunctionMembraneFunctionAdapter
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<SubGV,
               typename DGF::Traits::RangeFieldType,
               DGF::Traits::RangeType::dimension * SubGV::dimension,
               // # components = # components of membrane flux image (= number of species) times number of space dimensions
               //typename DGF::Traits::RangeType>,
               Dune::FieldVector<typename DGF::Traits::RangeFieldType,
                DGF::Traits::RangeType::dimension * SubGV::dimension
               > >,
     Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV, DGF, PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<SubGV,                           // grid view type
               typename DGF::Traits::RangeFieldType,                     // range field type (double)
               DGF::Traits::RangeType::dimension * SubGV::dimension, // number of components of image
               //typename DGF::Traits::RangeType                           // image type (Dune::FieldVector<double, 1>)
               //Dune::FieldVector<
               Dune::FieldVector<typename DGF::Traits::RangeFieldType,
                 DGF::Traits::RangeType::dimension * SubGV::dimension
               >
               //  DGF::Traits::RangeType::dimension> // Crude hack to comply with IonFluxGridFunction: Wrap this into a fieldvector of dimension 1
               > Traits;

  //typedef typename DGF::Traits Traits;

  typedef typename DGF::Traits::RangeFieldType RF;
  typedef typename DGF::Traits::RangeType RT;


  //! constructor
  Ax1BoundaryFunctionMembraneFunctionAdapter (const SubGV& gv_, DGF& dgf_, PHYSICS& physics_)
    : gridView(gv_), dgf(dgf_), physics(physics_)
  {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    assert(dgf.getGridView().contains(gridView.grid().multiDomainEntity(e)));

    const int dimDomain = SubGV::dimension;

    //debug_jochen << "Visiting element @ " << e.geometry().center() << std::endl;

    //const typename PHYSICS::Element& me = gridView.grid().multiDomainEntity(e);
    //for(typename PHYSICS::ElementIntersectionIterator miit = me.ileafbegin(); miit != me.ileafend(); ++miit)
    for(typename PHYSICS::SubDomainElementIntersectionIterator siit = gridView.ibegin(e);
        siit != gridView.iend(e); ++siit)
    {
      if(physics.isMembraneInterface(*siit))
      {
        //debug_jochen << "Visiting intersection @ " << siit->geometry().center() << std::endl;

        Dune::PDELab::IntersectionGeometry<typename PHYSICS::SubDomainElementIntersection> ig(*siit, -1);

        // FieldVector of size 3
        typename DGF::Traits::RangeType membFlux;
        // Assume the value at this intersection is valid for the whole element, therefore x doesn't matter
        dgf.evaluate(ig, ig.geometry().local(ig.geometry().center()),membFlux);

        // We evaluate the flux from a membrane element, therefore signs are flipped!
        membFlux *= -1;

        //debug_jochen << "* " << membFlux << std::endl;
        for(int i=0; i<y.size(); i++)
        {
          int j = i / dimDomain;

          // Implicitly convert double (membFlux[j]) to Dune::FieldVector<double,1> (y[j])
          // TODO Make this work for general grids
          if(dimDomain == 2)
          {
            if(i % dimDomain != 0)
            {
              y[i] = membFlux[j]; // Set y-component of this element's ion flux to the membrane flux
            }
          } else { // dim = 1
            y[i] = membFlux[j];
          }
          //debug_jochen << "y[" << i << "] = " << y[i] << std::endl;

        }
        return;
      }
    }

    DUNE_THROW(Dune::Exception,
        "No membrane intersection found for this element. ""Cannot evaluate underlying boundary gridfunction!");
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  const SubGV& gridView;
  DGF& dgf;
  PHYSICS& physics;
};


#endif /* DUNE_AX1_BOUNDARYFUNCTION_MEMBRANEFUNCTION_ADAPTER_HH */
