/*
 * ax1_multidomaingridfunction.hh
 *
 *  Created on: Mar 23, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_MULTIDOMAINGRIDFUNCTION_HH
#define DUNE_AX1_MULTIDOMAINGRIDFUNCTION_HH

#include <dune/ax1/common/constants.hh>

//! \brief
  /**
   * This class puzzles together the grid function on two subdomains to yield a grid function on the
   * multidomain grid
   *
   * The class allows the two grid functions defined on subdomains to have different types
   */
template<typename GV, typename PHYSICS, typename DGF1, typename DGF2 = DGF1>
class Ax1MultiDomainGridFunction
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,
               typename DGF1::Traits::RangeFieldType,
               DGF1::Traits::RangeType::dimension,
               typename DGF1::Traits::RangeType>,
               Ax1MultiDomainGridFunction<GV, DGF1, DGF2, PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,                           // grid view type
               typename DGF1::Traits::RangeFieldType,                     // range field type (double)
               DGF1::Traits::RangeType::dimension,                        // number of components of image (1)
               typename DGF1::Traits::RangeType                           // image type (Dune::FieldVector<double, 1>)
               > Traits;

  //typedef typename DGF::Traits Traits;

  typedef typename DGF1::Traits::RangeFieldType RF;
  typedef typename DGF1::Traits::RangeType RT;


  //! constructor
  Ax1MultiDomainGridFunction (const GV& gv_, DGF1& dgf1_, DGF2& dgf2_, PHYSICS& physics_)
    : gridView(gv_), dgf1(dgf1_), dgf2(dgf2_), physics(physics_)
  {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    //debug_jochen << "Ax1MultiDomainGridFunction::evaluate @ x = " << e.geometry().center() << std::endl;

    if(not physics.isMembrane(e))
    {
      typename PHYSICS::SubDomainElementPointer sdep = dgf1.getGridView().grid().subDomainEntityPointer(e);
      assert(dgf1.getGridView().contains(*sdep));
      dgf1.evaluate(*sdep,x,y);
    } else {
      typename PHYSICS::SubDomainElementPointer sdep = dgf2.getGridView().grid().subDomainEntityPointer(e);
      assert(dgf2.getGridView().contains(*sdep));
      dgf2.evaluate(*sdep,x,y);
    }
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  const GV& gridView;
  DGF1& dgf1;
  DGF2& dgf2;
  PHYSICS& physics;
};


//! \brief
/**
 * This class puzzles together the grid function on two subdomains to yield a grid function on the
 * multidomain grid
 *
 * The class allows the two grid functions defined on subdomains to have different types
 */
template<typename GV, typename PHYSICS, typename DGF1, typename DGF2 = DGF1>
class Ax1ElectrolytePlusMembraneIntersectionToMultiDomainGridFunction
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,
               typename DGF1::Traits::RangeFieldType,
               DGF1::Traits::RangeType::dimension,
               typename DGF1::Traits::RangeType>,
               Ax1ElectrolytePlusMembraneIntersectionToMultiDomainGridFunction<GV, DGF1, DGF2, PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,                           // grid view type
               typename DGF1::Traits::RangeFieldType,                     // range field type (double)
               DGF1::Traits::RangeType::dimension,                        // number of components of image (1)
               typename DGF1::Traits::RangeType                           // image type (Dune::FieldVector<double, 1>)
               > Traits;

  //typedef typename DGF::Traits Traits;

  typedef typename DGF1::Traits::RangeFieldType RF;
  typedef typename DGF1::Traits::RangeType RT;


  //! constructor
  Ax1ElectrolytePlusMembraneIntersectionToMultiDomainGridFunction (
      const GV& gv_, DGF1& dgf1_, DGF2& dgf2_, PHYSICS& physics_)
    : gridView(gv_), dgf1(dgf1_), dgf2(dgf2_), physics(physics_)
  {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    //debug_jochen << "Ax1MultiDomainGridFunction::evaluate @ x = " << e.geometry().center() << std::endl;

    if(not physics.isMembrane(e))
    {
      typename PHYSICS::SubDomainElementPointer sdep = dgf1.getGridView().grid().subDomainEntityPointer(e);
      assert(dgf1.getGridView().contains(*sdep));
      dgf1.evaluate(*sdep,x,y);
    } else {
      // For membrane elements: Serach for the matching membrane intersection for which dgf2 is defined; this
      // is the membrane intersection in radial direction at the membrane-cytosol interface, i.e. the one
      // below this element! Use a combination of Physics methods to find it and evaluate it for each of
      // possibly multiple) membrane element layers for the same x coordinate!

      const typename PHYSICS::ElementIntersectionIterator& iit = physics.getNextMembraneInterface(e);
      int iIndex = physics.getIntersectionIndex(*iit);

      //debug_jochen << "Found membrane intersection #" << iIndex << " @" << iit->geometry().center() << std::endl;

      // Something went wrong here, the found membrane intersection is not a primary one!
      if(! physics.isIntersectionInMembraneInterfacesVector(*iit))
        DUNE_THROW(Dune::Exception, "Could not find a matching membrane interface for Intersection GF!");

      typename DGF2::Traits::RangeType yIntersection(0.0);
      dgf2.evaluate(*iit,iit->geometryInInside().local(x),yIntersection);
      for(int i=0; i<y.size(); i++)
      {
        assert(Traits::GridViewType::dimension == 2);
        // This works only under the assumption dimGrid = 2 and membrane orientation is horizontal!
        if(i % 2 != 0)
        {
          y[i] = yIntersection[i / 2];
        } else {
          y[i] = 0.0;
        }
      }

      return;

    }
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  const GV& gridView;
  DGF1& dgf1;
  DGF2& dgf2;
  PHYSICS& physics;
};


#endif /* DUNE_AX1_MULTIDOMAINGRIDFUNCTION_HH */
