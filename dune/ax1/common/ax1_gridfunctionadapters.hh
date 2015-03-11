/*
 * ax1_extrapolationgridfunction.hh
 *
 *  Created on: Aug 27, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_COARSEGRIDINTERPOLATIONGRIDFUNCTION
#define DUNE_AX1_COARSEGRIDINTERPOLATIONGRIDFUNCTION

/**!
 * \note NEW VERSION of the grid data transfer functions defined below
 *
 * \brief This gridfunction acts as if it was defined on a gridview FineGV, but actually only interpolates the
 * values from a coarser gridview CoarseGV. Useful when transferring the system equilibrium state computed
 * on a coarse grid to a much larger (number-of-DOFs-wise) grid
 *
 * This method assumes that coarse and fine grid elements are order in the same way: First elements in x-direction,
 * then in y-direction!
 */
template<typename FineGV, typename CoarseGF, typename PHYSICS>
class Ax1CoarseToFineGridTransferGridFunction
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<FineGV,
                                       typename CoarseGF::Traits::RangeFieldType,
                                       CoarseGF::Traits::RangeType::dimension,
                                       typename CoarseGF::Traits::RangeType>,
                                       Ax1CoarseToFineGridTransferGridFunction<FineGV,CoarseGF,PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<FineGV,   // grid view type
                                            typename CoarseGF::Traits::RangeFieldType, // image field type (double)
                                            CoarseGF::Traits::RangeType::dimension,                                        // number of components of image (k)
                                            typename CoarseGF::Traits::RangeType // image type (Dune::FieldVector<double, k>)
                                           > Traits;

  typedef typename CoarseGF::Traits::GridViewType CoarseGV;
  typedef typename CoarseGV::template Codim<0>::Entity CoarseEntity;
  typedef typename CoarseGV::template Codim<0>::EntityPointer CoarseEntityPointer;
  typedef typename CoarseGV::template Codim<0>::
    template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator CoarseElementLeafIterator;

  //! constructor
  Ax1CoarseToFineGridTransferGridFunction (const FineGV& fineGV_, const std::vector<Dune::shared_ptr<CoarseGF> >& coarseGF_,
      const PHYSICS& physics_, int nElementsXCoarse_, int nElementsYCoarse_,
      std::vector<int>& coarseGroups_, bool doInterpolate_ = false, bool useGridConvergenceMode_ = false)
    : fineGV(fineGV_),
      coarseGF(coarseGF_),
      coarseGV(coarseGF_[0]->getGridView()),
      physics(physics_),
      nElementsXCoarse(nElementsXCoarse_),
      nElementsYCoarse(nElementsYCoarse_),
      coarseGroups(coarseGroups_),
      nElementsXFine(physics_.nElements(0)),
      nElementsYFine(physics_.nElements(1)),
      nElementsFine(0),
      lastElementIndex(std::numeric_limits<int>::lowest()),
      countEvaluated(0),
      isThisGridFunctionCompatible(true),
      // There is not overlap anyway in the coarse grid, as nElementsXCoarse is assumed to be 1!
      ceit(coarseGV.template begin<0,Dune::PartitionIteratorType::Interior_Partition>()),
      // Initialize pointer to current coarse entity
      doInterpolate(doInterpolate_),
      useGridConvergenceMode(useGridConvergenceMode_),
      verbose(false)
  {
    if(!useGridConvergenceMode)
    {
      assert(nElementsXCoarse ==  1);

      // Loop over first row of elements and count, including overlap entities!
      typedef typename FineGV::template Codim<0>::Iterator FineElementLeafIterator;
      FineElementLeafIterator feit = fineGV.template begin<0>();
      typename Traits::DomainType x_start = feit->geometry().center();
      for(; feit != fineGV.template end<0>(); ++feit)
      {
        if(feit->geometry().center()[1] > x_start[1])
          break;

        if(feit->partitionType() != Dune::PartitionType::InteriorEntity)
        {
          nElementsXFine++;
          //debug_jochen << feit->geometry().center() << " -> nElementsXFine = " << nElementsXFine;
        }
      }

      assert(nElementsXFine <= physics.nElements(0) + 2*fineGV.overlapSize(0));

      debug_verb << "nElementsXCoarse = " << nElementsXCoarse << std::endl;
      debug_verb << "nElementsYCoarse = " << nElementsYCoarse << std::endl;
      debug_verb << "nElementsXFine = " << nElementsXFine << std::endl;
      debug_verb << "nElementsYFine = " << nElementsYFine << std::endl;

      if(nElementsYFine % nElementsYCoarse != 0)
      {
        nElementsYFine -= physics.getParams().nMembraneElements();
        nElementsYCoarse -= physics.getParams().nMembraneElements();

        if(nElementsYFine % nElementsYCoarse != 0)
        {
          debug_warn << "This gridfunction is not compatible for transferring values from fine to coarse grid!"
              << std::endl;

          isThisGridFunctionCompatible = false;
          nElementsYFine += physics.getParams().nMembraneElements();
          nElementsYCoarse += physics.getParams().nMembraneElements();
        }
      }

      nElementsFine = (nElementsXFine / nElementsXCoarse) * (nElementsYFine / nElementsYCoarse);
      debug_verb << "Using " << nElementsFine << " fine elements for restriction!" << std::endl;
    } else {
      nElementsFine = std::numeric_limits<int>::max();
    }
  }

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    if(! isThisGridFunctionCompatible)
    {
      DUNE_THROW(Dune::Exception, "This GF is not compatible for transferring value from coarse to fine grid:"
          << "nElementsXCoarse = " << nElementsXCoarse
          << ", nElementsYCoarse = " << nElementsYCoarse
          << ", nElementsXFine = " << nElementsXFine
          << ", nElementsYFine = " << nElementsYFine);
    }

    assert(ceit != (coarseGV.template end<0,Dune::PartitionIteratorType::Interior_Partition>()));

    y = 0.0;

    // New element? Count it!
    if(physics.getElementIndex(e) != lastElementIndex)
    {
      countEvaluated++;
      lastElementIndex = physics.getElementIndex(e);
    }

    // Check if searched y coordinate is actually contained in the current coarse element ceit
    typename Traits::DomainType xglobal = e.geometry().global(x);
    typename Traits::DomainType xlocal_coarse = ceit->geometry().local(xglobal);

    bool containsYCoordinate = (xlocal_coarse[1]-1e-6 <= 1 && xlocal_coarse[1]+1e-6 >= 0);

    // After we evaluated the coarse entity nElementsXFine times, we need to go to the next coarse entity!
    // This might also be the case if the desired y coordinate is not contained in the coarse entity
    if(countEvaluated > nElementsFine || !containsYCoordinate)
    {
      //debug_jochen << "Finished interpolating from coarse entity @" << ceit->geometry().center()
      //    << ", used it a total of " << (countEvaluated -1) << " times!" << std::endl;

      countEvaluated = 1;
      // Move coarse entity pointer to the next row (next y position);
      ++ceit;
    }

    // I don't get it. Why did this uncommented older version only use the x coordinate, but not y?
//    // This will give a value between 0 and 1 and can be used as a local coordinate in the coarse entity
//    typename Traits::DomainFieldType xnew = (e.geometry().global(x)[0] - physics.getParams().xMin()) /
//        (physics.getParams().xMax() - physics.getParams().xMin());
//
//    typename Traits::DomainType xcoarse(x);
//    xcoarse[0] = xnew;


    // xnew has now local coordinates with respect to ce!
    typename Traits::DomainType xnew = ceit->geometry().local(e.geometry().global(x));

    typename Traits::DomainType diffCheck = e.geometry().global(x);
    diffCheck -= ceit->geometry().global(xnew);
    //debug_jochen << "coarse/fine grid position: " << ce.geometry().global(xnew) << " / "
    //    << e.geometry().global(x) << std::endl;
    assert(diffCheck.two_norm() < 1e-6);

    assert(coarseGF.size() <= 2);

    // Do interpolation here, coarseGF is now a vector of gridfunctions!
    // Standard evaluation, no interpolation
    if(coarseGF.size() == 1)
    {
      // Evaluate the coarse gridfunction
      coarseGF[0]->evaluate(*ceit,xnew,y);
    } else {
      // When on membrane: interpolate values!

      // Hardcoded version for maximum 2 membrane groups
      double perm1 = physics.getGroupPermittivity(coarseGroups[0]);
      double perm2 = physics.getGroupPermittivity(coarseGroups[1]);

      double d_perm = perm2-perm1;

      // How the interpolation should generally work (at the moment hardcoded for 2 gridfunctions):
      // Check x-coordinate to determine between which two membrane groups this element lies; then it is also easy to
      // pick exactly those two gridfunctions that we want to interpolate here, that means we can even allow more than 2
      // coarse GF in the vector!
      // - Get group index of membrane intervals left and right
      // - Get (possibly smoothed) membrane permittivity for this x-coordinate, regardless if we are on the membrane or not
      // - Interpolate the two picked GFs (using the left/right group indices) using their membrane group permittivity and
      //   my membrane group permittivity

      // Find the membrane element with the same x-coordinate than the current element.
      // We need the (membrane) permittivity at that very element!
      double correspondingMembranePermittivity = -1;
      for(typename PHYSICS::MIterator mit = physics.mBegin(); mit != physics.mEnd(); ++mit)
      {
        // Compare x coordinates
        if(std::abs(mit->geometry().center()[0] - e.geometry().center()[0]) < 1e-6)
        {
          correspondingMembranePermittivity = physics.getPermittivity(*mit->inside());
          break;
        }
      }
      assert(correspondingMembranePermittivity > -1);

      double scale = (correspondingMembranePermittivity - perm1) / d_perm;

      if(verbose && physics.isMembrane(e))
      {
        debug_jochen << "Element @ " << e.geometry().center() << std::endl;
        debug_jochen << "My membrane group permittivity = " << correspondingMembranePermittivity << " => scale = " << scale << std::endl;
      }

      if(! doInterpolate)
      {
        scale = std::round(scale);
      }

      assert(scale >= 0.0 && scale <= 1.0);

      // Do linear interpolation on y
      typename Traits::RangeType y_interp1(0.0);
      coarseGF[0]->evaluate(*ceit,xnew,y_interp1);

      y = y_interp1;

      typename Traits::RangeType y_interp2(0.0);
      coarseGF[1]->evaluate(*ceit,xnew,y_interp2);
      y_interp2 -= y_interp1;
      y_interp2 *= scale;

      y += y_interp2;
    }

    if(verbose && physics.isMembrane(e))
      debug_jochen << "Fine element #" << physics.getElementIndex(e) << "[" << countEvaluated << "], "
        << "transfer value " << y << " @ coarse: " << ceit->geometry().global(xnew)
        << " --> fine: " << e.geometry().global(x) << std::endl;

  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return fineGV;
  }

  bool isCompatible() const
  {
    return isThisGridFunctionCompatible;
  }

  void setVerbose(bool verbose_)
  {
    verbose = verbose_;
  }

private:

  const FineGV& fineGV;
  const std::vector<Dune::shared_ptr<CoarseGF> >& coarseGF;
  const CoarseGV& coarseGV;
  const PHYSICS& physics;

  int nElementsXCoarse;
  int nElementsYCoarse;
  int nElementsXFine;
  int nElementsYFine;
  int nElementsFine;

  mutable int lastElementIndex;
  mutable int countEvaluated;
  bool isThisGridFunctionCompatible;
  mutable CoarseElementLeafIterator ceit;

  std::vector<int> coarseGroups;
  bool doInterpolate;
  bool useGridConvergenceMode;
  bool verbose;
};


/**!
 * \note NEW VERSION of the grid data transfer functions defined below
 *
 * \brief This gridfunction acts as if it was defined on a gridview CoarseGV, but actually only interpolates the
 * values from a finer gridview CoarseGV. Useful when transferring the element information from a finer grid
 * to a coarse grid (possibly loaded from a previous simulation, but not tagged with element subdomains)
 *
 * This method assumes that coarse and fine grid elements are order in the same way: First elements in x-direction,
 * then in y-direction!
 */
template<typename CoarseGV, typename FineGF, typename PHYSICS>
class Ax1FineToCoarseGridTransferGridFunction
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<CoarseGV,
                                       typename FineGF::Traits::RangeFieldType,
                                       FineGF::Traits::RangeType::dimension,
                                       typename FineGF::Traits::RangeType>,
                                       Ax1FineToCoarseGridTransferGridFunction<CoarseGV,FineGF,PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<CoarseGV,   // grid view type
                                            typename FineGF::Traits::RangeFieldType, // image field type (double)
                                            FineGF::Traits::RangeType::dimension,                                        // number of components of image (k)
                                            typename FineGF::Traits::RangeType // image type (Dune::FieldVector<double, k>)
                                           > Traits;

  typedef typename FineGF::Traits::GridViewType FineGV;
  typedef typename FineGV::template Codim<0>::Entity FineEntity;
  typedef typename FineGV::template Codim<0>::EntityPointer FineEntityPointer;
  typedef typename FineGV::template Codim<0>::
    template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator FineElementLeafIterator;

  //! constructor
  Ax1FineToCoarseGridTransferGridFunction (const CoarseGV& coarseGV_, const FineGF& fineGF_,
      const PHYSICS& physics_, int nElementsXCoarse_, int nElementsYCoarse_)
    : coarseGV(coarseGV_),
      fineGF(fineGF_),
      fineGV(fineGF_.getGridView()),
      physics(physics_),
      nElementsXCoarse(nElementsXCoarse_),
      nElementsYCoarse(nElementsYCoarse_),
      nElementsXFine(physics_.nElements(0)),
      nElementsYFine(physics_.nElements(1)),
      nElementsFine(0),
      lastElementIndex(std::numeric_limits<int>::lowest()),
      countEvaluated(0),
      isThisGridFunctionCompatible(true),
      // It is more or less arbitrary which partition to choose, as we never really get to
      // see the whole grid anyway on this processor!
      feit(fineGV.template begin<0,Dune::PartitionIteratorType::Interior_Partition>())
      // Initialize pointer to current coarse entity
  {
    assert(nElementsXCoarse == 1);

    //nElementsYFine = fineGV.size(0) / nElementsXFine;

//    if(nElementsXFine % nElementsXCoarse != 0
//        || nElementsYFine % nElementsYCoarse != 0)
//    {
//      DUNE_THROW(Dune::Exception,  "Number of elements of fine and coarse grid are not compatible!"
//          << " nElementsXFine = " << nElementsXFine << ", nElementsXCoarse = " << nElementsXCoarse
//          << " nElementsYFine = " << nElementsYFine << ", nElementsYCoarse = " << nElementsYCoarse);
//    }
//
//    nElementsXFine /= nElementsXCoarse;
//    nElementsYFine /= nElementsYCoarse;

    debug_verb << "nElementsXCoarse = " << nElementsXCoarse << std::endl;
    debug_verb << "nElementsYCoarse = " << nElementsYCoarse << std::endl;
    debug_verb << "nElementsXFine = " << nElementsXFine << std::endl;
    debug_verb << "nElementsYFine = " << nElementsYFine << std::endl;

    if(nElementsYFine % nElementsYCoarse != 0)
    {
      nElementsYFine -= physics.getParams().nMembraneElements();
      nElementsYCoarse -= physics.getParams().nMembraneElements();

      if(nElementsYFine % nElementsYCoarse != 0)
      {
        debug_warn << "This gridfunction is not compatible for transferring values from fine to coarse grid!"
            << std::endl;

        isThisGridFunctionCompatible = false;
        nElementsYFine += physics.getParams().nMembraneElements();
        nElementsYCoarse += physics.getParams().nMembraneElements();
      }
    }

    nElementsFine = (nElementsXFine / nElementsXCoarse) * (nElementsYFine / nElementsYCoarse);
    debug_verb << "Using " << nElementsFine << " fine elements for restriction!" << std::endl;
  }

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    if(! isThisGridFunctionCompatible)
    {
      DUNE_THROW(Dune::Exception, "This GF is not compatible for transferring value from fine to coarse grid:"
          << " nElementsXCoarse = " << nElementsXCoarse
          << ", nElementsYCoarse = " << nElementsYCoarse
          << ", nElementsXFine = " << nElementsXFine
          << ", nElementsYFine = " << nElementsYFine);
    }

    assert(feit != (fineGV.template end<0,Dune::PartitionIteratorType::Interior_Partition>()));

    y = 0.0;

    typename Traits::DomainType xglobal = e.geometry().global(x);
    //debug_jochen << "Coarse element @" << xglobal << std::endl;

    int countEvaluated = 0;
    // In order to interpolate the fine grid values onto the coarse grid, we have to average the
    // nElementsFine values from the fine grid corresponding to the nElementsXCoarse from the
    // coarse grid (currently, always 1); note tha nElementsFine is an upper bound for the number
    // of fine elements; in regions, which contain less fine elements (e.g. membrane which is only
    // refined in x-direction), we have to be careful. The following loop tackles this by additionally
    // checking the y-coordinate and only considers fine elements which containt this coordinate. I.e.,
    // this procedure is not an actual averaging of all fine elements contained in the coarse elements,
    // but only
    for(int i=0; i<nElementsFine; i++)
    {
      // Check if this fine element is actually contained in the coarse elements
      typename Traits::DomainType xlocal_coarse = e.geometry().local(feit->geometry().center());

      bool containsYCoordinate = (xlocal_coarse[1]-1e-6 <= 1 && xlocal_coarse[1]+1e-6 >= 0);

      //debug_jochen << "Fine element @" << feit->geometry().center() << ", xlocal_coarse = "
      //   << xlocal_coarse << std::endl;

      if(containsYCoordinate)
      {
        countEvaluated++;

        //debug_jochen << "Using fine element @" << feit->geometry().center()
        //   << std::endl;

        typename Traits::RangeType y_fine;

        // The coordinate x does not really matter here! We just evaluate at the center
        typename Traits::DomainType xlocal_fine(0.5);
        fineGF.evaluate(*feit,xlocal_fine,y_fine);
        y += y_fine;
        ++feit;
      } else {
        // If the current fine entity (y position) is not contained the coarse entity, distinguish two cases:
        if(xlocal_coarse[1] < 0)
        {
          // (1) The fine entity lies below the desired coarse entity; increment!
          ++feit;
        } else {
          // (2) The fine entity lies above the desired coarse entity; break loop and return!
          break;
        }
      }
    }


    // Average value from all fine elements
    y /= countEvaluated;

    //debug_jochen << "Coarse element @" << e.geometry().global(x) << ", used " << countEvaluated
    //   << "fine elements --> transfer value " << y << std::endl;

  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return coarseGV;
  }

  bool isCompatible() const
  {
    return isThisGridFunctionCompatible;
  }

private:

  const CoarseGV& coarseGV;
  const FineGF& fineGF;
  const FineGV& fineGV;
  const PHYSICS& physics;

  int nElementsXCoarse;
  int nElementsYCoarse;
  int nElementsXFine;
  int nElementsYFine;
  int nElementsFine;

  mutable int lastElementIndex;
  mutable int countEvaluated;
  bool isThisGridFunctionCompatible;
  mutable FineElementLeafIterator feit;
};


/**!
 * This gridfunction acts as if it was defined on a gridview FineGV, but actually only interpolates the
 * values from a coarser gridview CoarseGV. Useful when transferring the system equilibrium state computed
 * on a coarse grid to a much larger (number-of-DOFs-wise) grid
 *
 * This is a terribly unperformant implementation, as it is iterated over the CoarseGV for every function call.
 * Therefore it should not be used extensively; best only once at the beginning of the simulation!
 */
template<typename FineGV, typename CoarseGF, typename PHYSICS>
class Ax1CoarseGridInterpolationGridFunction
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<FineGV,
                                       typename CoarseGF::Traits::RangeFieldType,
                                       CoarseGF::Traits::RangeType::dimension,
                                       typename CoarseGF::Traits::RangeType>,
                                       Ax1CoarseGridInterpolationGridFunction<FineGV,CoarseGF,PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<FineGV,   // grid view type
                                            typename CoarseGF::Traits::RangeFieldType, // image field type (double)
                                            CoarseGF::Traits::RangeType::dimension,                                        // number of components of image (k)
                                            typename CoarseGF::Traits::RangeType // image type (Dune::FieldVector<double, k>)
                                           > Traits;

  typedef typename CoarseGF::Traits::GridViewType CoarseGV;
  typedef typename CoarseGV::template Codim<0>::Entity CoarseEntity;
  typedef typename CoarseGV::template Codim<0>::EntityPointer CoarseEntityPointer;
  typedef typename CoarseGV::template Codim<0>::Iterator CoarseElementLeafIterator;


  //! constructor
  Ax1CoarseGridInterpolationGridFunction (const FineGV& fineGV_,
      const CoarseGF& coarseGF_, const PHYSICS& physics_)
    : fineGV(fineGV_),
      coarseGF(coarseGF_),
      coarseGV(coarseGF_.getGridView()),
      physics(physics_),
      zero(0.0),
      one(1.0)
  {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    // Loop over coarseGV's elements and find that one that contains e
    for(CoarseElementLeafIterator ceit = coarseGV.template begin<0>();
        ceit != coarseGV.template end<0>(); ++ceit)
    {
      CoarseEntity& ce = *ceit;

      // Check if the coarse entity contains the position x
      if(contains(ce, e.geometry().global(x)))
      {
        // Convert the fine-entity local coordinate to a coarse-entity local coordinate!

        //typename Traits::DomainType e_diag = e.geometry().corner(e.geometry().corners()-1);
        //e_diag -= e.geometry().corner(0);

        //typename Traits::DomainType ce_diag = ce.geometry().corner(ce.geometry().corners()-1);
        //ce_diag -= ce.geometry().corner(0);

        // xnew has now local coordinates with respect to ce!
        typename Traits::DomainType xnew = ce.geometry().local(e.geometry().global(x));

        typename Traits::DomainType diffCheck = e.geometry().global(x);
        diffCheck -= ce.geometry().global(xnew);
        //debug_jochen << "coarse/fine grid position: " << ce.geometry().global(xnew) << " / "
        //    << e.geometry().global(x) << std::endl;
        assert(diffCheck.two_norm() < 1e-6);

        // Evalute the coarse gridfunction with respect to the new local coordinate
        coarseGF.evaluate(ce,xnew,y);

//        if(std::abs(ceit->geometry().global(x)[1]-502.5) < 5)
//          debug_jochen << "Fine element #" << physics.getElementIndex(e) << ", "
//            << "transfer value " << y << " at " << ceit->geometry().global(x)
//            << " --> " << e.geometry().global(x) << std::endl;

        // Coarse grid element containing e was found => return!
        return;
      }
    }

    DUNE_THROW(Dune::Exception, "Could not find a matching coarse entity for the given entity!"
        << " Fine grid position: " << e.geometry().global(x));
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return fineGV;
  }

private:

  // Quick and dirty containment check
  inline bool contains(const CoarseEntity& ce, const typename Traits::DomainType& xglobal) const
  {
    //typename Traits::DomainType e_lowerLeftCorner = e.geometry().corner(0);
    //typename Traits::DomainType e_upperRightCorner = e.geometry().corner(e.geometry().corners()-1);

    //typename Traits::DomainType ce_lowerLeftCorner = ce.geometry().corner(0);
    //typename Traits::DomainType ce_upperRightCorner = ce.geometry().corner(e.geometry().corners()-1);


    //debug_verb << "Checking coarse entity [" << ce_lowerLeftCorner << ", " << ce_upperRightCorner << "]  <->  ["
    //        << e_lowerLeftCorner << ", " << e_upperRightCorner << "] : " << success << std::endl;

    //return (Tools::lessOrEqualThan(ce_lowerLeftCorner, e_lowerLeftCorner, 1e-6)
    //         &&  Tools::greaterOrEqualThan(ce_upperRightCorner, e_upperRightCorner, 1e-6));



    typename Traits::DomainType ceLocal = ce.geometry().local(xglobal);

    return (Tools::lessOrEqualThan(ceLocal, one, 1e-6)
         &&  Tools::greaterOrEqualThan(ceLocal, zero, 1e-6));
  }

  const FineGV& fineGV;
  const CoarseGF& coarseGF;
  const CoarseGV& coarseGV;
  const PHYSICS& physics;

  const typename Traits::DomainType zero;
  const typename Traits::DomainType one;
};


/**!
 * This gridfunction acts as if it was defined on a gridview CoarseGV, but actually only restricts the
 * values from a finer gridview FineGV. Useful when transferring multidomain/subdomain information
 * from a fine grid to a coarse grid (number-of-DOFs-wise)
 *
 * The restriction is done by taking the value from that fine grid entity that matches the given coordinate.
 * If it matches multiple entities (i.e., the coordinate lies on a fine grid intersection), the arithmetic
 * average is taken.
 *
 * This is a terribly unperformant implementation, as it is iterated over the FineGV for every function call.
 * Therefore it should not be used extensively; best only once at the beginning of the simulation!
 *
 * The constructor parameter maxFineElements can be used to speed things up a bit: The user can specify
 * a maximum number of fine elements that may contain a position of the coarse grid. In structured Cartesian
 * 2D grids, this would for example 4 at the nodes, 2 at edges. The default value is 8. (Is there a case where
 * there ever might be more than 8? I don't think so).
 */
template<typename CoarseGV, typename FineGF, typename PHYSICS>
class Ax1FineGridRestrictionGridFunction
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<CoarseGV,
                                       typename FineGF::Traits::RangeFieldType,
                                       FineGF::Traits::RangeType::dimension,
                                       typename FineGF::Traits::RangeType>,
    Ax1FineGridRestrictionGridFunction<CoarseGV,FineGF,PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<CoarseGV,   // grid view type
                                            typename FineGF::Traits::RangeFieldType, // image field type (double)
                                            FineGF::Traits::RangeType::dimension,                                        // number of components of image (k)
                                            typename FineGF::Traits::RangeType // image type (Dune::FieldVector<double, k>)
                                           > Traits;

  typedef typename FineGF::Traits::GridViewType FineGV;
  typedef typename FineGV::template Codim<0>::Entity FineEntity;
  typedef typename FineGV::template Codim<0>::EntityPointer FineEntityPointer;
  typedef typename FineGV::template Codim<0>::Iterator FineElementLeafIterator;


  //! constructor
  Ax1FineGridRestrictionGridFunction (const CoarseGV& coarseGV_,
      const FineGF& fineGF_, const PHYSICS& physics_, int maxFineElements_ = 8)
    : coarseGV(coarseGV_),
      fineGF(fineGF_),
      fineGV(fineGF_.getGridView()),
      physics(physics_),
      maxFineElements(maxFineElements_),
      zero(0.0),
      one(1.0)
  {
    // TODO Check if the both gridview are consistent (i.e. coarseGV subset fineGV)
  }

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    int count = 0;

    // Loop over fineGV's elements and find that one(s) that contains e
    for(FineElementLeafIterator feit = fineGV.template begin<0>();
        feit != fineGV.template end<0>(); ++feit)
    {
      FineEntity& fe = *feit;

      //debug_jochen << "Trying fine grid element @ " << feit->geometry().center() << std::endl;

      // Check if the fine entity fe contains the position x at which e shall be evaluated
      if(contains(fe, e.geometry().global(x)))
      {
        // Convert the fine-entity local coordinate to a coarse-entity local coordinate!
        typename Traits::DomainType xnew = fe.geometry().local(e.geometry().global(x));

        count++;

        // xnew has now local coordinates with respect to ce!
        typename Traits::DomainType diffCheck = e.geometry().global(x);
        diffCheck -= fe.geometry().global(xnew);
        //debug_jochen << "coarse/fine grid position: " << e.geometry().global(x) << " / "
        //    << fe.geometry().global(xnew) << std::endl;
        assert(diffCheck.two_norm() < 1e-6);

        // Evaluate the coarse gridfunction with respect to the new local coordinate
        typename FineGF::Traits::RangeType ynew(0.0);
        fineGF.evaluate(fe,xnew,ynew);
        y += ynew;

        // The maximum number of fine elements that may contain the searched position was
        // found, break search loop!
        if(count >= maxFineElements) break;
      }
    }

    // No fine grid entity was found!
    if(count == 0)
    {
      DUNE_THROW(Dune::Exception, "Could not find a matching fine entity for the given entity!"
          << " Coarse grid position: " << e.geometry().global(x));
    } else {
      // Finalize averaging by dividing by the number of evaluated fine grid entities
      y /= count;
    }

    //debug_jochen << "count = " << count << std::endl;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return fineGV;
  }

  //! Because I can.
  inline void makeMeASandwich () const
  {
    std::cout << "Here you are, Sir!" << std::endl;
  }

private:

  // Quick and dirty containment check
  inline bool contains(const FineEntity& fe, const typename Traits::DomainType& xglobal) const
  {
    //typename Traits::DomainType e_lowerLeftCorner = e.geometry().corner(0);
    //typename Traits::DomainType e_upperRightCorner = e.geometry().corner(e.geometry().corners()-1);

    //typename Traits::DomainType ce_lowerLeftCorner = ce.geometry().corner(0);
    //typename Traits::DomainType ce_upperRightCorner = ce.geometry().corner(e.geometry().corners()-1);


    //debug_verb << "Checking coarse entity [" << ce_lowerLeftCorner << ", " << ce_upperRightCorner << "]  <->  ["
    //        << e_lowerLeftCorner << ", " << e_upperRightCorner << "] : " << success << std::endl;

    //return (Tools::lessOrEqualThan(ce_lowerLeftCorner, e_lowerLeftCorner, 1e-6)
    //         &&  Tools::greaterOrEqualThan(ce_upperRightCorner, e_upperRightCorner, 1e-6));

    typename Traits::DomainType feLocal = fe.geometry().local(xglobal);

    return (Tools::lessOrEqualThan(feLocal, one, 1e-6)
         &&  Tools::greaterOrEqualThan(feLocal, zero, 1e-6));
  }

  const CoarseGV& fineGV;
  const FineGF& fineGF;
  const FineGV& coarseGV;
  const PHYSICS& physics;
  const int maxFineElements;

  const typename Traits::DomainType zero;
  const typename Traits::DomainType one;
};


template<typename GV,typename PHYSICS>
class Ax1RefinementMapGridFunction
  :  public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV, int, 1,
                                                                            Dune::FieldVector<int,1> >,
                                                                            Ax1RefinementMapGridFunction<GV,PHYSICS> >
{
  typedef Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV, int, 1,
                                         Dune::FieldVector<int,1> >,
                                         Ax1RefinementMapGridFunction<GV,PHYSICS> > BaseT;

public:
    typedef typename BaseT::Traits Traits;

    Ax1RefinementMapGridFunction(const GV& gv_, const PHYSICS& physics_,
        int nElementsXCoarse_, int nElementsYCoarse_)
    : gv(gv_),
      physics(physics_),
      nElementsXCoarse(nElementsXCoarse_),
      nElementsYCoarse(nElementsYCoarse_),
      nElementsXFine(physics.nElements(gv,0)),
      nElementsYFine(physics.nElements(gv,1)),
      nMembEl(physics.getParams().nMembraneElements()),
      factorX(nElementsXFine / nElementsXCoarse),
      factorY((nElementsYFine-nMembEl) / (nElementsYCoarse-nMembEl)),
      membraneVisited(false)
    {
      debug_jochen << "nElementsXFine = " << nElementsXFine << ", nElementsXCoarse = " << nElementsXCoarse
          << " => factorX = " << factorX << std::endl;
      debug_jochen << "nElementsYFine = " << nElementsYFine << ", nElementsYCoarse = " << nElementsYCoarse
          << " => factorY = " << factorY << std::endl;
    }


    // Evalute grid function; return 0 if the element is part of the membrane!
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                                typename Traits::RangeType& y) const
    {
      y = -1;
      int elemIndex = physics.getElementIndex(e);

      int col = (elemIndex % nElementsXFine) / factorX;
      int row;
      if(physics.isMembrane(e))
      {
        membraneVisited = true;
        row = elemIndex / (nElementsXFine*factorY);
      } else {
        int membOffset = (membraneVisited) ? (nElementsXFine * (factorY - nMembEl)) : 0;

        row = (membOffset+elemIndex) / (nElementsXFine*factorY);
      }
      //debug_jochen << " col=" << col << ", row=" << row << std::endl;
      y = row * nElementsXCoarse + col;
    }

    const typename Traits::GridViewType& getGridView() const
    {
      return gv;
    }

private:
    const GV& gv;
    const PHYSICS& physics;
    int nElementsXCoarse;
    int nElementsYCoarse;
    int nElementsXFine;
    int nElementsYFine;
    int nMembEl;
    int factorX;
    int factorY;
    mutable bool membraneVisited;
};

// Dirty hack to make the above Ax1RefinementMapGridFunction also work on a subdomain GV;
// this should be merged again by using template magic, but no time for this right now
template<typename SubGV,typename PHYSICS>
class Ax1RefinementMapGridFunction_Subdomain
  :  public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<SubGV, int, 1,
                                                                            Dune::FieldVector<int,1> >,
                                                                            Ax1RefinementMapGridFunction<SubGV,PHYSICS> >
{
  typedef Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<SubGV, int, 1,
                                         Dune::FieldVector<int,1> >,
                                         Ax1RefinementMapGridFunction<SubGV,PHYSICS> > BaseT;

public:
    typedef typename BaseT::Traits Traits;

    Ax1RefinementMapGridFunction_Subdomain(const SubGV& gv_, const PHYSICS& physics_,
        int nElementsXCoarse_, int nElementsYCoarse_)
    : gv(gv_),
      physics(physics_),
      nElementsXCoarse(nElementsXCoarse_),
      nElementsYCoarse(nElementsYCoarse_),
      nElementsXFine(physics.nElements(gv,0)),
      nElementsYFine(physics.nElements(gv,1)),
      nMembEl(0),
      factorX(nElementsXFine / nElementsXCoarse),
      factorY((nElementsYFine-nMembEl) / (nElementsYCoarse-nMembEl)),
      membraneVisited(false)
    {
      debug_jochen << "nElementsXFine = " << nElementsXFine << ", nElementsXCoarse = " << nElementsXCoarse
          << " => factorX = " << factorX << std::endl;
      debug_jochen << "nElementsYFine = " << nElementsYFine << ", nElementsYCoarse = " << nElementsYCoarse
          << " => factorY = " << factorY << std::endl;
    }


    // Evalute grid function; return 0 if the element is part of the membrane!
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                                typename Traits::RangeType& y) const
    {
      y = -1;
      int elemIndex = physics.getSubDomainElementIndex(e);

      int col = (elemIndex % nElementsXFine) / factorX;
      int row;
      if(physics.isMembrane(e))
      {
        membraneVisited = true;
        row = elemIndex / (nElementsXFine*factorY);
      } else {
        int membOffset = (membraneVisited) ? (nElementsXFine * (factorY - nMembEl)) : 0;

        row = (membOffset+elemIndex) / (nElementsXFine*factorY);
      }
      //debug_jochen << " col=" << col << ", row=" << row << std::endl;
      y = row * nElementsXCoarse + col;
    }

    const typename Traits::GridViewType& getGridView() const
    {
      return gv;
    }

private:
    const SubGV& gv;
    const PHYSICS& physics;
    int nElementsXCoarse;
    int nElementsYCoarse;
    int nElementsXFine;
    int nElementsYFine;
    int nMembEl;
    int factorX;
    int factorY;
    mutable bool membraneVisited;
};


  template <typename PHYSICS, bool isSubdomain>
  struct ElementPointerSwitch
  {
    typedef typename PHYSICS::ElementPointer type;
  };

  template<typename PHYSICS>
  struct ElementPointerSwitch<PHYSICS,true>
  {
    typedef typename PHYSICS::SubDomainElementPointer type;
  };

  template<typename PHYSICS, typename Grid, bool isSubdomain>
  struct get_element_pointer
  {
    static typename ElementPointerSwitch<PHYSICS,isSubdomain>::type
    get(const Grid& grid, typename PHYSICS::ElementPointer& ep)
    {
      return ep;
    }
  };

  template<typename PHYSICS, typename Grid>
  struct get_element_pointer<PHYSICS,Grid,true>
  {
    static typename ElementPointerSwitch<PHYSICS,true>::type
    get(const Grid& grid, typename PHYSICS::ElementPointer& ep)
    {
      return grid.template subDomainEntityPointer<0>(*ep);
    }
  };


/**!
 * \note Even newer version of all the grid data transfer functions defined above
 *
 * \brief This gridfunction acts as if it was defined on a gridview FineGV, but actually only interpolates the
 * values from a coarser gridview CoarseGV. Useful when transferring the system equilibrium state computed
 * on a coarse grid to a much larger (number-of-DOFs-wise) grid
 *
 * This method does not take any assumptions about the layout of coarse and fine grid other than there is
 * a unique mapping from a fine entity to a coarse entity, which is the case if the fine grid is a (bisection)
 * refinement of the coarse grid. This mapping has to be provided by the user!
 */
template<typename FineGV, typename CoarseGF, typename PHYSICS, bool isSubdomainGridview = false>
class Ax1CoarseToFineGridMap
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<FineGV,
                                       typename CoarseGF::Traits::RangeFieldType,
                                       CoarseGF::Traits::RangeType::dimension,
                                       typename CoarseGF::Traits::RangeType>,
                                       Ax1CoarseToFineGridMap<FineGV,CoarseGF,PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<FineGV,   // grid view type
                                            typename CoarseGF::Traits::RangeFieldType, // image field type (double)
                                            CoarseGF::Traits::RangeType::dimension,                                        // number of components of image (k)
                                            typename CoarseGF::Traits::RangeType // image type (Dune::FieldVector<double, k>)
                                           > Traits;

  typedef typename CoarseGF::Traits::GridViewType CoarseGV;
  typedef typename CoarseGV::template Codim<0>::Entity CoarseEntity;
  typedef typename CoarseGV::template Codim<0>::EntityPointer CoarseEntityPointer;
  typedef typename CoarseGV::template Codim<0>::
    template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator CoarseElementLeafIterator;

  typedef Dune::SingleCodimSingleGeomTypeMapper<FineGV, 0> ElementMapper;

  //! constructor
  Ax1CoarseToFineGridMap (const FineGV& fineGV_, const std::vector<Dune::shared_ptr<CoarseGF> >& coarseGF_,
      const PHYSICS& physics_, const std::map<int,CoarseEntityPointer>& coarseToFineElementMap_,
      std::vector<int>& coarseGroups_, bool doInterpolate_ = false)
    : fineGV(fineGV_),
      coarseGF(coarseGF_),
      coarseGV(coarseGF_[0]->getGridView()),
      physics(physics_),
      coarseToFineElementMap(coarseToFineElementMap_),
      coarseGroups(coarseGroups_),
      countEvaluated(0),
      // Initialize pointer to current coarse entity
      doInterpolate(doInterpolate_),
      verbose(false),
      elemMapper(fineGV)
  {
  }

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    int elementIndex = elemMapper.map(e);

    // This is the coarse entity in which the current fine entity is contained
    CoarseEntityPointer cep = coarseToFineElementMap.at(elementIndex);

    // xnew has now local coordinates with respect to cep!
    typename Traits::DomainType xnew = cep->geometry().local(e.geometry().global(x));

    // Paranoia check
    typename Traits::DomainType diffCheck = e.geometry().global(x);
    diffCheck -= cep->geometry().global(xnew);
    //debug_jochen << "  - coarse/fine grid position: " << cep->geometry().global(xnew) << " / "
    //    << e.geometry().global(x) << std::endl;
    assert(diffCheck.two_norm() < 1e-6);

    assert(coarseGF.size() <= 2);

    // Do interpolation here, coarseGF is now a vector of gridfunctions!
    // Standard evaluation, no interpolation
    if(coarseGF.size() == 1)
    {
      // Evaluate the coarse gridfunction
      coarseGF[0]->evaluate(*cep,xnew,y);
    } else {
      // When on membrane: interpolate values!

      // Hardcoded version for maximum 2 membrane groups
      double perm1 = physics.getGroupPermittivity(coarseGroups[0]);
      double perm2 = physics.getGroupPermittivity(coarseGroups[1]);

      double d_perm = perm2-perm1;

      // How the interpolation should generally work (at the moment hardcoded for 2 gridfunctions):
      // Check x-coordinate to determine between which two membrane groups this element lies; then it is also easy to
      // pick exactly those two gridfunctions that we want to interpolate here, that means we can even allow more than 2
      // coarse GF in the vector!
      // - Get group index of membrane intervals left and right
      // - Get (possibly smoothed) membrane permittivity for this x-coordinate, regardless if we are on the membrane or not
      // - Interpolate the two picked GFs (using the left/right group indices) using their membrane group permittivity and
      //   my membrane group permittivity

      // Find the membrane element with the same x-coordinate than the current element.
      // We need the (membrane) permittivity at that very element!
      double correspondingMembranePermittivity = -1;
      for(typename PHYSICS::MIterator mit = physics.mBegin(); mit != physics.mEnd(); ++mit)
      {
        // Compare x coordinates
        if(std::abs(mit->geometry().center()[0] - e.geometry().center()[0]) < 1e-6)
        {
          correspondingMembranePermittivity = physics.getPermittivity(*mit->inside());
          break;
        }
      }
      assert(correspondingMembranePermittivity > -1);

      double scale = (correspondingMembranePermittivity - perm1) / d_perm;

      if(verbose && physics.isMembrane(e))
      {
        debug_jochen << "Element @ " << e.geometry().center() << std::endl;
        debug_jochen << "My membrane group permittivity = " << correspondingMembranePermittivity << " => scale = " << scale << std::endl;
      }

      if(! doInterpolate)
      {
        scale = std::round(scale);
      }

      assert(scale >= 0.0 && scale <= 1.0);

      // Do linear interpolation on y
      typename Traits::RangeType y_interp1(0.0);
      coarseGF[0]->evaluate(*cep,xnew,y_interp1);

      y = y_interp1;

      typename Traits::RangeType y_interp2(0.0);
      coarseGF[1]->evaluate(*cep,xnew,y_interp2);
      y_interp2 -= y_interp1;
      y_interp2 *= scale;

      y += y_interp2;
    }

    if(verbose && physics.isMembrane(e))
      debug_jochen << "Fine element #" << physics.getElementIndex(e) << "[" << countEvaluated << "], "
        << "transfer value " << y << " @ coarse: " << cep->geometry().global(xnew)
        << " --> fine: " << e.geometry().global(x) << std::endl;

  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return fineGV;
  }

  void setVerbose(bool verbose_)
  {
    verbose = verbose_;
  }

private:



  const FineGV& fineGV;
  const std::vector<Dune::shared_ptr<CoarseGF> >& coarseGF;
  const CoarseGV& coarseGV;
  const PHYSICS& physics;
  const std::map<int,CoarseEntityPointer>& coarseToFineElementMap;

  mutable int countEvaluated;

  std::vector<int> coarseGroups;
  bool doInterpolate;
  bool verbose;
  ElementMapper elemMapper;
};


#endif /* DUNE_AX1_COARSEGRIDINTERPOLATIONGRIDFUNCTION */
