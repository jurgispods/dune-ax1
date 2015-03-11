/*
 * acme2_cyl_geometrywrapper.hh
 *
 *  Created on: Oct 24, 2012
 *      Author: jpods
 */

#ifndef DUNE_ACME2_CYL_GEOMETRYWRAPPER_HH
#define DUNE_ACME2_CYL_GEOMETRYWRAPPER_HH

#include <dune/grid/common/geometry.hh>

//#define USE_CYL_GEOMETRY_WRAPPER

#ifndef USE_CYL_GEOMETRY_WRAPPER

/**
 * Cylinder Geometry, derives from a Cartesian Dune::Geometry
 */
template<typename Geometry>
class Acme2CylGeometry : public Geometry
{

public:

  typedef Geometry BaseT;

  typedef typename BaseT::ctype ctype;
  static const int dimension = BaseT::dimension;
  static const int dimensionworld = BaseT::dimensionworld;
  static const int mydimension = BaseT::mydimension;
  static const int coorddimension = BaseT::coorddimension;

private:

  typedef Dune::FieldVector<ctype,coorddimension> GlobalCoords;
  typedef Dune::FieldVector<ctype,mydimension> LocalCoords;

public:


  Acme2CylGeometry(const Geometry& geometry)
    : BaseT(geometry)
      //left_bottom(0.0),
      //right_top(1.0),
      //r1(BaseT::global(LocalCoords(0.0))[1]),
      //r2(BaseT::global(LocalCoords(1.0))[1])
  {
    dune_static_assert(dimension == 2, "This implementation for cylinder coordinates works in 2D only");
    dune_static_assert(dimensionworld == 2, "This implementation for cylinder coordinates works in 2D only");
    dune_static_assert(coorddimension == 2, "This implementation for cylinder coordinates works in 2D only");
    dune_static_assert(mydimension <= 2, "mydimension cannot be greater than 2");
  }


  ctype integrationElement(const LocalCoords& local) const
  {
    // Factor 2*pi*r
    const ctype factor = 2 * con_pi * BaseT::global(local)[1];

    return factor * BaseT::integrationElement(local);
  }

  ctype volume() const
  {
    // Factor pi*(r2+r1)
    const ctype factor = con_pi * (BaseT::corner(BaseT::corners()-1)[1]
                             + BaseT::corner(0)[1]);

    return factor * BaseT::volume();
    //return con_pi * (r2+r1) * BaseT::volume();
  }


  /**
   * Method to check the orientation of a codim-1 entity
   * @return
   */
  bool isOrientedAxially() const
  {
    if(mydimension != 1) return false;

    ctype r_diff = this->corner(this->corners()-1)[1];
    r_diff -= this->corner(0)[1];

    if(std::abs(r_diff) < 1e-6)
      return true;
    else
      return false;
  }

private:
  //const LocalCoords left_bottom;
  //const LocalCoords right_top;
  //const ctype r1;
  //const ctype r2;

};

/**
 * A specialization of an ElementGeometry for the cylinder case
 */
template<typename EG>
class Acme2CylElementGeometry : public EG
{
  public:
    typedef EG BaseT;
    typedef typename BaseT::Entity Entity;
    typedef typename BaseT::Geometry HostGeometry;
    typedef Acme2CylGeometry<HostGeometry> Geometry;

    Acme2CylElementGeometry(const EG& eg_)
    : BaseT(eg_),
      cylGeo(eg_.geometry())
    {}

    /**
     * Override base class' geometry() method and return an Acme2CylGeometry instead
     */
    const Geometry& geometry () const
    {
      return cylGeo;
    }

  private:
    Geometry cylGeo;
};

namespace detail
{
  template<typename T, bool isMultiDomainClass>
  struct BaseTypes
  {
     typedef void IntersectionType;
     typedef void EntityWrapperType;
  };

   template<typename T>
   struct BaseTypes<T,true>
   {
     typedef typename T::Intersection IntersectionType;
     typedef typename T::EntityWrapper EntityWrapperType;
   };

   // SFINAE, bitch!
   template <typename T>
   struct has_typedef_EntityWrapper
   {
     // Types "yes" and "no" are guaranteed to have different sizes,
     // specifically sizeof(yes) == 1 and sizeof(no) == 2.
     typedef char yes[1];
     typedef char no[2];

     template <typename C>
     static yes& test(typename C::EntityWrapper*);

     template <typename>
     static no& test(...);

     // If the "sizeof" the result of calling test<T>(0) would be equal to the sizeof(yes),
     // the first overload worked and T has a nested type named foobar.
     static const bool value = sizeof(test<T>(0)) == sizeof(yes);
   };
}



/**
 * A specialization of an IntersectioGeometry for the cylinder case
 */
template<typename IG>
class Acme2CylIntersectionGeometry : public IG
{
  public:
    typedef IG BaseT;

    typedef typename BaseT::Geometry HostGeometry;
    typedef Acme2CylGeometry<HostGeometry> Geometry;
    typedef typename BaseT::LocalGeometry LocalGeometry;
    typedef typename BaseT::Entity Entity;
    typedef typename BaseT::EntityPointer EntityPointer;
    typedef typename BaseT::ctype ctype;

    enum { dimension=BaseT::dimension };
    enum { dimensionworld=BaseT::dimensionworld };

    // This will not work, as BaseT is of complete type, but SkeletonIntersectionWrapper is not!
//    typedef typename detail::BaseTypes<BaseT,Dune::IsBaseOf<
//        BaseT,Dune::PDELab::MultiDomain::SkeletonIntersectionWrapper> >::IntersectionType Intersection;
//    typedef typename detail::BaseTypes<BaseT,Dune::IsBaseOf<
//        BaseT,Dune::PDELab::MultiDomain::SkeletonIntersectionWrapper> >::EntityWrapperType EntityWrapper;

    // Enable MultiDomain typedef if they are available! (void otherwise)
    typedef typename detail::BaseTypes<
        BaseT,detail::has_typedef_EntityWrapper<BaseT>::value>::EntityWrapperType EntityWrapper;
    typedef typename detail::BaseTypes<
        BaseT,detail::has_typedef_EntityWrapper<BaseT>::value>::IntersectionType Intersection;


    //! \todo Please doc me!
    Acme2CylIntersectionGeometry (const IG& ig_)
        : BaseT(ig_),
          cylGeo(ig_.geometry())
    {}

    /**
     * Override base class' geometry() method and return an Acme2CylGeometry instead
     */
    const Geometry& geometry () const
    {
      return cylGeo;
    }

  private:
    Geometry cylGeo;
};

/**
 * Template magic for switching between 2D Cartesian and 2D Cylinder coordinates
 */
struct Acme2CylGeometrySwitch
{
    /**
     * Some TMP magic to switch between the used ElementGeometry wrappers
     */
    template<typename EG, bool useCylinderCoords = USE_CYLINDER_COORDINATES>
    struct ElementGeometrySwitch
    {
      typedef EG type;
    };

    template<typename EG>
    struct ElementGeometrySwitch<EG,true>
    {
      typedef Acme2CylElementGeometry<EG> type;
    };

    template<typename IG, bool useCylinderCoords = USE_CYLINDER_COORDINATES>
    struct IntersectionGeometrySwitch
    {
      typedef IG type;
    };

    template<typename IG>
    struct IntersectionGeometrySwitch<IG,true>
    {
      typedef Acme2CylIntersectionGeometry<IG> type;
    };

    template<typename GEO, bool useCylinderCoords = USE_CYLINDER_COORDINATES>
    struct GeometrySwitch
    {
      typedef GEO type;
    };

    template<typename GEO>
    struct GeometrySwitch<GEO,true>
    {
      typedef Acme2CylGeometry<GEO> type;
    };
};


#else

/**
 * Geometry wrapper, default version forwards every call to the wrapped geometry
 */
template<typename Geometry>
class Acme2CylGeometryWrapper
{

public:

  typedef typename Geometry::ctype ctype;
  static const int dimension = Geometry::dimension;
  static const int dimensionworld = Geometry::dimensionworld;
  static const int mydimension = Geometry::mydimension;
  static const int coorddimension = Geometry::coorddimension;

private:

  typedef Dune::FieldVector<ctype,coorddimension> GlobalCoords;
  typedef Dune::FieldVector<ctype,mydimension> LocalCoords;

public:


  Acme2CylGeometryWrapper(const Geometry& wrappedGeometry)
    : _wrappedGeometry(wrappedGeometry),
      r1(_wrappedGeometry.global(LocalCoords(0.0))[1]),
      r2(_wrappedGeometry.global(LocalCoords(1.0))[1])
  {
    dune_static_assert(dimension == 2, "This implementation for cylinder coordinates works in 2D only");
    dune_static_assert(dimensionworld == 2, "This implementation for cylinder coordinates works in 2D only");
    dune_static_assert(coorddimension == 2, "This implementation for cylinder coordinates works in 2D only");
    dune_static_assert(mydimension <= 2, "mydimension cannot be greater than 2");
  }

  Dune::GeometryType type() const {
    return _wrappedGeometry.type();
  }

  int corners() const {
    return _wrappedGeometry.corners();
  }

  bool affine() const {
    return _wrappedGeometry.affine();
  }

  GlobalCoords corner(int i) const {
    return _wrappedGeometry.corner(i);
  }

  GlobalCoords global(const LocalCoords& local) const {
    return _wrappedGeometry.global(local);
  }

  LocalCoords local(const GlobalCoords& global) const {
    return _wrappedGeometry.local(global);
  }

  bool checkInside(const LocalCoords& local) const {
    return _wrappedGeometry.checkInside(local);
  }

  ctype integrationElement(const LocalCoords& local) const {

    ctype r = _wrappedGeometry.global(local)[1];
    ctype h = _wrappedGeometry.global(local)[0];

    // Factor 2*pi*(r-r0)
    ctype factor = 2 * con_pi * r;

    if(isOrientedAxially())
    {
      debug_jochen << " => AXIAL, r=" << r << ", h=" << h << std::endl;
      debug_jochen << "    factor " << factor << std::endl;
      debug_jochen << "    original integration element " << _wrappedGeometry.integrationElement(local) << std::endl;
      debug_jochen << "    --> integration element " << (factor * _wrappedGeometry.integrationElement(local)) << std::endl;
    }

    return factor * _wrappedGeometry.integrationElement(local);
  }

  ctype volume() const {
    // Factor pi*(r2+r1) ??!?
    //ctype factor = con_pi * (_wrappedGeometry.corner(_wrappedGeometry.corners()-1)[1]
    //                           + _wrappedGeometry.corner(0)[1]);
    const ctype factor = con_pi * (r2+r1);

    return factor * _wrappedGeometry.volume();
  }

  GlobalCoords center() const {
    return _wrappedGeometry.center();
  }

  const Dune::FieldMatrix<ctype,mydimension,coorddimension>&
  jacobianTransposed(const LocalCoords& local) const {
    return _wrappedGeometry.jacobianTransposed(local);
  }

  const Dune::FieldMatrix<ctype,coorddimension,mydimension>&
  jacobianInverseTransposed(const LocalCoords& local) const {
    return _wrappedGeometry.jacobianInverseTransposed(local);
  }

  /**
   * Method to check the orientation of a codim-1 entity
   * @return
   */
  bool isOrientedAxially() const
  {
    if(mydimension != 1) return false;

    ctype r_diff = _wrappedGeometry.corner(_wrappedGeometry.corners()-1)[1];
    r_diff -= _wrappedGeometry.corner(0)[1];

    if(std::abs(r_diff) < 1e-6)
      return true;
    else
      return false;
  }

private:
  const Geometry _wrappedGeometry;
  const ctype r1;
  const ctype r2;

};

/**
 * A wrapper for PDELab's ElementGeometry class
 */
template<typename EG>
class Acme2CylElementGeometryWrapper
{
  public:
    typedef typename EG::Entity Entity;
    typedef typename EG::Geometry HostGeometry;
    typedef Acme2CylGeometryWrapper<HostGeometry> Geometry;

    Acme2CylElementGeometryWrapper(const EG& eg_)
    : eg(eg_), geometryWrapper(eg.geometry())
    {}

    Geometry geometry () const
    {
      return geometryWrapper;
    }

    const Entity& entity () const
    {
      return eg.entity();
    }

    const Entity& hostEntity () const
    {
      return eg.entity();
    }

    //TODO Add missing multidomain methods


  private:
    const EG& eg;
    Geometry geometryWrapper;
};



/**
 * A wrapper for multidomains's IntersectionGeometry classes
 */
template<typename IG>
class Acme2CylIntersectionGeometryWrapper
{
  public:
    typedef typename IG::Intersection Intersection;

    typedef typename IG::Geometry HostGeometry;
    typedef Acme2CylGeometryWrapper<HostGeometry> Geometry;

    typedef typename IG::LocalGeometry LocalGeometry;

    typedef typename IG::Entity Entity;

    typedef typename IG::EntityPointer EntityPointer;

    // Multidomain intersection wrapper typedef
    typedef typename IG::EntityWrapper EntityWrapper;

    typedef typename IG::ctype ctype;
    enum { dimension=IG::dimension };
    enum { dimensionworld=IG::dimensionworld };

      //! \todo Please doc me!
    Acme2CylIntersectionGeometryWrapper (const IG& ig_)
        : ig(ig_), geometryWrapper(ig.geometry())
    {}

    //! \todo Please doc me!
    int insideDomainIndex() const
    {
      return ig.insideDomainIndex();
    }

    //! \todo Please doc me!
    int outsideDomainIndex() const
    {
      ig.outsideDomainIndex();
    }

    //! return true if intersection is with interior or exterior boundary (see the cases above)
    bool boundary () const
    {
      return ig.boundary();
    }

    /**
     \brief Identifier for boundary segment from macro grid.

     One can attach a boundary Id to a boundary segment on the macro
     grid. This Id will also be used for all fragments of these
     boundary segments.

     The numbering is defined as:
     - Id==0 for all intersections without boundary()==false
     - Id>=0 for all intersections without boundary()==true

     The way the Identifiers are attached to the grid may differ
     between the different grid implementations.

    */
    int boundaryId () const
    {
      return ig.boundaryId();
    }

    //! @brief return true if intersection is shared with another element.
    bool neighbor () const
    {
      return ig.neighbor();
    }

    /*! @brief geometrical information about this intersection in local
    coordinates of the inside() entity.

    This method returns a Geometry object that provides a mapping from
    local coordinates of the intersection to local coordinates of the
    inside() entity.
    */
    LocalGeometry geometryInInside () const
    {
      return ig.geometryInInside();
    }

    /*! @brief geometrical information about this intersection in local
    coordinates of the outside() entity.

    This method returns a Geometry object that provides a mapping from
    local coordinates of the intersection to local coordinates of the
    outside() entity.
    */
    LocalGeometry geometryInOutside () const
    {
      return ig.geometryInOutside();
    }

    /*! @brief geometrical information about this intersection in global coordinates.

    This method returns a Geometry object that provides a mapping from
    local coordinates of the intersection to global (world) coordinates.
    */
    Geometry geometry () const
    {
      //return ig.geometry();
      return geometryWrapper;
    }

    //! Local number of codim 1 entity in the inside() Entity where intersection is contained in
    int indexInInside () const
    {
      return ig.indexInInside ();
    }

    //! Local number of codim 1 entity in outside() Entity where intersection is contained in
    int indexInOutside () const
    {
      return ig.indexInOutside ();
    }

    /*! @brief Return an outer normal (length not necessarily 1)

    The returned vector may depend on local position within the intersection.
    */
    Dune::FieldVector<ctype, dimensionworld> outerNormal (const Dune::FieldVector<ctype, dimension-1>& local) const
    {
      return ig.outerNormal(local);
    }

    /*! @brief return outer normal scaled with the integration element
    @copydoc outerNormal
    The normal is scaled with the integration element of the intersection. This
    method is redundant but it may be more efficient to use this function
    rather than computing the integration element via intersectionGlobal().
    */
    Dune::FieldVector<ctype, dimensionworld> integrationOuterNormal (const Dune::FieldVector<ctype, dimension-1>& local) const
    {
      DUNE_THROW(Dune::NotImplemented, "Implement me!");
      //FIXME use correct integration element!
      //return ig.integrationOuterNormal(local);
    }

    /*! @brief Return unit outer normal (length == 1)

    The returned vector may depend on the local position within the intersection.
    It is scaled to have unit length.
    */
    Dune::FieldVector<ctype, dimensionworld> unitOuterNormal (const Dune::FieldVector<ctype, dimension-1>& local) const
    {
      return ig.unitOuterNormal(local);
    }

    /*! @brief Return unit outer normal (length == 1)

    The returned vector may depend on the local position within the intersection.
    It is scaled to have unit length.
    */
    Dune::FieldVector<ctype, dimensionworld> centerUnitOuterNormal () const
    {
    return ig.centerUnitOuterNormal();
    }

    /*! @brief return EntityPointer to the Entity on the inside of this
    intersection. That is the Entity where we started this .
    */
    EntityPointer inside() const
    {
      //TODO Add warning or assure the user is not trying to call this method,
      // as using the inside's geometry integration element / volume is not valid!
      return ig.inside();
    }

    /*! @brief return EntityPointer to the Entity on the inside of this
    intersection. That is the Entity where we started this .
    */
    EntityPointer insideHostEntity() const
    {
      DUNE_THROW(Dune::Exception,"This should never be called.");
      return ig.inside();
    }

    /*! @brief return EntityPointer to the Entity on the outside of this
    intersection. That is the neighboring Entity.

    @warning Don't call this method if there is no neighboring Entity
    (neighbor() returns false). In this case the result is undefined.
    */
    EntityPointer outside() const
    {
      //TODO Add warning or assure the user is not trying to call this method,
      // as using the outside's geometry integration element / volume is not valid!
      return ig.outside();
    }

    //! \todo Please doc me!
    const Intersection& intersection () const
    {
      return ig.intersection();
    }

    unsigned int intersectionIndex() const
    {
      return ig.intersectionIndex();
    }

    //TODO Add missing multidomain methods

  private:
    const IG& ig;
    Geometry geometryWrapper;
};

/**
 * Template magic for switching between 2D Cartesian and 2D Cylinder coordinates
 */
struct Acme2CylGeometrySwitch
{
    /**
     * Some TMP magic to switch between the used ElementGeometry wrappers
     */
    template<typename EG, bool useCylinderCoords = USE_CYLINDER_COORDINATES>
    struct ElementGeometrySwitch
    {
      typedef EG type;
    };

    template<typename EG>
    struct ElementGeometrySwitch<EG,true>
    {
      typedef Acme2CylElementGeometryWrapper<EG> type;
    };

    template<typename IG, bool useCylinderCoords = USE_CYLINDER_COORDINATES>
    struct IntersectionGeometrySwitch
    {
      typedef IG type;
    };

    template<typename IG>
    struct IntersectionGeometrySwitch<IG,true>
    {
      typedef Acme2CylIntersectionGeometryWrapper<IG> type;
    };

    template<typename GEO, bool useCylinderCoords = USE_CYLINDER_COORDINATES>
    struct GeometrySwitch
    {
      typedef GEO type;
    };

    template<typename GEO>
    struct GeometrySwitch<GEO,true>
    {
      typedef Acme2CylGeometryWrapper<GEO> type;
    };
};

// Add helper template for the rare case where an explicit Acme2CylGeometry is requested
template<typename Geometry>
struct Acme2CylGeometry : public Acme2CylElementGeometryWrapper<Geometry>
{
};

#endif


#endif /* DUNE_ACME2_CYL_GEOMETRYWRAPPER_HH */


