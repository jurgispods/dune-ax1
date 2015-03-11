/*
 * ax1_subgridview.hh
 *
 *  Created on: Feb 14, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_SUBGRIDVIEW_HH
#define DUNE_AX1_SUBGRIDVIEW_HH

#include <dune/geometry/type.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>

template<typename SubGV, typename HostElementMapper>
class Ax1SubGridView
{
  public:
    typedef Ax1SubGridView<SubGV,HostElementMapper> ThisType;

    //typedef typename SubGV::GridViewImp GridViewImp;
    typedef SubGV GridViewImp;

    /** \brief Traits class */
    typedef typename SubGV::Traits Traits;

    /** \brief type of the grid */
    typedef typename Traits::Grid Grid;

    /** \brief type of the index set */
    typedef typename Traits::IndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Traits::Intersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Traits::IntersectionIterator IntersectionIterator;

    /** \brief type of the collective communication */
    typedef typename Traits::CollectiveCommunication CollectiveCommunication;

    /** \brief A struct that collects all associated types of one implementation
                     from the Traits class.
    */
    template< int cd >
    struct Codim {
     /** \brief type of iterator returned by the grid view */
     typedef typename Traits :: template Codim<cd> :: Iterator Iterator;

     /** \brief type of corresponding entity pointer */
     typedef typename Traits :: template Codim<cd> :: EntityPointer EntityPointer;

     /** \brief type of corresponding entity */
     typedef typename Traits :: template Codim<cd> :: Entity Entity;

     /** \brief type of the geometry implementation */
     typedef typename Traits :: template Codim<cd> :: Geometry Geometry;

     /** \brief type of the implementation for local geometries */
     typedef typename Traits :: template Codim<cd> :: LocalGeometry LocalGeometry;

     /** \brief Define types needed to iterate over entities of a given partition type */
     template< Dune::PartitionIteratorType pit >
     struct Partition
     {
       /** \brief iterator over a given codim and partition type */
       typedef typename Traits :: template Codim< cd >
         :: template Partition< pit > :: Iterator Iterator;
     };
    }; //: public Traits :: template Codim<cd> {};

    enum {  //! \brief Export if this grid view is conforming */
           conforming = Traits :: conforming };

    /** \brief type used for coordinates in grid */
    typedef typename Grid::ctype ctype;

    enum { //! \brief The dimension of the grid
           dimension = Grid :: dimension };

    enum { //! \brief The dimension of the world the grid lives in.
           dimensionworld = Grid :: dimensionworld };

  public:
    /** \brief constructor (engine concept) */
    Ax1SubGridView (const SubGV &impl, HostElementMapper& hostElementMapper_ )
    : impl_( impl ),
     hostElementMapper(hostElementMapper_)
    {}

    /** \brief Copy constructor */
    Ax1SubGridView ( const ThisType &other, HostElementMapper& hostElementMapper_)
    : impl_( other.impl_ ),
     hostElementMapper(hostElementMapper_)
    {}

    /** \brief assignment operator */
    ThisType &operator= ( const ThisType &other )
    {
     impl_ = other.impl_;
     hostElementMapper = other.hostElementMapper;
     return *this;
    }

  public:
    /** \brief obtain a const reference to the underlying hierarchic grid */
    const Grid &grid () const
    {
      return asImp().grid();
    }

    /** \brief obtain the index set */
    const IndexSet &indexSet () const
    {
      return asImp().indexSet();
    }

    /** \brief obtain numer of entities in a given codimension */
    int size ( int codim ) const
    {
      return asImp().size( codim );
    }

    /** \brief obtain number of entities with a given geometry type */
    int size ( const Dune::GeometryType &type ) const
    {
      return asImp().size( type );
    }

    /** @brief Return true if the given entity is contained in this grid view
     * @todo Currently we call the implementation on the IndexSet.  This may lead to suboptimal efficiency.
     *
     * \note If the input element e is not an element of the grid, then
     *       the result of contains() is undefined.
    */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      if(asImp().indexSet().contains(e))
      {
        if(hostElementMapper.contains(indexSet().index(e)))
        {
          //Check for consistent geometry of subgrid element and host element
          return (hostElementMapper.map(indexSet().index(e))->geometry().center() == e.geometry().center());
        }
      }
      return false;
    }

    /** \brief obtain begin iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator begin () const
    {
      return asImp().template begin<cd>();
    }

    /** \brief obtain end iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator end () const
    {
      return asImp().template end<cd>();
    }

    /** \brief obtain begin iterator for this view */
    template< int cd , Dune::PartitionIteratorType pitype >
    typename Codim< cd > :: template Partition< pitype > :: Iterator
    begin () const
    {
      return asImp().template begin<cd,pitype>();
    }

    /** \brief obtain end iterator for this view */
    template< int cd, Dune::PartitionIteratorType pitype >
    typename Codim< cd > :: template Partition< pitype > :: Iterator
    end () const
    {
      return asImp().template end<cd,pitype>();
    }

    /** \brief obtain begin intersection iterator with respect to this view */
    IntersectionIterator
    ibegin ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return asImp().ibegin(entity);
    }

    /** \brief obtain end intersection iterator with respect to this view */
    IntersectionIterator
    iend ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return asImp().iend(entity);
    }

    /** \brief obtain collective communication object */
    const CollectiveCommunication &comm () const
    {
      return asImp().comm();
    }

    /** \brief Return size of the overlap region for a given codim on the grid view.  */
    int overlapSize(int codim) const
    {
     return asImp().overlapSize(codim);
    }

    /** \brief Return size of the ghost region for a given codim on the grid view.  */
    int ghostSize(int codim) const
    {
     return asImp().ghostSize(codim);
    }

    /** communicate data on this view */
    template< class DataHandleImp, class DataType >
    void communicate ( Dune::CommDataHandleIF< DataHandleImp, DataType > &data,
                       Dune::InterfaceType iftype,
                       Dune::CommunicationDirection dir ) const
    {
      asImp().communicate(data,iftype,dir);
    }

    // ====================== BEGIN dune-ax1 custom subgrid stuff ============================================

    const typename Codim<0>::Entity& getHostElement(const typename Codim<0>::Entity& e) const
    {
      // TODO check if element is contained?
      return *hostElementMapper.map(indexSet().index(e));
    }

    const typename Codim<0>::Entity& getHostElement(const unsigned int subIndex) const
    {
      // TODO check if element is contained?
      return *hostElementMapper.map(subIndex);
    }

    // ====================== END dune-ax1 custom subgrid stuff =============================================

#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
  public:
#else
  protected:
    // give the GridDefaultImplementation class access to the realImp
    friend class Dune::GridDefaultImplementation< Grid::dimension, Grid::dimensionworld, typename Grid::ctype, typename Grid::GridFamily >;
#endif
    // type of underlying implementation, for internal use only
    typedef GridViewImp Implementation;

    //! return reference to the real implementation
    Implementation &impl () { return impl_; }
    //! return reference to the real implementation
    const Implementation &impl () const { return impl_; }

  protected:
    Implementation impl_;

    GridViewImp& asImp ()
    {
      return impl_;
    }

    const GridViewImp& asImp () const
    {
      return impl_;
    }

  private:
    HostElementMapper& hostElementMapper;

};

class NoSubGridView
{
  public:
    template<typename Element>
    const Element& getHostElement(const Element& e)
    {}
};

#endif /* DUNE_AX1_SUBGRIDVIEW_HH */
