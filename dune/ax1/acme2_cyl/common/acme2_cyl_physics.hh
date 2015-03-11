#ifndef DUNE_AX1_ACME2CYL_PHYSICS_HH
#define DUNE_AX1_ACME2CYL_PHYSICS_HH

#include <vector>
#include <iomanip>

#include <dune/pdelab/gridfunctionspace/intersectionindexset.hh>

#include <dune/ax1/common/ax1_subgrid_hostelement_mapper.hh>
#include <dune/ax1/common/ax1_lfs_tools.hh>
#include <dune/ax1/common/ax1_groupmapper.hh>
#include <dune/ax1/common/ax1_interfaceiterator.hh>
#include <dune/ax1/common/electrolyte.hh>
#include <dune/ax1/common/element_subdomain_mapper.hh>
#include <dune/ax1/common/membrane.hh>
#include <dune/ax1/common/membrane_interface_mapper.hh>
#include <dune/ax1/common/intersectiongridfunction.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_geometrywrapper.hh>

#include <dune/ax1/common/channelgridfunction.hh>

template<typename PHYSICS>
struct is_intersection_in_vector
{
  is_intersection_in_vector(const typename PHYSICS::ElementIntersection& is_, const PHYSICS& physics_)
  : is(is_),
    physics(physics_)
  {}

  bool operator() (const typename PHYSICS::ElementIntersection& other)
  {
    //debug_jochen << "Comparing intersections for equality: "
    //    << physics.getIntersectionIndex(is) << " - "
    //    << physics.getIntersectionIndex(other) << std::endl;
    return physics.getIntersectionIndex(is) == physics.getIntersectionIndex(other);
  }

  const typename PHYSICS::ElementIntersection& is;
  const PHYSICS& physics;
};

template<typename GV, typename T, typename ConfigTraits>
class Physics {
public:

  typedef Physics<GV,T,ConfigTraits> This;

  // General typedefs
  typedef GV GridView;
  typedef T FieldType;
  typedef ConfigTraits Traits;

  typedef typename GV::Grid::SubDomainGrid               SubGrid;
  typedef typename GV::Grid::SubDomainGrid::LeafGridView SubGV;

  // Grid-specific typedefs
  typedef typename GV::template Codim<0>::Entity Element;
  typedef typename GV::template Codim<0>::EntityPointer ElementPointer;
  typedef typename GV::template Codim<0>
    ::template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator ElementIterator;
  typedef typename GV::template Codim<0>::Iterator ElementIterator_All;
  typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

  typedef typename GV::Intersection ElementIntersection;
  typedef typename Element::LeafIntersectionIterator ElementIntersectionIterator;
  
  // Dune::SubGrid typedefs
  typedef typename SubGV::template Codim<0>::Entity SubDomainElement;
  typedef typename SubGV::template Codim<0>::EntityPointer SubDomainElementPointer;
  typedef typename SubGV::template Codim<0>
    ::template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator SubDomainElementIterator;
  typedef typename SubGV::template Codim<0>::Iterator SubDomainElementIterator_All;
  typedef Dune::SingleCodimSingleGeomTypeMapper<SubGV, 0> SubDomainElementMapper;

  typedef typename SubDomainElement::LeafIntersectionIterator SubDomainElementIntersectionIterator;
  typedef typename SubDomainElementIntersectionIterator::Intersection SubDomainElementIntersection;

  typedef typename GV::Grid::LocalIdSet::IdType IdType;
  typedef typename SubGV::Grid::LocalIdSet::IdType SubDomainIdType;

  typedef std::vector<ElementIntersectionIterator> MembraneInterfaces;
  typedef std::map<int,ElementIntersectionIterator> MembraneInterfaceMap;

  //typedef typename MembraneInterfaces::const_iterator MIterator;
  typedef Ax1InterfaceIterator<const MembraneInterfaces, ElementIntersectionIterator, const ElementIntersection> MIterator;

  typedef SubGridHostElementMapper<ElementPointer,ElementMapper> HostElementMapper;

  typedef typename Membrane<T>::ChannelSet ChannelSet;

  typedef Dune::FieldVector<FieldType,GV::dimension> CoordType;
  typedef std::unordered_map<CoordType, double, std::hash<CoordType>, field_vector_equal>
        NodeVolumesMap;


  Physics ( const GV& gv_, const SubGV& elecGV_, const SubGV& membGV_,
           ElementSubdomainMapper& elementSubdomainMapper_,
           Ax1ElementGroupMapper& elementGroupMapper_,
           Electrolyte<T>& electro_, Membrane<T>& memb_, Acme2CylParameters& params_,
           const T& timeScale_, const T& lengthScale_)
  : gv(gv_),
    elecGV(elecGV_),
    membGV(membGV_),
    elemSubdomainMapper(elementSubdomainMapper_),
    elemGroupMapper(elementGroupMapper_),
    electro(electro_),
    membrane(memb_),
    params(params_),
    TIME_SCALE(timeScale_),
    LENGTH_SCALE(lengthScale_),
    elementMapper(gv),
    elecElementMapper(elecGV),
    membElementMapper(membGV),
    nGridElementsPerAxis(GV::dimension, -1),
    nOffsetPerAxis(GV::dimension, -1),
    nSubdomainElementsPerAxis(2),
    membraneInterfaces(0),
    membraneInterfaceMap_primary(0),
    membraneInterfaceMap_bidirectional(0),
    membraneInterfaceMap_self(0),
    partition_x(gv.comm().size()),
    partition_y(gv.comm().size()),
    offset_x(gv.comm().size()),
    offset_y(gv.comm().size()),
    groupNames(params.membrane.getSubKeys())
  {
    // First two groups are always the electrolyte solutions;
    // Insert those at the beginning (intracellular: 0 / extracellular: 1)
    groupNames.insert(groupNames.begin(), "solution_ex");
    groupNames.insert(groupNames.begin(), "solution_in");
    calcNElements();
    calcNodeVolumes();

    initPartition();

    membraneInterfaces = new MembraneInterfaces();
    membraneInterfaceMap_primary = new MembraneInterfaceMap();
    membraneInterfaceMap_bidirectional = new MembraneInterfaceMap();
    membraneInterfaceMap_self = new MembraneInterfaceMap();
    setupMembraneInterfaceMap();

    setupChannels();

    // This sets membrane group-specific permittivities
    initPermittivities();
  }

  //! destructor: delete heap-allocated vector/map memory
  ~Physics()
  {
    //delete[] membraneInterfaces;
    //delete[] membraneInterfaceMap;
    //delete[] membraneInterfaceMap_bidirectional;
  }

  /*! \brief This method tests if an element is a membrane element or not.
   */
  bool isMembrane(const Element& e) const
  {
    //return membGV.contains(e);
    return isMembrane(getSubdomainIndex(getElementIndex(e)));
  }

  /*! \brief This method tests if a subdomain element is a membrane element or not.
   */
  bool isMembrane(const SubDomainElement& e) const
  {
    //return membGV.contains(e);
    return isMembrane(getSubdomainIndex(getElementIndex(e)));
  }

  //! For host grid Intersections: Check if it is a membrane interface
  bool isMembraneInterface(const ElementIntersection& is) const
  {
    if(is.boundary()) // real domain boundary; membrane might be located at the domain boundary
    {
      // This will return true only for oriented horizontally intersections at an outer domain boundary;
      // This compares only with yMemb[0], so it should work only on the lower boundary
      return (std::abs(is.geometry().center()[1]-params.yMemb()[0]) < 1e-6
              && isMembrane(*is.inside())); // TODO Better set back to 'return false' here
    }
    if(! is.neighbor())
      return false; // processor boundary
//    if(is.outside()->partitionType() != Dune::PartitionType::InteriorEntity)
//    {
//      debug_jochen << "Intersection @ " << is.geometry().center()
//          << " with adjacent elements interior (" << is.inside()->geometry().center()
//          << ") and outside (" <<  is.outside()->geometry().center()
//          << ")" << std::endl;
//      return false; // interior/overlap boundary
//    }
    // Exactly one of the neighboring elements is a membrane element
    return (isMembrane(*is.inside()) != isMembrane(*is.outside()));
  }

  //! For SubGridIntersections: Check if it is a membrane interface
  bool isMembraneInterface(const SubDomainElementIntersection& is) const
  {
    // According to our current assumptions, the membrane interface is always located at
    // the boundary of a subdomain!
    if(is.boundary())
    {
      // Check if the corresponding multidomain intersection is a membrane interface
      return isMembraneInterface(gv.grid().multiDomainIntersection(is));
    }
    return false;
  }

  //! \brief This method returns the 'membrane index', the index for one compartment (or patch)
  //! of membrane. As the membrane might be multiple elements thick, but contains only one set
  //! set of channels per element in x-direction, the channels are associated to membrane-electrolyte
  //! interfaces. For a unique mapping, only the lower (cytosol-membrane) interfaces are identified
  //! with the channel sites. This method returns the index of the intersection representing the
  //! membrane patch, whether the handed over intersection is on the 'correct' side of the membrane or
  //! not. That means, it the given parameter 'is' is on the other side (membrane-extracellular interface),
  //! this method will return the index of the opposite membrane interface on the lower membrane surface.
  // TODO Isn't it the other way round? primary = membrane-extracellular interface?
  int getMembraneIndex(const ElementIntersection& is) const
  {
    assert(isMembraneInterface(is));

    // Get the intersection index of the membrane interface the channels are associated with
    // (the correct intersection is provided by membraneInterfaceMap) and return the membrane
    // index belonging to this primary intersection
    return getChannelSet().getMembraneIndex(
        getIntersectionIndex(*membraneInterfaceMap_primary->at(getIntersectionIndex(is))));
  }

  int getMembraneIndex(const SubDomainElementIntersection& is) const
  {
    return getMembraneIndex(gv.grid().multiDomainIntersection(is));
  }

  int getMembraneNumber(const ElementIntersection& is) const
  {
    const std::vector<FieldType> yMemb = params.yMemb();

    for(int i=0; i<yMemb.size(); i++)
    {
      if(std::abs(is.geometry().corner(0)[1]-yMemb[i]) < params.dMemb()+1e-6)
        return i;
    }

    DUNE_THROW(Dune::Exception, "Could not find a matching membrane coordinate for intersection @"
        << is.geometry().center());
  }


  const ElementIntersectionIterator& getNextMembraneInterface(const Element& e) const
  {
    if(isMembrane(e))
    {
      //debug_jochen << "Searching for next primary membrane intersection for membrane element @"
      //    << e.geometry().center() << std::endl;
      ElementIntersectionIterator iit_start = gv.ibegin(e);
      for(ElementIntersectionIterator iit = gv.ibegin(e); iit != gv.iend(e); ++iit)
      {
        //debug_jochen << "=checking if intersection " << iit->geometry().center() << " is a candidate..." << std::endl;
        // Do not enter recursion: Membrane element already found!
        if(isMembraneInterface(*iit))
        {
          //debug_jochen << "=found membrane intersection " << iit->geometry().center() << "!" << std::endl;
          //debug_jochen << "=returning membrane intersection "
          //    << membraneInterfaceMap->at(getIntersectionIndex(*iit))->geometry().center() << "!" << std::endl;
          // Now call membrane interface map in order to always return the primary intersection!
          return membraneInterfaceMap_primary->at(getIntersectionIndex(*iit));
        }
        Acme2CylGeometry<typename ElementIntersectionIterator::Intersection::Geometry> cylGeo(iit->geometry());
        // Save candidate for next recursion call; search direction orthogonal to membrane orientation!
        if(cylGeo.isOrientedAxially())
        {
          iit_start = iit;
          //debug_jochen << "Next candidate intersection @" << iit_start->geometry().center() << std::endl;
        }
      }
      //debug_jochen << "Calling getNext... for intersection @" << iit_start->geometry().center() << std::endl;
      //debug_jochen << " with inside @" << iit_start->inside()->geometry().center() << std::endl;
      //debug_jochen << " and outside @" << iit_start->outside()->geometry().center() << std::endl;
      return getNextMembraneInterface(iit_start);
    }
    DUNE_THROW(Dune::Exception, "Given subdomain element is not a membrane element!");
  }


  const ElementIntersectionIterator& getNextMembraneInterface(const SubDomainElement& se) const
  {
    if(isMembrane(se))
    {
      return getNextMembraneInterface(membGV.grid().multiDomainEntity(se));
    }
    DUNE_THROW(Dune::Exception, "Given subdomain element is not a membrane element!");
  }

  template<typename Pred>
  MIterator my_find_if(MIterator first, MIterator last, Pred pred) const
  {
    while (first!=last) {
      debug_jochen << "Checking is @" << first->geometry().center() << std::endl;
       if (pred(*first)) return first;
       ++first;
     }
     return last;
  }

  bool isIntersectionInMembraneInterfacesVector(const ElementIntersection& is) const
  {
    is_intersection_in_vector<This> pred(is, *this);

//    debug_jochen << "Checking my " << std::endl;
//    bool is_in = my_find_if(mBegin(), mEnd(), pred) != mEnd();
//    debug_jochen << " is_in: " << is_in << std::endl;

//    debug_jochen << "Checking std " << std::endl;
    bool is_in = std::find_if(mBegin(), mEnd(), pred) != mEnd();
//    debug_jochen << " is_in_std: " << is_in << std::endl;

//    if(is_in_std != is_in)
//      DUNE_THROW(Dune::Exception, "KAKCEEEEEEEEEEEEEEEE");

    return is_in;
  }

  //! \brief get subdomain index
  int getSubdomainIndex (const Element& e) const
  {
    return getSubdomainIndex(getElementIndex(e));
  }

  //! \brief get subdomain index
  int getSubdomainIndex (const SubDomainElement& e) const
  {
    return getSubdomainIndex(getElementIndex(e));
  }

  //! \brief get subdomain index
  int getSubdomainIndex (const ElementIntersection& is) const
  {
    return getSubdomainIndex(*is.inside());
  }

  //! \brief get subdomain index
  int getSubdomainIndex (const SubDomainElementIntersection& is) const
  {
    return getSubdomainIndex(*is.inside());
  }

  //! \brief get group index
  int getGroupIndex (const Element& e) const
  {
    return getGroupIndex(getElementIndex(e));
  }

  //! \brief get subdomain index
  int getGroupIndex (const SubDomainElement& e) const
  {
    return getGroupIndex(getElementIndex(e));
  }

  int getGroupIndexByName(const std::string& groupName) const
  {
    for(int i=0; i<groupNames.size(); i++)
    {
      if(groupNames[i] == groupName)
        return i;
    }
    DUNE_THROW(Dune::Exception, "Could not find an element group for given name '" << groupName << "'!");
  }

  //! \brief get group name for a given element
  std::string getGroupName(const Element& e) const
  {
    return groupNames[getGroupIndex(e)];
  }

  //! \brief get group name for a given element
  std::string getGroupName(const SubDomainElement& se) const
  {
    return groupNames[getGroupIndex(se)];
  }

  //! \brief get group name for a given group index
  std::string getGroupName(const int groupIndex) const
  {
   return groupNames[groupIndex];
  }

  std::string getGroupInfo() const
  {
    std::stringstream info;
    info << "Groups: ";
    for(int i=0; i<groupNames.size(); i++)
    {
      info << "[" << i << "] = " << groupNames[i] << " ";
    }
    return info.str();
  }

  //! \brief get element index
  int getElementIndex (const Element& e) const
  {
    return elementMapper.map(e);
  }
  

  //! \brief get host element index of a subgrid entity
  int getElementIndex (const SubDomainElement& e) const
  {
    return getElementIndex(gv.grid().multiDomainEntity(e));
  }

  //! Get index of intersection by getting index of corresponding Codim-1 entity!
  //! This is following Jö's mail "Re_ [Dune] ID for an interface" from 16.05.2011, Dune mailing list
  int getIntersectionIndex(const ElementIntersection& is) const
  {
    return gv.indexSet().subIndex(*is.inside(), is.indexInInside(), 1);
  }

  int getIntersectionIndex(const SubDomainElementIntersection& is) const
  {
    return getIntersectionIndex(gv.grid().multiDomainEntity(is));
  }

  ElementIntersectionIterator getIntersection(const int index) const
  {
    return membraneInterfaceMap_self->at(index);
  }


  //! Get index of intersection by getting index of corresponding Codim-1 entity!
  //! This is following Jö's mail "Re_ [Dune] ID for an interface" from 16.05.2011, Dune mailing list
  typename GV::Grid::LocalIdSet::IdType getIntersectionID(const ElementIntersection& is) const
  {
    return gv.grid().localIdSet().subId(*is.inside(), is.indexInInside(), 1);
  }

  //! Get index of intersection by getting index of corresponding Codim-1 entity!
  //! This is following Jö's mail "Re_ [Dune] ID for an interface" from 16.05.2011, Dune mailing list
  typename SubGV::Grid::LocalIdSet::IdType getIntersectionID(const SubDomainElementIntersection& is) const
  {
    if(getSubDomainNumber(*is.inside()) == GridDomains::DOMAIN_MEMB)
      return membGV.grid().localIdSet().subId(*is.inside(), is.indexInInside(), 1);
    else
      return elecGV.grid().localIdSet().subId(*is.inside(), is.indexInInside(), 1);
  }


  int getSubDomainElementIndex(const Element& e) const
  {

    if(not isMembrane(e))
    {
      const SubDomainElement& se = *elecGV.grid().subDomainEntityPointer(e);
      assert(elecGV.contains(se));
      return getSubDomainElementIndex(se);
    } else {
      const SubDomainElement& se = *membGV.grid().subDomainEntityPointer(e);
      assert(membGV.contains(se));
      return getSubDomainElementIndex(se);
    }
  }

  int getSubDomainElementIndex(const SubDomainElement& se) const
  {
    assert(elecGV.contains(se) != membGV.contains(se));

    if(elecGV.contains(se))
      return elecElementMapper.map(se);
    else
      return membElementMapper.map(se);
  }

  int getSubDomainNumber(const Element& e) const
  {
    if(not isMembrane(e))
    {
      const SubDomainElement& se = *elecGV.grid().subDomainEntityPointer(e);
      assert(elecGV.contains(se));
      return getSubDomainNumber(se);
    } else {
      const SubDomainElement& se = *membGV.grid().subDomainEntityPointer(e);
      assert(membGV.contains(se));
      return getSubDomainNumber(se);
    }
  }

  //! get multidomaingrid subdomain number, one of { DOMAIN_ELEC = 0, DOMAIN_MEMB = 1}
  int getSubDomainNumber(const SubDomainElement& se) const
  {
    assert(elecGV.contains(se) != membGV.contains(se));

    if(elecGV.contains(se))
    {
      return elecGV.grid().domain();
    } else {
      return membGV.grid().domain();
    }
  }

  template<typename GF>
  int getDomain(GF& gf) const
  {
    // Specialization for IntersectionGridFunctions: Those are assumed to live on the membrane!
    if(GF::Traits::dimDomain == GF::Traits::GridViewType::dimension-1)
    {
      return GridDomains::DOMAIN_MEMB_INTERFACE;
    }

    return getDomain(gf.getGridView());
  }

  // Got multidomain gridview
  int getDomain(const GV& gv) const
  {
    return GridDomains::DOMAIN_ALL;
  }

  // Got subdomain gridview
  int getDomain(const SubGV& gv) const
  {
    return gv.grid().domain();
  }


  // ================================ BEGIN SubGrid <-> HostGrid stuff =====================================
  ElementIntersectionIterator getOppositeMembraneIntersection(const ElementIntersection& is) const
  {
    return membraneInterfaceMap_bidirectional->at(getIntersectionIndex(is));
  }

  ElementIntersectionIterator getPrimaryMembraneIntersection(const ElementIntersection& is) const
  {
    return membraneInterfaceMap_primary->at(getIntersectionIndex(is));
  }

  //! \brief Calculate the membrane surface area for a given membrane intersection by returning the
  //! surface area of the (only) intersection with the extracellular space
  T getMembraneSurfaceArea(const ElementIntersection& is) const
  {
    if(not isMembraneInterface(is))
    {
      DUNE_THROW(Dune::Exception, "HUHN-HAHN!!!!!!!");
    }

    // Get the 'primary' membrane intersection the channels are associated with; this is either
    // is itself or the opposite interface!
    const ElementIntersectionIterator iit_primary = membraneInterfaceMap_primary->at(getIntersectionIndex(is));

    // TODO Maybe remove this expensive assertion later on; compiling with -NDEBUG should do!
    assert(isIntersectionInMembraneInterfacesVector(*iit_primary));

    typedef typename ElementIntersectionIterator::Intersection::Geometry GEO_ORIG;
    typename Acme2CylGeometrySwitch::GeometrySwitch<GEO_ORIG>::type geo(iit_primary->geometry());

    return geo.volume();
  }

  //! \brief Calculated potential jump between given intersection and opposite membrane interface
  //! intersection using a discrete grid function for the calculated potential
  template<typename DGF_POT>
  void getMembranePotential(const ElementIntersection& is,
      const DGF_POT& dgfPot, typename DGF_POT::Traits::RangeType& membPot) const
  {
    if(not isMembraneInterface(is))
    {
      DUNE_THROW(Dune::Exception, "HUHN-HAHN!!!!!!!");
    }

    // Little tricky: Get the intersection which has the conventional orientation (inside element is membrane)
    int membInterfaceMapIndex = getIntersectionIndex(is);
    ElementIntersectionIterator iit = membraneInterfaceMap_self->at(membInterfaceMapIndex);

    typename DGF_POT::Traits::RangeType pot_Inside = 0.0;
    typename DGF_POT::Traits::RangeType pot_Outside = 0.0;
    membPot = 0.0;

    ElementIntersectionIterator iit_opposite = getOppositeMembraneIntersection(*iit);

    switch(getSubdomainIndex(*iit->outside()))
    {
      case CYTOSOL:
      {
        dgfPot.evaluate(*iit->inside(), iit->geometryInInside().center(), pot_Inside);
        dgfPot.evaluate(*iit_opposite->inside(), iit_opposite->geometryInInside().center(), pot_Outside);
        break;
      }
      case ES:
      {
        dgfPot.evaluate(*iit->inside(), iit->geometryInInside().center(), pot_Outside);
        dgfPot.evaluate(*iit_opposite->inside(), iit_opposite->geometryInInside().center(), pot_Inside);
        break;
      }
      default:
        DUNE_THROW(Dune::Exception, "Outside element of membrane interface is neither CYTOSOL nor ES!");
    }


    membPot = pot_Inside;
    membPot -= pot_Outside;
    //debug_verb << "Membrane potential: [" << pot_Inside << "  " << pot_Outside << "] --> " << membPot
    //    << " [" << convertTo_mV(membPot) << " mV]"<< std::endl;
  }

  //! \brief Calculated potential jump between given intersection and opposite membrane interface
  //! intersection using a discrete grid function for the calculated potential
  template<typename DGF_POT>
  void getMembranePotentialJump(const ElementIntersection& is, DGF_POT& dgfPot,
      typename DGF_POT::Traits::RangeType& potJump) const
  {
    if(not isMembraneInterface(is))
    {
      DUNE_THROW(Dune::Exception, "HUHN-HAHN!!!!!!!");
    }

    // Little tricky: Get the intersection which has the conventional orientation (inside element is membrane)
    int membInterfaceMapIndex = getIntersectionIndex(is);
    ElementIntersectionIterator iit = membraneInterfaceMap_self->at(membInterfaceMapIndex);

    potJump = 0.0;

    typename DGF_POT::Traits::RangeType pot1(0.0), pot2(0.0);
    typename DGF_POT::Traits::DomainType x1(0.0), x2(0.0);
    switch(getSubdomainIndex(*iit->outside()))
    {
      case CYTOSOL:
      case ES:
      {
        x1 = iit->outside()->geometry().global(iit->geometryInOutside().center());
        dgfPot.evaluate(*iit->outside(), iit->geometryInOutside().center(), pot1);

        ElementIntersectionIterator iit_opposite = getOppositeMembraneIntersection(*iit);

        assert(getIntersectionIndex(*iit_opposite) != getIntersectionIndex(*iit));

        x2 = iit_opposite->outside()->geometry().global(iit_opposite->geometryInOutside().center());
        dgfPot.evaluate(*iit_opposite->outside(), iit_opposite->geometryInOutside().center(), pot2);

        break;
      }
      default:
        DUNE_THROW(Dune::Exception, "Outside element of membrane interface is neither CYTOSOL nor ES!");
    }

    potJump = pot2;
    potJump -= pot1;

    // Multiply by sign of y difference (= y-axis direction)
    int ySign = Dune::sign(x1[1] - x2[1]);

    potJump *= ySign;

    // ySign positive: pot1 is top, pot2 is bottom
    if(ySign > 0)
    {
      debug_verb << "Potential bottom/top " << pot2 << "/" << pot1 << std::endl;
    } else {
      debug_verb << "Potential bottom/top " << pot1 << "/" << pot2 << std::endl;
    }
    debug_verb << "Potential jump (oriented): " << potJump
        << (potJump < 0 ? "↓" : "↑") << std::endl;
  }

  //! \brief New version of above method which is called from operator. Assume the membrane element is on the outside!
  template<typename DGF_POT>
  void getMembranePotentialJump(const ElementIntersection& is, DGF_POT& dgfPot,
      typename DGF_POT::Traits::RangeType& potJump, const typename DGF_POT::Traits::RangeType& uPotInside) const
  {
    if(not isMembraneInterface(is))
    {
      DUNE_THROW(Dune::Exception, "HUHN-HAHN!!!!!!!");
    }

    // Little tricky: Get the intersection which has the conventional orientation (inside element is membrane)
    int membInterfaceMapIndex = getIntersectionIndex(is);
    ElementIntersectionIterator iit = membraneInterfaceMap_self->at(membInterfaceMapIndex);

    potJump = 0.0;

    typename DGF_POT::Traits::RangeType pot1(0.0), pot2(0.0);
    typename DGF_POT::Traits::DomainType x1(0.0), x2(0.0);
    switch(getSubdomainIndex(*iit->outside()))
    {
      case CYTOSOL:
      case ES:
      {
        x1 = iit->outside()->geometry().global(iit->geometryInOutside().center());
        pot1 = uPotInside;

        // Use the previously obtained membrane intersection to get the opposite one from the stored map
        ElementIntersectionIterator iit_opposite = getOppositeMembraneIntersection(*iit);

        assert(getIntersectionIndex(*iit_opposite) != getIntersectionIndex(*iit));

        x2 = iit_opposite->outside()->geometry().global(iit_opposite->geometryInOutside().center());
        dgfPot.evaluate(*iit_opposite->outside(), iit_opposite->geometryInOutside().center(), pot2);

        break;
      }
      default:
        DUNE_THROW(Dune::Exception, "Outside element of membrane interface is neither CYTOSOL nor ES!");
    }

    potJump = pot2;
    potJump -= pot1;

    // Multiply by sign of y difference (= y-axis direction)
    int ySign = Dune::sign(x1[1] - x2[1]);

    potJump *= ySign;

    // ySign positive: pot1 is top, pot2 is bottom
    if(ySign > 0)
    {
      debug_verb << "Potential bottom/top " << pot2 << "/" << pot1 << std::endl;
    } else {
      debug_verb << "Potential bottom/top " << pot1 << "/" << pot2 << std::endl;
    }
    debug_verb << "Potential jump (oriented): " << potJump
        << (potJump < 0 ? "↓" : "↑") << std::endl;
  }


  //! deprecated
  //! \brief Get concentrations on both ends of the membrane
  template<typename DGF_CON>
  void getMembraneConcentrationJump(const ElementIntersection& is, DGF_CON& dgfCon,
    typename DGF_CON::Traits::RangeType& conCytosol, typename DGF_CON::Traits::RangeType& conExtracellular) const
  {
    if(not isMembraneInterface(is))
    {
      DUNE_THROW(Dune::Exception, "HUHN-HAHN!!!!!!!");
    }

    // Little tricky: Get the intersection which has the conventional orientation (inside element is membrane)
    int membInterfaceMapIndex = getIntersectionIndex(is);
    ElementIntersectionIterator iit = membraneInterfaceMap_self->at(membInterfaceMapIndex);

    ElementIntersectionIterator iit_opposite = getOppositeMembraneIntersection(*iit);

    assert(getIntersectionIndex(*iit_opposite) != getIntersectionIndex(*iit));

    switch(getSubdomainIndex(*iit->outside()))
    {
      // We need to evaluate concentrations on the outside of the intersection, as they are 0 on the membrane
      case CYTOSOL:
      {
        dgfCon.evaluate(*iit->outside(), iit->geometryInOutside().center(), conCytosol);
        dgfCon.evaluate(*iit_opposite->outside(), iit_opposite->geometryInOutside().center(), conExtracellular);
        break;
      }
      case ES:
      {
        dgfCon.evaluate(*iit->outside(), iit->geometryInOutside().center(), conExtracellular);
        dgfCon.evaluate(*iit_opposite->outside(), iit_opposite->geometryInOutside().center(), conCytosol);
        break;
      }
      default:
        DUNE_THROW(Dune::Exception, "Outside element of membrane interface is neither CYTOSOL nor ES!");
    }

    debug_verb << "Concentration inside/outside ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_verb << "[" << ION_NAMES[j] << "] = " << conCytosol[j] << "/" << conExtracellular[j] << "  ";
    }
    debug_verb << std::endl;

    debug_verb << "Concentration jump: ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_verb << "[" << ION_NAMES[j] << "] = " << (conCytosol[j] - conExtracellular[j]) << "  ";
    }
    debug_verb << std::endl;
  }

  //! deprecated
  //! \brief Get concentrations on both ends of the membrane
  template<typename DGF_CON>
  void getMembraneConcentrationJump(const ElementIntersection& is, DGF_CON& dgfCon,
    typename DGF_CON::Traits::RangeType& conCytosol, typename DGF_CON::Traits::RangeType& conExtracellular,
    const typename DGF_CON::Traits::RangeType& uConInside) const
  {
    if(not isMembraneInterface(is))
    {
      DUNE_THROW(Dune::Exception, "HUHN-HAHN!!!!!!!");
    }

    // Little tricky: Get the intersection which has the conventional orientation (inside element is membrane)
    int membInterfaceMapIndex = getIntersectionIndex(is);
    ElementIntersectionIterator iit = membraneInterfaceMap_self->at(membInterfaceMapIndex);

    ElementIntersectionIterator iit_opposite = getOppositeMembraneIntersection(*iit);

    assert(getIntersectionIndex(*iit_opposite) != getIntersectionIndex(*iit));

    switch(getSubdomainIndex(*iit->outside()))
    {
      // We need to evaluate concentrations on the outside of the intersection, as they are 0 on the membrane
      case CYTOSOL:
      {
        conCytosol = uConInside;
        dgfCon.evaluate(*iit_opposite->outside(), iit_opposite->geometryInOutside().center(), conExtracellular);
        break;
      }
      case ES:
      {
        conExtracellular = uConInside;
        dgfCon.evaluate(*iit_opposite->outside(), iit_opposite->geometryInOutside().center(), conCytosol);
        break;
      }
      default:
        DUNE_THROW(Dune::Exception, "Outside element of membrane interface is neither CYTOSOL nor ES!");
    }

    debug_verb << "Concentration CY/ES ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_verb << "[" << ION_NAMES[j] << "] = " << conCytosol[j] << "/" << conExtracellular[j] << "  ";
    }
    debug_verb << std::endl;

    debug_verb << "Concentration jump: ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_verb << "[" << ION_NAMES[j] << "] = " << (conCytosol[j] - conExtracellular[j]) << "  ";
    }
    debug_verb << std::endl;
  }


  //! \brief Calculated concentration jump across this membrane element. It calculates the concentration on both
  //! ends (i.e., neighboring electrolyte boundary) and saves the oriented difference in the parameter conJump.
  //! The sign of potJump gives the orientation of the jump in positive y-axis direction
  template<typename DGF_CON>
  void getMembraneConcentrationRatio(const ElementIntersection& is, DGF_CON& dgfCon,
    typename DGF_CON::Traits::RangeType& conRatio) const
  {
    if(not isMembraneInterface(is))
    {
      DUNE_THROW(Dune::Exception, "HUHN-HAHN!!!!!!!");
    }

    // Little tricky: Get the intersection which has the conventional orientation (inside element is membrane)
    int membInterfaceMapIndex = getIntersectionIndex(is);
    ElementIntersectionIterator iit = membraneInterfaceMap_self->at(membInterfaceMapIndex);

    int countMembraneInterfaces = 0;
    typename DGF_CON::Traits::RangeType con1(0.0), con2(0.0);
    typename DGF_CON::Traits::DomainType x1(0.0), x2(0.0);
    conRatio = 0.0;

    switch(getSubdomainIndex(*iit->outside()))
    {
      // We need to evaluate concentrations on the outside of the intersection, as they are 0 on the membrane
      case CYTOSOL:
      case ES:
      {
        x1 = iit->outside()->geometry().global(iit->geometryInOutside().center());
        dgfCon.evaluate(*iit->outside(), iit->geometryInOutside().center(), con1);

        // We are using the unidirectional map here: Assume this method was called from the lower
        // membrane interface stored in vector 'membraneInterfaces'!
        ElementIntersectionIterator iit_opposite = getOppositeMembraneIntersection(*iit);

        assert(getIntersectionIndex(*iit_opposite) != getIntersectionIndex(*iit));

        x2 = iit_opposite->outside()->geometry().global(iit_opposite->geometryInOutside().center());
        dgfCon.evaluate(*iit_opposite->outside(), iit_opposite->geometryInOutside().center(), con2);

        break;
      }
      default:
        DUNE_THROW(Dune::Exception, "Outside element of membrane interface is neither CYTOSOL nor ES!");
    }

    // Sign of y difference (= y-axis direction)
    int ySign = Dune::sign(x1[1] - x2[1]);
    debug_verb << "x1 = " << x1 << ", x2 = " << x2
      << " => sign = "<< ySign  << (ySign < 0 ? "↓" : "↑") << std::endl;

    // ySign positive: con1 is top, con2 is bottom
    if(ySign > 0)
    {
      conRatio = con2;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        conRatio[j] /= con1[j];
      }
    }
    // ySign is negative: con2 is top, con1 is bottom
    if(ySign < 0)
    {
      conRatio = con1;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        conRatio[j] /= con2[j];
      }
    }
    // conRatio is always con(bottom) / con(top), i.e. the ratio from the bottom membrane interface POV

    debug_verb << "Concentration bottom/top ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {// ySign positive: con1 is top, con2 is bottom
      if(ySign > 0)
      {
        debug_verb << "[" << ION_NAMES[j] << "] = " << con2[j] << "/" << con1[j] << "  ";
      } else {
        debug_verb << "[" << ION_NAMES[j] << "] = " << con1[j] << "/" << con2[j] << "  ";
      }
    }
    debug_verb << std::endl;

    debug_verb << "Concentration ratio (oriented): ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_verb << "[" << ION_NAMES[j] << "] = " << conRatio[j]
        << (conRatio[j] < 1.0 ? "↓" : "↑") << "  ";
    }
    debug_verb << std::endl;
  }

  //! \brief Calculated concentration jump across this membrane element. It calculates the concentration on both
  //! ends (i.e., neighboring electrolyte boundary) and saves the oriented difference in the parameter conJump.
  //! The sign of potJump gives the orientation of the jump in positive y-axis direction
  template<typename DGF_CON>
  void getMembraneConcentrationRatio(const ElementIntersection& is, DGF_CON& dgfCon,
    typename DGF_CON::Traits::RangeType& conRatio, const typename DGF_CON::Traits::RangeType& uConInside) const
  {
    if(not isMembraneInterface(is))
    {
      DUNE_THROW(Dune::Exception, "HUHN-HAHN!!!!!!!");
    }

    int countMembraneInterfaces = 0;
    typename DGF_CON::Traits::RangeType con1(0.0), con2(0.0);
    typename DGF_CON::Traits::DomainType x1(0.0), x2(0.0);
    conRatio = 0.0;

    // Little tricky: Get the intersection which has the conventional orientation (inside element is membrane)
    int membInterfaceMapIndex = getIntersectionIndex(is);
    ElementIntersectionIterator iit = membraneInterfaceMap_self->at(membInterfaceMapIndex);

    switch(getSubdomainIndex(*iit->outside()))
    {
      // We need to evaluate concentrations on the outside of the intersection, as they are 0 on the membrane
      case CYTOSOL:
      case ES:
      {
        x1 = iit->outside()->geometry().global(iit->geometryInOutside().center());
        con1 = uConInside;

        // We are using the unidirectional map here: Assume this method was called from the lower
        // membrane interface stored in vector 'membraneInterfaces'!
        ElementIntersectionIterator iit_opposite = getOppositeMembraneIntersection(*iit);

        assert(getIntersectionIndex(*iit_opposite) != getIntersectionIndex(*iit));

        x2 = iit_opposite->outside()->geometry().global(iit_opposite->geometryInOutside().center());
        dgfCon.evaluate(*iit_opposite->outside(), iit_opposite->geometryInOutside().center(), con2);

        break;
      }
      default:
        DUNE_THROW(Dune::Exception, "Outside element of membrane interface is neither CYTOSOL nor ES!");
    }

    // Sign of y difference (= y-axis direction)
    int ySign = Dune::sign(x1[1] - x2[1]);
    debug_verb << "x1 = " << x1 << ", x2 = " << x2
      << " => sign = "<< ySign  << (ySign < 0 ? "↓" : "↑") << std::endl;

    // ySign positive: con1 is top, con2 is bottom
    if(ySign > 0)
    {
      conRatio = con2;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        conRatio[j] /= con1[j];
      }
    }
    // ySign is negative: con2 is top, con1 is bottom
    if(ySign < 0)
    {
      conRatio = con1;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        conRatio[j] /= con2[j];
      }
    }
    // conRatio is always con(bottom) / con(top), i.e. the ratio from the bottom membrane interface POV

    debug_verb << "Concentration bottom/top ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      // ySign positive: con1 is top, con2 is bottom
      if(ySign > 0)
      {
        debug_verb << "[" << ION_NAMES[j] << "] = " << con2[j] << "/" << con1[j] << "  ";
      } else {
        debug_verb << "[" << ION_NAMES[j] << "] = " << con1[j] << "/" << con2[j] << "  ";
      }
    }
    debug_verb << std::endl;

    debug_verb << "Concentration ratio (oriented): ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_verb << "[" << ION_NAMES[j] << "] = " << conRatio[j]
        << (conRatio[j] < 1.0 ? "↓" : "↑") << "  ";
    }
    debug_verb << std::endl;
  }

  // ================================ END SubGrid <-> HostGrid stuff =======================================

  // ATTENTION: This methof only considers interior elements! This should not be used actually!
  const std::vector<T> serializeChannelStates() const
  DUNE_DEPRECATED_MSG("Use DynamicChannelSet::serializeChannelStates instead!")
  {
    int vSize = 0;

    // Determine total number of state values in this channelset
    for(int k=0; k<membrane.getChannelSet().size(); ++k)
    {
      vSize += membrane.getChannelSet().getChannel(k).numGatingParticles();
    }
    vSize *= nSubdomainElementsPerAxis[DOMAIN_MEMB][0];

    // Initialize state vector
    std::vector<T> states(vSize,0.0);

    debug_jochen << "Initialized channel state vector to size " << vSize << std::endl;

    //typename std::vector<T>::iterator it = states.begin();

    int index = 0;
    // Loop over membrane elements; only consider interior elements for output
    for(SubDomainElementIterator memb_it = membGV.template begin<0,Dune::PartitionIteratorType::Interior_Partition>();
          memb_it != membGV.template end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++memb_it)
    {
      assert(memb_it->partitionType() == Dune::PartitionType::InteriorEntity);

      // Loop over channels of this membrane element
      for(int k=0; k<membrane.getChannelSet().size(); ++k)
      {
        // Loop over this channel's gating particles
        for(int j=0; j<membrane.getChannelSet().getChannel(k).numGatingParticles(); j++)
        {
          // Get channel k, gating particle j, at membrane element i
          states[index] = membrane.getChannelSet().getGatingParticle(k,j,getElementIndex(*memb_it));
          index++;
        }
      }
    }

    assert(states.size() == vSize);
    return states;
  }

  //! \brief convert dimensionless potential to units of mV
  T convertTo_mV (const T& potential) const
  {
    if (params.mV_output())
    {
      return potential * con_k * electro.getTemperature() / con_e * 1.0e3;
    }
    else
    {
      return potential;
    }
  }

  //! \brief convert dimensionless potential to units of mV
  template<typename V>
  V convertTo_mV (const V& potential) const
  {
    V potential_mV(potential);
    if (params.mV_output())
    {
      for(int i=0; i<potential_mV.size(); i++)
      {
        potential_mV[i] = potential[i] * con_k * electro.getTemperature() / con_e * 1.0e3;
      }
    }
    return potential_mV;
  }

  //! \brief convert mV potential to dimensionless potential
  T convertFrom_mV (const T& potential) const
  {
    return potential * con_e / (1.0e3 * con_k * electro.getTemperature());
  }


  void getDebyeLength(std::vector<T>& debyeLength, const std::vector<T>& ionicStrength) const
  {
    const int dim = GV::dimension;
    const int intorder = 2;

    debyeLength.resize(2);

    T ionicStrength_in = ionicStrength[0];
    T ionicStrength_out = ionicStrength[1];
    debug_verb << "Intracellular ionic strength / N_A: " << ionicStrength_in << std::endl;
    debug_verb << "Extracellular ionic strength / N_A: " << ionicStrength_out << std::endl;

    T permittivity = electro.getPermittivity();

    debyeLength[0] = std::sqrt( con_eps0 * permittivity * con_k * electro.getTemperature() /
                ( 2.0 * con_e * con_e * electro.getStdCon() * ionicStrength_in) );
    debyeLength[1] = std::sqrt( con_eps0 * permittivity * con_k * electro.getTemperature() /
                ( 2.0 * con_e * con_e * electro.getStdCon() * ionicStrength_out) );

  }

  void getDebyeLength(std::vector<T>& debyeLength) const
  {
    std::vector<T> ionicStrength(2, 0.0);

    // Calculate ionic strength for intracellular/extracellular electrolytes
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      ionicStrength[0] += 0.5 * params.solution_in.get(ION_NAMES[j], 0.0) * getValence(j) * getValence(j);
      ionicStrength[1] += 0.5 * params.solution_ex.get(ION_NAMES[j], 0.0) * getValence(j) * getValence(j);
    }
    getDebyeLength(debyeLength, ionicStrength);
  }

  void gridInfo() const
  {
    const SubGV& elecGV = gv.grid().subDomain(0).leafView();
    const SubGV& membGV = gv.grid().subDomain(1).leafView();

    debug_verb << "================ Grid info: =========================" << std::endl;
    debug_verb << "Subgrid ELEC (" << elecGV.size(0) << " elements)" << std::endl;
    for(SubDomainElementIterator_All elec_it
        = elecGV.template begin<0>(); elec_it != elecGV.template end<0>(); ++elec_it)
    {
      int subElementIndex = elecGV.indexSet().index(*elec_it);
      debug_verb << "subGridElement[elec] #" << subElementIndex
         << " -> element #" << getElementIndex(*elec_it)
         << " [ " <<  elec_it->geometry().corner(0)
         << " - " << elec_it->geometry().center()
         << " - " <<  elec_it->geometry().corner(1) << "]"
         << ", partition: " << elec_it->partitionType()
         << std::endl;
    }
    debug_verb << "Subgrid MEMB (" << membGV.size(0) << " elements)" << std::endl;
    for(SubDomainElementIterator_All memb_it = membGV.template begin<0>();
            memb_it != membGV.template end<0>(); ++memb_it)
    {
      int subElementIndex = membGV.indexSet().index(*memb_it);
      debug_verb << "subGridElement[memb] #" << subElementIndex
         << " -> element #" << getElementIndex(*memb_it)
         << " [ " <<  memb_it->geometry().corner(0)
         << " - " << memb_it->geometry().center()
         << " - " <<  memb_it->geometry().corner(1) << "]"
         << ", partition: " << memb_it->partitionType()
         << std::endl;
    }

    debug_verb << std::endl;
    debug_verb << "subGV.contains() test:" << std::endl;
    // Test of member function "contains()"
    for(SubDomainElementIterator_All elec_it
        = elecGV.template begin<0>(); elec_it != elecGV.template end<0>(); ++elec_it)
    {
      assert(elecGV.contains(*elec_it));
    }
    for(SubDomainElementIterator_All memb_it
        = membGV.template begin<0>(); memb_it != membGV.template end<0>(); ++memb_it)
    {
      assert(! elecGV.contains(*memb_it));
    }
    for(SubDomainElementIterator_All elec_it
        = elecGV.template begin<0>(); elec_it != elecGV.template end<0>(); ++elec_it)
    {
      assert(! membGV.contains(*elec_it));
    }
    for(SubDomainElementIterator_All memb_it
        = membGV.template begin<0>(); memb_it != membGV.template end<0>(); ++memb_it)
    {
      assert(membGV.contains(*memb_it));
    }
    debug_verb << "Test successful!" << std::endl;
    debug_verb << "=====================================================" << std::endl;
  }

  void channelInfo() const
  {
    getChannelSet().info();
  }



  //! \brief debug output
  void info() const
  {
    debug_verb << "-------- Physics info -----------" << std::endl;
    debug_verb << "Element groups:" << std::endl;
    for(int i=0; i<groupNames.size(); ++i)
    {
      debug_verb << "G[" << i << "] = " << groupNames[i] << std::endl;
    }
    debug_verb << std::endl;
    const std::map<int,int>& membraneIndices = getChannelSet().getMembraneIndices();
    debug_verb << "Primary membrane interfaces:" << std::endl;
    for(MIterator mit = mBegin(); mit != mEnd(); ++mit)
    {
      // Inside (membrane) element belong to this membrane interface
      ElementPointer ep = mit->inside();
      int iIndex = getIntersectionIndex(*mit);
      int mIndex = membraneIndices.at(iIndex);
      int groupIndex = getGroupIndex(*ep);
      debug_verb << "I[" << iIndex << "] --> MI[" << mIndex << "]  -  G["
          << groupIndex << "] = " << getGroupName(groupIndex)
          << std::endl;
    }
    debug_verb << "---------------------------------" << std::endl;
  }


  //TODO Make this work generically
  int numberSubDomainInterfaces(const GV& gv) const
  {
    return params.useMembrane() ? 2 : 0;
  }

  //TODO Make this work generically
  int numberSubDomainInterfaces(const SubGV& gv) const
  {
    return 0;
  }

  Membrane<T>& getMembrane()
  {
    return membrane;
  }

  Electrolyte<T>& getElectrolyte()
  {
   return electro;
  }

  const ChannelSet& getChannelSet() const
  {
    return membrane.getChannelSet();
  }

  ChannelSet& getChannelSet()
  {
    return membrane.getChannelSet();
  }

  //! \brief get number of ion species
  int numOfSpecies () const
  {
    return electro.numOfSpecies();
  }

  //! \brief get permittivity
  T getPermittivity (const Element& e) const
  {
    if ( not isMembrane(e) )
    {
      //debug_jochen << "Element @ " << e.geometry().center() << ", permittivity = " << electro.getPermittivity() << std::endl;

      // electrolyte
      return electro.getPermittivity();
    } else {
      const SubDomainElementPointer sep = *membGV.grid().subDomainEntityPointer(e);
      int sdElemIndex = membElementMapper.map(*sep);
      //debug_jochen << "Element @ " << e.geometry().center() << ", permittivity = " << permittivities[sdElemIndex] << std::endl;

      // Membrane permittivities are stored in a separate vector
      return permittivities[sdElemIndex];
    }
  }

  //! \brief get membrane permittivity from membrane intersection
  T getMembranePermittivity (const ElementIntersection& is) const
  {
    assert(isMembraneInterface(is));

    if(isMembrane(*is.inside()))
    {
      return getPermittivity(*is.inside());
    } else {
      return getPermittivity(*is.outside());
    }
  }


  //! \brief get permittivity
  T getPermittivity (const SubDomainElement& se) const
  {
    if ( not isMembrane(se) )
    {
      // electrolyte
      return electro.getPermittivity();
    } else {
      int sdElemIndex = membElementMapper.map(se);
      // Membrane permittivities are stored in a separate vector
      return permittivities[sdElemIndex];
    }
  }

  //! \brief get group permittivity
  T getGroupPermittivity (const int groupIndex) const
  {
    // TODO Write this more general!
    if ( groupIndex < 2 )
    {
      // electrolyte
      return electro.getPermittivity();
    } else {
      // membrane: Try to find membrane group-specific permittivity in config file,
      // use default membrane permittivity as a fallback
      double perm = params.membrane.sub(getGroupName(groupIndex)).get("permittivity",
         -1);

      double default_dMemb = params.dMemb();
      double dMemb = params.membrane.sub(getGroupName(groupIndex)).get("d_memb", default_dMemb);

      perm *= (default_dMemb / dMemb);
      return perm;
    }
  }

  //! \brief get parameter tree (representing config file)
  const Acme2CylParameters& getParams() const
  {
    return params;
  }

  //! \brief get diffusion coefficient
  T getDiffCoeff (int ionSpecies, int elemIndex) const
  {
    int subdomainIndex = getSubdomainIndex(elemIndex);
    return getDiffCoeff(subdomainIndex, ionSpecies, elemIndex);
  }

  //! \brief get valence
  T getValence (const int ionSpecies) const
  {
    return electro.getValence( ionSpecies );
  }

  //! \brief get name of ion species
  std::string getIonName (const int ionSpecies) const
  {
    return electro.getIonName( ionSpecies );
  }

  //! \brief get squared length constant for poisson equation
  T getPoissonConstant() const
  {
    return electro.getPoissonConstant();
  }

  T getTimeStep() const
  {
    return dt;
  }

  void setTimeStep(T dt_)
  {
    dt = dt_;
  }

  //! \brief get the electrolyte
  const Electrolyte<T>& getElectrolyte() const
  {
    return electro;
  }

  //! \brief get the membrane
  const Membrane<T>& getMembrane() const
  {
    return membrane;
  }

  const T getTimeScale() const
  {
    return TIME_SCALE;
  }

  const T getLengthScale() const
  {
    return LENGTH_SCALE;
  }

  ElementSubdomainMapper& getElementSubdomainMapper()
  {
    return elemSubdomainMapper;
  }

  const GV& gridView() const
  {
    return gv;
  }

  //! This returns the number of interior elements on this process
  int nElements() const
  {
    return params.nElementsThisProcess();
  }

  //! This returns the number of interior elements on this process
  //! in the specified coordinate direction
  int nElements(int axis) const
  {
    return nGridElementsPerAxis[axis];
  }

  //! This returns the number of interior elements in the given gridview
  //! on this process
  int nElements(const GV& gv) const
  {
    return nElements();
  }

  //! This returns the number of interior elements in the given gridview
  //! on this process, specified coordinate direction
  int nElements(const GV& gv, const int axis) const
  {
    return nElements(axis);
  }

  //! This returns the number of interior elements in the given subdomain gridview
  //! on this process
  int nElements(const SubGV& subgv) const
  {
    return nSubdomainElementsPerAxis[getDomain(subgv)][0] * nSubdomainElementsPerAxis[getDomain(subgv)][1];
  }

  //! This returns the number of interior elements in the given subdomain gridview
  //! on this process in the specified coordinate direction
  int nElements(const SubGV& subgv, const int axis) const
  {
    return nSubdomainElementsPerAxis[getDomain(subgv)][axis];
  }

  //! This returns the number of interior elements in the given (sub)domain
  //! \param domainNumber must be one of DOMAIN_ELEC, DOMAIN_MEMB (see enum GridDomains in constants)
  int nSubdomainElements(const int subdomainNumber) const
  {
    return nSubdomainElementsPerAxis[subdomainNumber][0] * nSubdomainElementsPerAxis[subdomainNumber][1];
  }

  //! This returns the number of interior elements in the given (sub)domain
  //! \param domainNumber must be one of DOMAIN_ELEC, DOMAIN_MEMB (see enum GridDomains in constants)
  int nSubdomainElements(const int subdomainNumber, const int axis) const
  {
    return nSubdomainElementsPerAxis[subdomainNumber][axis];
  }

  template<typename GF>
  int nElementsForGF(const GF& gf) const
  {
    int domain = getDomain(gf);
    if(domain == GridDomains::DOMAIN_ALL)
      return nElements();
    // Special case membrane interface: Return total number of membrane elements divided by number of element layers
    if(domain == GridDomains::DOMAIN_MEMB_INTERFACE)
      return nSubdomainElements(GridDomains::DOMAIN_MEMB) / params.nMembraneElements();
    else
      return nSubdomainElements(domain);
  }

  template<typename GF>
  int nElementsForGF(const GF& gf, const int axis) const
  {
    int domain = getDomain(gf);
    if(domain == GridDomains::DOMAIN_ALL)
      return nElements(axis);
    // Special case membrane interface: Return total number of membrane elements divided by number of element layers
    // (y-direction only, in x-direction the number of elements and primary interfaces is identical)
    if(domain == GridDomains::DOMAIN_MEMB_INTERFACE)
    {
      if(axis == 1)
        return nSubdomainElements(GridDomains::DOMAIN_MEMB, axis) / params.nMembraneElements();
      else
        return nSubdomainElements(GridDomains::DOMAIN_MEMB, axis);
    }
    else
      return nSubdomainElements(domain,axis);
  }

  //! offset (in units "number of elements" that this process is shifted with respect to the macro grid,
  //! in the specified coordinate direction
  int nOffset(int axis) const
  {
    return nOffsetPerAxis[axis];
  }

  //! number of vertices on this process
  int nVertices() const
  {
    return gv.size(GV::dimension);
  }

  //! number of vertices in a given coordinate direction
  int nVertices(const int axis) const
  {
    assert(axis < GV::dimension);

    int nVerticesY = nElements(1)+1; // no partitioning in y-direction!
    int nVerticesX = nVertices() / nVerticesY;

    if(axis == 0) return nVerticesX;
    else return nVerticesY;
  }

  int nVerticesCytosol() const
  {
    // number of cytosol vertices = (#vertices in y direction) * (#vertices in x-direction on this process)
    return (params.getYOffsetMembrane() * nVertices(0));
  }

  //! \brief This is the number of inner membrane elements, without those on membrane interfaces; together with
  //! nVerticesCytosol() and nVerticesExtracellular() this gives a real partition of all grid vertices, i.e.
  //! nVerticesCytosol() + nVerticesMembrane() + nVerticesExtracellular() == gv.size(dimGrid)
  int nVerticesMembrane() const
  {
    // number of membrane vertices = (#inner membrane vertices in y-direction) * (#vertices in x-direction on this process)
    return (params.nMembraneElements()-1) * nVertices(0);
  }

  int nVerticesExtracellular() const
  {
    // number of extracellular vertices = (#remaining vertices in y-direction) * (#vertices in x-direction on this process)
    return ((nVertices(1) - params.getYOffsetMembrane() - (params.nMembraneElements()-1)) * nVertices(0));
  }

  const std::vector<int>& nElements_AllProcessors(int axis) const
  {
    assert(axis < 2);

    if(axis==0)
      return partition_x;
    else
      return partition_y;
  }

  const std::vector<int>& nOffset_AllProcessors(int axis) const
  {
    assert(axis < 2);

    if(axis==0)
      return offset_x;
    else
      return offset_y;
  }

  template<typename MultiGFS>
  void setupDOFPermutation(const MultiGFS& multigfs, std::vector<std::size_t>& permutation, std::vector<std::size_t>& inv_permutation) const
  {
    Dune::PDELab::EntityIndexCache<MultiGFS> eIndexCache(multigfs);
    typedef typename MultiGFS::Traits::GridViewType MultiGV;

    permutation.resize(multigfs.size());
    inv_permutation.resize(multigfs.size());

    std::size_t index = 0;
    typedef typename MultiGV::template Codim<MultiGV::dimension>::Iterator VertexIterator;
    for(VertexIterator vit = multigfs.gridView().template begin<GV::dimension>();
        vit != multigfs.gridView().template end<GV::dimension>(); ++vit)
    {
      eIndexCache.update(*vit);

      // Loop over all DOFs belonging to this vertex
      for(int i=0; i<eIndexCache.size(); i++)
      {
        permutation[index] = eIndexCache.containerIndex(i).back();
        inv_permutation[eIndexCache.containerIndex(i).back()] = index;
        index++;
      }
    }
  }

  template<typename DGF_POT>
  std::vector<FieldType> getBoundaryValues(const DGF_POT& dgfPot,
      typename DGF_POT::Traits::DomainFieldType& xReference,
      typename DGF_POT::Traits::DomainFieldType& yReference) const
  {
    typedef typename DGF_POT::Traits::DomainFieldType DF;
    typedef typename GV::template Codim<0>::Entity Entity;
    typedef typename GV::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GV::template Codim<0>::Iterator ElementLeafIterator;

    const GV& gv = dgfPot.getGridView();

    std::vector<FieldType> bValues = {0.0, 0.0, 0.0, 0.0};
    int count = 0;
    for (ElementLeafIterator eit = gv.template begin<0,Dune::PartitionIteratorType::Interior_Partition> ();
        eit != gv.template end<0,Dune::PartitionIteratorType::Interior_Partition> (); ++eit)
    {
      typename DGF_POT::Traits::ElementType& e = *eit;
      typename DGF_POT::Traits::DomainType firstCorner = e.geometry().corner(0);
      typename DGF_POT::Traits::DomainType lastCorner = e.geometry().corner(e.geometry().corners()-1);

      bool containsReferencePositionX = firstCorner[0]-1e-6 < xReference && lastCorner[0]+1e-6 > xReference;
      bool containsReferencePositionY = firstCorner[1]-1e-6 < yReference && lastCorner[1]+1e-6 > yReference;

      if(containsReferencePositionY)
      {
        typename DGF_POT::Traits::DomainType x(0.5);
        typename DGF_POT::Traits::RangeType y;

        x[0] = 0.0;
        // Left
        if(e.geometry().global(x)[0] - 1e-6 < params.xMin())
        {
          dgfPot.evaluate(e, x, y);
          debug_jochen << e.geometry().global(x) << "[LEFT]: " << y[0] << std::endl;
          bValues[1] = y[0];
          count++;
        }
        x[0] = 1.0;
        // Right
        if(e.geometry().global(x)[0] + 1e-6 > params.xMax())
        {
          dgfPot.evaluate(e, x, y);
          debug_jochen << e.geometry().global(x) << "[RIGHT]: " << y[0] << std::endl;
          bValues[2] = y[0];
          count++;
        }
      }
      if(containsReferencePositionX)
      {
        typename DGF_POT::Traits::DomainType x(0.5);
        typename DGF_POT::Traits::RangeType y;

        x[1] = 0.0;
        // Bottom
        if(e.geometry().global(x)[1] - 1e-6 < params.yMin())
        {
          dgfPot.evaluate(e, x, y);
          debug_jochen << e.geometry().global(x) << "[BOTTOM]: " << y[0] << std::endl;
          bValues[0] = y[0];
          count++;
        }
        x[1] = 1.0;
        // Top
        if(e.geometry().global(x)[1] + 1e-6 > params.yMax())
        {
          dgfPot.evaluate(e, x, y);
          debug_jochen << e.geometry().global(x) << "[TOP]: " << y[0] << std::endl;
          bValues[3] = y[0];
          count++;
        }
      }
    }
    if(count < 4)
      DUNE_THROW(Dune::Exception, "Error finding boundary values!");

    return bValues;
  }

  T getBulkConcentrationByElement(const int ionSpecies, const int elemIndex) const
  {
    return getBulkConcentration(getSubdomainIndex(elemIndex), ionSpecies);
  }

  T getIntracellularBulkConcentration(const int ionSpecies) const
  {
    return getBulkConcentration(CYTOSOL, ionSpecies);
  }

  T getExtracellularBulkConcentration(const int ionSpecies) const
  {
    return getBulkConcentration(ES, ionSpecies);
  }

  T getNodeVolume(const CoordType& pos) const
  {
    assert(nodeVolumesMap.count(pos) > 0);

    return nodeVolumesMap.at(pos);
  }

  MIterator mBegin() const
  {
    //return mbegin;
    //return membraneInterfaces.begin();
    return MIterator(*membraneInterfaces, 0);
  }

  MIterator mEnd() const
  {
    //return mend;
    //return membraneInterfaces.end();
    return MIterator(*membraneInterfaces, membraneInterfaces->size());
  }

  MIterator mInteriorBegin() const
  {
    //return mbegin;
    //return membraneInterfaces.begin();
    return MIterator(*membraneInterfaces, 0, Dune::PartitionIteratorType::Interior_Partition);
  }

  MIterator mInteriorEnd() const
  {
    //return mend;
    //return membraneInterfaces.end();
    return MIterator(*membraneInterfaces, membraneInterfaces->size(), Dune::PartitionIteratorType::Interior_Partition);
  }

  void printIntersectionInfo() const
  {
    debug_info << "Finished setting up membrane interfaces map, # of associations: "
        << membraneInterfaceMap_bidirectional->size() << std::endl;

    debug_jochen << "Complete bidirectional map: " << std::endl;
    for(typename MembraneInterfaceMap::const_iterator it = membraneInterfaceMap_bidirectional->begin();
        it != membraneInterfaceMap_bidirectional->end(); ++it)
    {
      debug_jochen << it->first << " - " << getIntersectionIndex(*it->second)
          << " @" << (*it->second).geometry().center()
          << " [inside @" << (*it->second).inside()->geometry().center() << "]" << std::endl;
    }

    debug_jochen << "Complete primary map: " << std::endl;
    for(typename MembraneInterfaceMap::const_iterator it = membraneInterfaceMap_primary->begin();
        it != membraneInterfaceMap_primary->end(); ++it)
    {
      debug_jochen << it->first << " - " << getIntersectionIndex(*it->second)
            << " @" << (*it->second).geometry().center()
            << " [inside @" << (*it->second).inside()->geometry().center() << "]" << std::endl;
    }
  }

  void testMembraneIntersectionMapping()
  {
    typedef typename ElementIntersectionIterator::Intersection::Geometry::GlobalCoordinate Coord;

    for(ElementIterator_All eit = gv.template begin<0>();
                eit != gv.template end<0>(); ++eit)
    {
      // Convention: Membrane interfaces have the membrane element on their 'inside'
      if(isMembrane(*eit))
      {
        for(ElementIntersectionIterator iit = gv.ibegin(*eit); iit != gv.iend(*eit); ++iit)
        {
          if(isMembraneInterface(*iit))
          {
            debug_jochen << "[" << getGroupName(*eit) << "] inside element" << std::endl;
            debug_jochen << "Trying to get opposite intersection for membrane intersection #" << getIntersectionIndex(*iit)
                << " @" << iit->geometry().center() << "..." << std::endl;
            ElementIntersectionIterator iit_opp = getOppositeMembraneIntersection(*iit);

            debug_jochen << "[0] That worked, opposite intersection #" << getIntersectionIndex(*iit_opp) << " is @"
                << iit_opp->geometry().center() << std::endl;
          }
        }
      } else { // else case is just a TEST, remove
        for(ElementIntersectionIterator iit = gv.ibegin(*eit); iit != gv.iend(*eit); ++iit)
        {
          if(isMembraneInterface(*iit))
          {
            debug_jochen << "[" << getGroupName(*eit) << "] inside element" << std::endl;
            try {
              debug_jochen << "Trying to get opposite intersection for membrane intersection #" << getIntersectionIndex(*iit)
                  << " @" << iit->geometry().center() << "..." << std::endl;
              ElementIntersectionIterator iit_opp = getOppositeMembraneIntersection(*iit);

              debug_jochen << "[1] That worked, opposite intersection #" << getIntersectionIndex(*iit_opp) << " is @"
                  << iit_opp->geometry().center() << std::endl;

              ElementIntersectionIterator iit_conv = getIntersection(getIntersectionIndex(*iit));
              debug_jochen << "[1] Conventional intersection #" << getIntersectionIndex(*iit_opp) << " found @"
                  << iit_conv->geometry().center() << std::endl;
            } catch(...) {
              debug_jochen << "Did not work! Trying different way..." << std::endl;

              // Mean hack: Get the index of the desired intersection (which is spatially the same one we already have,
              // but with the roles of inside and outside swapped) by assuming that in 2D and with horizontally oriented
              // membrane, the indexInInside is either 2 or 3. So just switch
              int myIndexInInside = iit->indexInInside();

              debug_jochen << "Index in inside: " << myIndexInInside << std::endl;

              if(myIndexInInside != 2 && myIndexInInside != 3)
                DUNE_THROW(Dune::Exception, "Expected indexInInside to be either 2 or 3!");

              int indexInInside = (myIndexInInside == 2 ? 3 : 2);
              int membInterfaceMapIndex = gv.indexSet().subIndex(*iit->outside(), indexInInside, 1);

              debug_jochen << "Calculated map index: " << membInterfaceMapIndex << std::endl;

              ElementIntersectionIterator iit_opp = membraneInterfaces->at(membInterfaceMapIndex);

              debug_jochen << "[2] That worked, opposite intersection #" << membInterfaceMapIndex << " is @"
                  << iit_opp->geometry().center() << std::endl;

              ElementIntersectionIterator iit_conv = getIntersection(getIntersectionIndex(*iit));
              debug_jochen << "[2] Conventional intersection #" << getIntersectionIndex(*iit_opp) << " found @"
                  << iit_conv->geometry().center() << std::endl;
            }

          }
        }
      }
    }
  }

  template<typename DGF_POT, typename DGF_CON>
  void testMembraneFunctions(const DGF_POT& dgfPot, const DGF_CON& dgfCon)
  {
    typedef typename ElementIntersectionIterator::Intersection::Geometry::GlobalCoordinate Coord;

    typename DGF_POT::Traits::RangeType potTestCY(-3);
    typename DGF_POT::Traits::RangeType potTestES(5);

    typename DGF_CON::Traits::RangeType conTestCY(100.);
    typename DGF_CON::Traits::RangeType conTestES(200.);

    debug_info << "========================= Physics::testMembraneFunctions ========================" << std::endl;

    for(ElementIterator_All eit = gv.template begin<0>();
                eit != gv.template end<0>(); ++eit)
    {
      // Convention: Membrane interfaces have the membrane element on their 'inside'
      if(isMembrane(*eit))
      {
        for(ElementIntersectionIterator iit = gv.ibegin(*eit); iit != gv.iend(*eit); ++iit)
        {
          if(isMembraneInterface(*iit))
          {
            try {
              typename DGF_POT::Traits::RangeType pot;
              typename DGF_CON::Traits::RangeType con;
              typename DGF_CON::Traits::RangeType con2;
              debug_info << "[" << getGroupName(*eit) << "] inside element, intersection index #"
                  << getIntersectionIndex(*iit) << std::endl;

              // getMembranePotential
              getMembranePotential(*iit, dgfPot, pot);
              debug_info << "  OLD getMembranePotential:" << pot << std::endl;

              // getMembranePotentialJump
              getMembranePotentialJump(*iit, dgfPot, pot);
              debug_info << "  OLD getMembranePotentialJump:" << pot << std::endl;

              // getMembraneConcentrationRatio
              getMembraneConcentrationRatio(*iit, dgfCon, con);
              debug_info << "  OLD getMembraneConcentrationRatio:" << con << std::endl;

              // getMembraneConcentrationJump
              getMembraneConcentrationJump(*iit, dgfCon, con, con2);
              debug_info << "  OLD getMembraneConcentrationJump:" << con << " / " << con2 << std::endl;
            } catch (std::exception& e) {
              debug_warn << "  Exception @" << iit->geometry().center() << std::endl;
              throw e;
            }
          }
        }
      } else { // else case is just a TEST, remove
        for(ElementIntersectionIterator iit = gv.ibegin(*eit); iit != gv.iend(*eit); ++iit)
        {
          if(isMembraneInterface(*iit))
          {
            try {
              typename DGF_POT::Traits::RangeType pot;
              typename DGF_CON::Traits::RangeType con;
              typename DGF_CON::Traits::RangeType con2;
              debug_info << "[" << getGroupName(*eit) << "] inside element, intersection index #"
                  << getIntersectionIndex(*iit) << std::endl;

              // getMembranePotential
              getMembranePotential(*iit, dgfPot, pot);
              debug_info << "  OLD getMembranePotential:" << pot << std::endl;
  //            if(getGroupName(*eit) == "solution_ex")
  //            {
  //              getMembranePotential(*iit, dgfPot, pot, potTestES);
  //            } else if(getGroupName(*eit) == "solution_in")
  //            {
  //              getMembranePotential(*iit, dgfPot, pot, potTestCY);
  //            }
  //            debug_info << "  NEW getMembranePotential:" << pot << std::endl;

              // getMembranePotentialJump
              getMembranePotentialJump(*iit, dgfPot, pot);
              debug_info << "  OLD getMembranePotentialJump:" << pot << std::endl;
              if(getGroupName(*eit) == "solution_ex")
              {
                getMembranePotentialJump(*iit, dgfPot, pot, potTestES);
              } else if(getGroupName(*eit) == "solution_in")
              {
                getMembranePotentialJump(*iit, dgfPot, pot, potTestCY);
              }
              debug_info << "  NEW getMembranePotentialJump:" << pot << std::endl;

              // getMembraneConcentrationRatio
              getMembraneConcentrationRatio(*iit, dgfCon, con);
              debug_info << "  OLD getMembraneConcentrationRatio:" << con << std::endl;
              if(getGroupName(*eit) == "solution_ex")
              {
                getMembraneConcentrationRatio(*iit, dgfCon, con, conTestES);
              } else if(getGroupName(*eit) == "solution_in")
              {
                getMembraneConcentrationRatio(*iit, dgfCon, con, conTestCY);
              }
              debug_info << "  NEW getMembraneConcentrationJump:" << con << std::endl;

              // getMembraneConcentrationJump
              getMembraneConcentrationJump(*iit, dgfCon, con, con2);
              debug_info << "  OLD getMembraneConcentrationJump:" << con << " / " << con2 << std::endl;
              if(getGroupName(*eit) == "solution_ex")
              {
                getMembraneConcentrationJump(*iit, dgfCon, con, con2, conTestES);
              } else if(getGroupName(*eit) == "solution_in")
              {
                getMembraneConcentrationJump(*iit, dgfCon, con, con2, conTestCY);
              }
              debug_info << "  NEW getMembraneConcentrationJump:" << con << " / " << con2 << std::endl;

            } catch (std::exception& e) {
              debug_warn << "  Exception @" << iit->geometry().center() << std::endl;
              throw e;
            }
          }
        }
      }
    }
  }


private:

  const ElementIntersectionIterator& getNextMembraneInterface(const ElementIntersectionIterator& iit_current) const
  {
    //debug_jochen << "getNext... received intersection @" << iit_current->geometry().center() << std::endl;
    //debug_jochen << " with inside @" << iit_current->inside()->geometry().center() << std::endl;
    //debug_jochen << " and outside @" << iit_current->outside()->geometry().center() << std::endl;

    ElementPointer ep = iit_current->outside();
    if(isMembrane(*ep))
    {
      //debug_jochen << "Trying to find membrane interface for membrane element @" << ep->geometry().center()
      //    << std::endl;

      ElementIntersectionIterator iit_next = gv.ibegin(*ep);
      for(ElementIntersectionIterator iit = gv.ibegin(*ep); iit != gv.iend(*ep); ++iit)
      {
        //debug_jochen << "-checking if intersection " << iit->geometry().center() << " is a candidate..." << std::endl;
        // End of recursion: Membrane interface was found!
        if(isMembraneInterface(*iit))
        {
          //debug_jochen << "-found membrane intersection " << iit->geometry().center() << "!" << std::endl;
          //debug_jochen << "-returning membrane intersection "
          //    << membraneInterfaceMap->at(getIntersectionIndex(*iit))->geometry().center() << "!" << std::endl;

          // Now call membrane interface map in order to always return the primary intersection!
          return membraneInterfaceMap_primary->at(getIntersectionIndex(*iit));
        }
        Acme2CylGeometry<typename ElementIntersectionIterator::Intersection::Geometry> cylGeo(iit->geometry());
        // Save candidate for next recursion call; search direction orthogonal to membrane orientation!
        if(cylGeo.isOrientedAxially())
        {
          iit_next = iit;
          //debug_jochen << " -new candidate for recursive call: "
          //  << iit_next->outside()->geometry().center()
          //  << std::endl;
        }
      }
      // Recursively call getNextMembraneInterface until all membrane element layers have been traversed;
      return getNextMembraneInterface(iit_next);
    }
    DUNE_THROW(Dune::Exception, "Given element is not a membrane element!");
  }

  // return bulk concentration from config file
  T getBulkConcentration(const int ionSpecies, const int subdomainIndex) const
  {
    switch(subdomainIndex)
    {
      case CYTOSOL:
        return params.solution_in.get(ION_NAMES[ionSpecies], -1.);
      case ES:
        return params.solution_ex.get(ION_NAMES[ionSpecies], -1.);
      case MEMBRANE:
        return 0.0;
      default:
        DUNE_THROW(Dune::Exception, "Could not find specified subdomain index '" << subdomainIndex
            << "'!");
    }
  }

  /*! \brief This method tests if an element with the given
   * subdomainIndex is a membrane element or not.
   *
   * This method should always be used instead of testing for
   * subdomainIndex <= 1 or something alike
   * since this specific index might change and then only has to
   * be modified here instead of everywhere else.
   */
  bool isMembrane(int subdomainIndex) const
  {
    if (subdomainIndex == MEMBRANE)
    {
      return true;
    } else {
      return false;
    }
  }

  //! \brief get subdomain index, one of { CYTOSOL = 0, ES = 1, MEMBRANE = 2} (see enum in constants.hh)
  int getSubdomainIndex (const int elementIndex) const
  {
    return elemSubdomainMapper.map( elementIndex );
  }

  //! \brief get subdomain index
  int getGroupIndex (const int elementIndex) const
  {
    return elemGroupMapper.map( elementIndex );
  }


  //! \brief get diffusion coefficient
  T getDiffCoeff (int subdomainIndex, int ionSpecies, int elemIndex) const
  {
    // We don't calculate concentrations inside the membrane!
    //assert(isMembrane(subdomainIndex) == false);
    if(isMembrane(subdomainIndex)) return 0.0;

    return electro.getDiffConst( ionSpecies );
  }

  void setupMembraneInterfaceMap()
  {
    typedef typename ElementIntersectionIterator::Intersection::Geometry::GlobalCoordinate Coord;

    // Collect all membrane interfaces in a single vector; important: loop over the whole domain, we need
    // multidomain intersections, not subdomain intersections!
    std::vector<ElementIntersectionIterator> membraneInterfaces_temp;
    for(ElementIterator_All eit = gv.template begin<0>();
                eit != gv.template end<0>(); ++eit)
    {
      // Convention: Membrane interfaces have the membrane element on their 'inside'
      if(isMembrane(*eit))
      {
        for(ElementIntersectionIterator iit = gv.ibegin(*eit); iit != gv.iend(*eit); ++iit)
        {
          if(isMembraneInterface(*iit))
          {
            membraneInterfaces_temp.push_back(iit);

//            debug_jochen << gv.indexSet().subIndex(*iit->inside(), iit->indexInInside(), 1) << std::endl;

            // Save all membrane-extracellular interfaces (only one side of the membrane) in vector
            if(!iit->boundary() && getSubdomainIndex(*iit->outside()) == ES)
              membraneInterfaces->push_back(iit);
          }
        }
      }
    }
    debug_info << "Found " << membraneInterfaces_temp.size() << " membrane interfaces ("
        << membraneInterfaces->size() << " on membrane-extracellular boundary)."<< std::endl;

    // Loop over found membrane interfaces and find opposite interfaces
    for(int i=0; i<membraneInterfaces->size(); i++)
    {
      ElementIntersectionIterator iit = (*membraneInterfaces)[i];
      //IdType id = getIntersectionID(*iit);
      int index = getIntersectionIndex(*iit);

//      debug_jochen << "This is membrane interface #" << index << " @" << iit->geometry().center()
//          << " with associated inside element @" << iit->inside()->geometry().center() << std::endl;

      Coord center_this(iit->geometry().center());

      int j=0;
      bool found = false;
      // Search opposite interface in vector of all membrane interfaces
      for(j=0; j<membraneInterfaces_temp.size(); j++)
      {
        Coord center_other(membraneInterfaces_temp[j]->geometry().center());

        // Compare position; x coordinates should match exactly, y coordinates should be dMemb apart
        bool xEqual = (std::abs(center_this[0]-center_other[0]) < 1e-6);
        bool yEqual = (std::abs(center_this[1]-center_other[1]) < 1e-6);
        bool yMatch = (std::abs(center_this[1]-center_other[1])-1e-6 < params.dMemb());

//        debug_jochen << "  Comparing membrane interfaces @" << center_this << " and " << center_other
//            << ", xEqual:" << xEqual << ", yEqual:" << yEqual << ", yMatch:" << yMatch << std::endl;

        if(xEqual && !yEqual && yMatch)
        {
          found = true;
          break;
        }
      }
      if(! found)
        DUNE_THROW(Dune::Exception, "Could not find a matching opposite membrane interface!");

      ElementIntersectionIterator iit_opposite = membraneInterfaces_temp[j];
      //IdType id_opposite = getIntersectionID(*iit_opposite);
      int index_opposite = getIntersectionIndex(*iit_opposite);

//      debug_verb << "Associating with membrane interface #" << index_opposite << " @"
//          << iit_opposite->geometry().center()
//          << " with associated inside element @" << iit_opposite->inside()->geometry().center()
//          << std::endl << std::endl;

      // There are two versions of the membrane interfaces map: The bidirectional variant saves
      // associations in both directions, i.e. primary->opposite, opposite->primary;
      // the other one always maps to the primary membrane interface, i.e. primary->primary, opposite->primary;
      // this is useful because channels are associated with the primary one, so in practice one
      // always needs to get the index of the primary membrane index, which is also present in
      // vector 'membraneInterfaces'
      std::pair<int,ElementIntersectionIterator> primaryToOpposite(index,iit_opposite);
      std::pair<int,ElementIntersectionIterator> oppositeToPrimary(index_opposite,iit);
      membraneInterfaceMap_bidirectional->insert(primaryToOpposite); // primary->opposite
      membraneInterfaceMap_bidirectional->insert(oppositeToPrimary); // opposite->primary

      std::pair<int,ElementIntersectionIterator> primaryToPrimary(index,iit); //identity primary->primary
      membraneInterfaceMap_primary->insert(primaryToPrimary); // primary->primary
      membraneInterfaceMap_primary->insert(oppositeToPrimary); // opposite->primary

      std::pair<int,ElementIntersectionIterator> selfPrimary(index,iit);
      std::pair<int,ElementIntersectionIterator> selfOpposite(index_opposite,iit_opposite);
      membraneInterfaceMap_self->insert(selfPrimary); // primary->primary
      membraneInterfaceMap_self->insert(selfOpposite); // opposite->opposite
    }

//    debug_jochen << "sizeof(membraneInterfaces): "
//        << (sizeof(*membraneInterfaces) + sizeof(ElementIntersectionIterator) * membraneInterfaces->capacity())
//        << std::endl;
//    debug_jochen << "sizeof(membraneInterfaceMap): "
//            << (sizeof(*membraneInterfaceMap) + sizeof(ElementIntersectionIterator) * membraneInterfaceMap->size())
//            << std::endl;
//    debug_jochen << "sizeof(membraneInterfaceMap_bidirectional): "
//            << (sizeof(*membraneInterfaceMap_bidirectional) + sizeof(ElementIntersectionIterator) * membraneInterfaceMap_bidirectional->size())
//            << std::endl;

    debug_info << "-- membraneInterfaceMap.size(): " << membraneInterfaceMap_primary->size() << std::endl;
    debug_info << "-- membraneInterfaceMap_bidirectional.size(): " << membraneInterfaceMap_bidirectional->size()
      << std::endl;
    debug_info << "-- membraneInterfaceMap_self.size(): " << membraneInterfaceMap_self->size() << std::endl;

    //printIntersectionInfo();
    //testIterator();
  }


  void setupChannels()
  {
    ChannelSet& channels = membrane.getChannelSet();

    // Loop over membrane interfaces
    for(MIterator mit = mBegin(); mit != mEnd(); ++mit)
    {
      int iIndex = getIntersectionIndex(*mit);
      int subdomainIndex = getSubdomainIndex(*mit);

      switch(subdomainIndex)
      {
        case CYTOSOL:
        case ES:
          DUNE_THROW(Dune::Exception, "Found a non-membrane entity inside the membrane subdomain grid!");
          break;
        case MEMBRANE:
        {
          channels.addMembraneElement(iIndex);
          break;
        }
        default:
          DUNE_THROW(Dune::Exception, "Element has an unknown subdomain index!");
      }
    }
    channels.resize();
    debug_info << "Added " << channels.getMembraneIndices().size() << " membrane elements to channelset!"
        << std::endl;

    channels.init(*this);
    debug_info << "Channels successfully initialized! " << std::endl;


  }

  //! \brief Print all membrane intersection saved in vector
  void testIterator()
  {
    // Now loop that shit
    debug_verb << "All intersection interfaces on this processor:" << std::endl;
    for(MIterator mit = mBegin(); mit != mEnd(); mit++)
    {
      debug_jochen << "-- " << mit->geometry().center() << std::endl;
    }
    debug_verb << "Interior intersection interfaces on this processor:" << std::endl;
    for(MIterator mit = mInteriorBegin(); mit != mInteriorEnd(); mit++)
    {
      debug_jochen << "-- " << mit->geometry().center() << std::endl;
    }
  }

  //! This method checks how many elements are contained in the complete grid
  //! with respect to a certain coordinate axis. For a 2D grid, nElements(0)
  //! would return the number of elements in x-direction and nElements(1)
  //! the number of elements in y-direction; works for 2D only!
  void calcNElements()
  {
    T yCurrent(params.yMin());
    int elemCount = 0;

    if(params.nElementsThisProcess() == 1)
    {
      nGridElementsPerAxis[0] = 1;
      nGridElementsPerAxis[1] = 1;
      return;
    }

    // Loop over elements
    for (ElementIterator eit = gv.template begin<0,Dune::PartitionIteratorType::Interior_Partition> ();
        eit != gv.template end<0,Dune::PartitionIteratorType::Interior_Partition> (); ++eit)
    {
      elemCount++;

      // Calculate the offset of the first element in this processor's gridview
      // (assume vertical grid partitioning and element ordering in ascending x-direction)
      if(elemCount == 1)
      {
        typename ElementIterator::Entity::Geometry::GlobalCoordinate center = eit->geometry().center();

        // Find x coordinate in global vector x
        const std::vector<T>& x_global = params.X();

        int offset_x = 0;
        for(; offset_x<x_global.size(); offset_x++)
        {
          if(x_global[offset_x] > center[0]) break;
        }
        assert(offset_x<x_global.size());

        nOffsetPerAxis[0] = offset_x-1;
        nOffsetPerAxis[1] = 0;

        debug_jochen << "p" << gridView().comm().rank() << " offset x: " << nOffsetPerAxis[0]
          << ", y: " << nOffsetPerAxis[1] << std::endl;

        assert(nOffsetPerAxis[0] < params.X().size());
        assert(nOffsetPerAxis[1] < params.Y().size());
      }

      //debug_jochen << "yCurrent = " << yCurrent << ", corner(0) = " << eit->geometry().corner(0)[1] << std::endl;
      if(std::abs(yCurrent - eit->geometry().corner(0)[1]) > 1e-8)
      {
        break;
      }
    }
    nGridElementsPerAxis[0] = elemCount - 1;
    nGridElementsPerAxis[1] = params.nElementsThisProcess() / nGridElementsPerAxis[0];
    debug_jochen << "nGridElementsPerAxis[0]: " << nGridElementsPerAxis[0] << std::endl;
    debug_jochen << "nGridElementsPerAxis[1]: " << nGridElementsPerAxis[1] << std::endl;
    debug_jochen << "params.nElementsThisProcess() " << params.nElementsThisProcess() << std::endl;
    assert(nGridElementsPerAxis[0] * nGridElementsPerAxis[1] == params.nElementsThisProcess());

    // Membrane elements
    int nMembraneLayers = 0;
    if(params.useMembrane())
    {
      nMembraneLayers = params.nMembraneElements() * params.nMembranes();

      nSubdomainElementsPerAxis[GridDomains::DOMAIN_MEMB].push_back(nGridElementsPerAxis[0]);
      nSubdomainElementsPerAxis[GridDomains::DOMAIN_MEMB].push_back(nMembraneLayers);
    } else {
      nSubdomainElementsPerAxis[GridDomains::DOMAIN_MEMB].push_back(0);
      nSubdomainElementsPerAxis[GridDomains::DOMAIN_MEMB].push_back(0);
    }

    // Elec elements
    nSubdomainElementsPerAxis[GridDomains::DOMAIN_ELEC].push_back(nGridElementsPerAxis[0]);
    nSubdomainElementsPerAxis[GridDomains::DOMAIN_ELEC].push_back(nGridElementsPerAxis[1]-nMembraneLayers);

    debug_jochen << "MEMB #elems x :" << nSubdomainElementsPerAxis[GridDomains::DOMAIN_MEMB][0] << std::endl;
    debug_jochen << "MEMB #elems y :" << nSubdomainElementsPerAxis[GridDomains::DOMAIN_MEMB][1] << std::endl;
    debug_jochen << "ELEC #elems x :" << nSubdomainElementsPerAxis[GridDomains::DOMAIN_ELEC][0] << std::endl;
    debug_jochen << "ELEC #elems y :" << nSubdomainElementsPerAxis[GridDomains::DOMAIN_ELEC][1] << std::endl;
  }

  void calcNodeVolumes()
  {
    debug_jochen << "Physics::calcNodeVolumes" << std::endl;
    double reference_volume = -1.0;

    // This is the threshold y-coordinate, from which on the volume scaling is started. Nodes with a y-coordinate
    // below yStart are not scaled!

    std::vector<double> yMembDefault = params.yMemb();
    // Do not use vector returned by params.yMemb(), this might be a 1-element vector containing a zero if we don't
    // use a membrane! Instead read the actual value from config file, this makes more sense in this case.
    std::vector<double> yMemb = params.general.get("y_memb", yMembDefault);
    double cellWidth = (params.nMembranes() > 1 ?
        yMemb[yMemb.size()-1]-yMemb[yMemb.size()-2] : yMemb[0]);

    // Same for dMemb: Use actual value from config file
    double dMembDefault = params.dMemb();
    double dMemb = params.general.get("d_memb", dMembDefault);

    double yStart = params.general.get("volumeScalingYThreshold",yMemb.back() + dMemb + cellWidth);

    double minAnisotropy = std::numeric_limits<double>::max();
    double maxAnisotropy = 0.0;

    // Iterate over elements and assign minimum volumes to nodes
    // Loop over All_Partition, overlap entities need these values as well!
    for (ElementIterator_All eit = gv.template begin<0>(); eit != gv.template end<0>(); ++eit)
    {
      typedef typename ElementIterator::Entity::Geometry GEO_ORIG;
      typedef typename Acme2CylGeometrySwitch::GeometrySwitch<GEO_ORIG>::type GEO;
      GEO geo(eit->geometry());

      // Calculate element anisotropy from diagonal
      typename GEO::GlobalCoordinate diam = geo.corner(geo.corners()-1);
      diam -= geo.corner(0);

      double anisotropy = diam[1] / diam[0];
      if(anisotropy < 1)
        anisotropy = 1. / anisotropy;

      minAnisotropy =  std::min(minAnisotropy, anisotropy);
      maxAnisotropy =  std::max(maxAnisotropy, anisotropy);

      for(int j=0; j<eit->geometry().corners(); j++)
      {
        CoordType nodePos = geo.corner(j);

        // Do volume scaling only for extracellular nodes with a certain distance from the membrane
        if(geo.center()[1] > yStart && params.doVolumeScaling())
        {
          // Set reference volume the first time an element with center[1] > yStart is visited
          if(reference_volume < 0)
          {
            reference_volume = geo.volume();

            debug_info << "Doing threshold volume scaling, threshold yStart=" << yStart
                << ", reference volume @" << geo.center() << ": " << reference_volume << std::endl;
          }

          //int index = -1;
          bool exists = (nodeVolumesMap.count(nodePos) > 0);

          // Now calculate volume for this node
  //        double r1 = std::abs(eit->geometry().corner(0)[1]);
  //        double r2 = std::abs(eit->geometry().corner(eit->geometry().corners()-1)[1]);
  //        double h = std::abs(eit->geometry().corner(
  //            eit->geometry().corners()-1)[0] - eit->geometry().corner(0)[0]);
  //        double volume = con_pi * r2 * r2 * h - con_pi * r1 * r1 * h;

          // Use generic volume from geometry (valid for both cylinder and 2d geometry)
          double volume = geo.volume();

          // If not, add it to the vector of node positions
          if(! exists)
          {
            // Use minimum volume of visited elements
            nodeVolumesMap[nodePos] = volume/reference_volume;

          } else {
            nodeVolumesMap[nodePos] = std::min(nodeVolumesMap[nodePos], volume/reference_volume);
          }
        } else {
          // Use constant volume scaling of 1 for node (identity)
          nodeVolumesMap[nodePos] = 1.0;
        }
      }

    }
    int nodeVolumesSize = nodeVolumesMap.size();
    // This processor-local check is correct for a 2D tensor grid only (really?)
    assert(nodeVolumesSize == gv.size(GV::dimension));

    nodeVolumesSize = gv.comm().sum(nodeVolumesSize);
    assert(nodeVolumesSize >= params.nNodes()); // There are duplicates at processor boundaries

    minAnisotropy =  gv.comm().min(minAnisotropy);
    maxAnisotropy =  gv.comm().max(maxAnisotropy);

    debug_info << "Min/max grid cell anisotropy: " << minAnisotropy << " / " << maxAnisotropy << std::endl;

//    // Test the mapping nodeIndex <-> nodeCoordinate <-> nodeVolume
//    for (typename GV::template Codim<GV::dimension>::Iterator nit = gv.template begin<GV::dimension>();
//        nit != gv.template end<GV::dimension>(); ++nit)
//    {
//      //int index = Ax1LFSTools::getNodeIndex(params.getNodePositions(), nit->geometry().center());
//      debug_jochen << "Node position " << nit->geometry().center()
//          //<< ", index: " << index
//          << ", volume: " << getNodeVolume(nit->geometry().center())
//          << std::endl;
//    }

  }

  void initPartition()
  {
    int nX_this = nGridElementsPerAxis[0];
    int nY_this = nGridElementsPerAxis[1];

    // Communicate number of elements in x-direction on each process
    gv.comm().allgather(&nX_this,1,&partition_x[0]);
    // Same for y (should be equal on all processors!)
    gv.comm().allgather(&nY_this,1,&partition_y[0]);

    debug_info << "Communicated grid partition, # elements on each process:" << std::endl;
    for(int i=0; i<partition_x.size(); i++)
    {
      if(i>0)
      {
        offset_x[i] = offset_x[i-1]+partition_x[i-1];
        offset_y[i] = offset_y[i-1]+partition_y[i-1];
      }
      debug_info << "-p" << i << ": x=" << partition_x[i] << " (offset " << offset_x[i] << ") | y="
          << partition_y[i] << " (offset " << offset_y[i] << ")" << std::endl;
    }
  }

  // FIXME Generalize this to multiple membrane element layers
  //! Setup vector of permittivities (one value for each membrane element)
  void initPermittivities()
  {
    permittivities.resize(membGV.size(0));

    bool smoothPermittivities = params.membrane.get("smoothPermittivities", false);

    const typename Acme2CylParameters::MembraneGroups& membGroups = params.getMembraneGroups();
    assert(membGroups.size() > 0);

    typename Acme2CylParameters::MembraneGroups::const_iterator git = membGroups.begin();

    double xStart = membGV.template begin<0>()->geometry().corner(0)[0];
    debug_jochen << "Calculating membrane permittivities, xStart = " << xStart << std::endl;

    // Advance to next interval for which the right interval border is greater than the x-coordinate of the
    // first vertex 'xStart' on this processor; make sure we don't run out of the membrane group vector of course.
    while(std::get<2>(*git) < xStart && git != membGroups.end()-1)
    {
      //debug_jochen << "Skipping group [" << std::get<1>(*git) << " " << std::get<2>(*git) << "]" << std::endl;
      git++;
    }

    // Coordinate of the currently regarded group transition
    bool right = (xStart+1e-6 > std::get<1>(*git)) ? true : false;
    double currentGroupTransition = right ? std::get<2>(*git) : std::get<1>(*git);
    double width = std::get<2>(*git) - std::get<1>(*git);
    std::string groupName = std::get<0>(*git);

    //debug_jochen << "Initial group transition: " << currentGroupTransition << " (" << (right ? "right" : "left") << ")" << std::endl;

    double otherWidth = (git+1 != membGroups.end() ? (std::get<2>(*(git+1)) - std::get<1>(*(git+1))) : width);
    double dxTransition = params.membrane.sub(groupName).get("dx_transition", std::min(width, otherWidth));
    int nTransition = params.membrane.sub(groupName).get("n_transition", 10);

    // This is the transition interval over which the permittivity discontinuity shall be smoothed
    double dTransition = dxTransition*nTransition;
    bool smoothTransition = smoothPermittivities && params.membrane.sub(groupName).get("smoothTransition", false);

    // Special case with only one membrane group: There is nothing to smooth here!
    if(membGroups.size() == 1)
    {
      dTransition = 0.0;
      smoothTransition = false;
    }

    double default_permittivity = membrane.getPermittivity();
    double default_dMemb = params.dMemb();

    debug_info << "Membrane permittivities: " << std::endl;
    debug_info << "----------------------------------" << std::endl;
    for(SubDomainElementIterator_All mit = membGV.template begin<0>(); mit != membGV.template end<0>(); ++mit)
    {
      typename SubDomainElement::Geometry::GlobalCoordinate center = mit->geometry().center();
      int sdElemIndex = getSubDomainElementIndex(*mit);

      assert(sdElemIndex < permittivities.size());

      // Flag indicating if current coordinate is exceeding the transition interval border
      bool updated = false;
      // Check if current position exceeds a 'left' interval border
      if(! right && (center[0] > currentGroupTransition + dTransition) && git != membGroups.end()-1)
      {
        currentGroupTransition = std::get<2>(*git);
        right = true;
        updated = true;

      // Check if current position exceeds a 'right' interval border
      } else if(right && (center[0] > currentGroupTransition) && git != membGroups.end()-1)
      {
        // Advance to next interval
        git++;

        groupName = std::get<0>(*git);
        width = std::get<2>(*git) - std::get<1>(*git);

        smoothTransition = smoothPermittivities && params.membrane.sub(groupName).get("smoothTransition", false);
        // Take right boundary right away
        if(!smoothTransition && git != membGroups.end()-1)
        {
          currentGroupTransition = std::get<2>(*git);
          right = true;
        } else {
          currentGroupTransition = std::get<1>(*git);
          right = false;
        }
        updated = true;
      }

      // If transition coordinate has been updated, also update transition width 'dTransition'
      if(updated)
      {
        // Is the transition to be smoothed and do neighboring membrane groups actually exist?
        if(smoothTransition && (right ? (git+1 != membGroups.end()) : (git != membGroups.begin())))
        {
          typename Acme2CylParameters::MembraneGroups::const_iterator git_other = (right ? git+1 : git-1);
          double otherWidth = std::get<2>(*git_other) - std::get<1>(*git_other);

          dxTransition = params.membrane.sub(groupName).get("dx_transition", std::min(width, otherWidth));
          nTransition = params.membrane.sub(groupName).get("n_transition", 10);
          dTransition = dxTransition*nTransition;
        } else {
          // No transition interval
          dTransition = 0.0;
        }

//        debug_jochen << "Update group transition to " << currentGroupTransition << " (" << (right ? "right" : "left")
//            << ")" << std::endl;
//        debug_jochen << " dxTransition = " << dxTransition << ", n_transition = " << nTransition
//            << " => dTransition = " << dTransition << std::endl;
      }

      // Is there a real smooth transition (with a transition interval > 0)?
      if(smoothTransition && dTransition > 0.0)
      {
        // Now determine permittivities at either side of the transition coordinate
        typename Acme2CylParameters::MembraneGroups::const_iterator git_left = (right ? git : git-1);
        typename Acme2CylParameters::MembraneGroups::const_iterator git_right = (right ? git+1 : git);

        std::string groupName_left = std::get<0>(*git_left);
        std::string groupName_right = std::get<0>(*git_right);

        double permittivity_left = params.membrane.sub(groupName_left).get("permittivity", default_permittivity);
        double permittivity_right = params.membrane.sub(groupName_right).get("permittivity", default_permittivity);

        // Group-specific membrane thickness
        double dMemb_left = params.membrane.sub(groupName_left).get("d_memb", default_dMemb);
        double dMemb_right = params.membrane.sub(groupName_right).get("d_memb", default_dMemb);

        // Scale permittivity by real membrane thickness to get effective permittivity
        permittivity_left *= (default_dMemb / dMemb_left);
        permittivity_right *= (default_dMemb / dMemb_right);

        // Check consistency of parameters: When smoothing permittivities, only one side of each membrane group transition
        // should have flag 'smoothTransition' set!
        bool smoothTransition_left = smoothPermittivities && params.membrane.sub(groupName_left).get("smoothTransition", false);
        bool smoothTransition_right = smoothPermittivities && params.membrane.sub(groupName_right).get("smoothTransition", false);

        if(smoothTransition_left && smoothTransition_right)
          DUNE_THROW(Dune::Exception, "When smoothing permittivities at membrane group transitions, only one side of two groups"
              << " should have 'smoothTransition = yes' in config file! Found flags set for groups '" << groupName_left << "' and '"
              << groupName_right << "'!" );

        //debug_jochen << "permittivity_left = " << permittivity_left << ", permittivity_right = " << permittivity_right << std::endl;

        // This should contain a value between 0 and 1
        double x_scale = 0;

        // If x coordinate is in the transition interval, calculate its local position x_scale \in [0 1] in this interval
        if(std::abs(center[0]-currentGroupTransition) < dTransition)
        {
          // Coordinate is in transition interval!
          if(right)
            x_scale = (center[0] - (currentGroupTransition-dTransition)) / (dTransition);
          else
            x_scale = (center[0] - currentGroupTransition) / (dTransition);
        } else {
          if(right) x_scale = 0.0;
          else x_scale = 1.0;
        }
        //debug_jochen << " => x_scale = " << x_scale << std::endl;
        assert(x_scale >= 0 && x_scale <= 1.0);

        // Smooth permittivities at group transitions by 'smooth step' strategy!
        permittivities[sdElemIndex] = permittivity_left + (permittivity_right-permittivity_left) * (x_scale*x_scale*(3-2*x_scale));

      } else {
        // Simply take the permittivity belonging to this membrane group

        int groupIndex = getGroupIndex(*mit);
        std::string groupName = getGroupName(groupIndex);

        double permittivity = params.membrane.sub(groupName).get("permittivity", default_permittivity);
        double dMemb = params.membrane.sub(groupName).get("d_memb", default_dMemb);

        permittivity *= (default_dMemb / dMemb);

        permittivities[sdElemIndex] = permittivity;
      }
      debug_info << "Membrane element [" << sdElemIndex << "] @" << center
              << ": permittivity = " << permittivities[sdElemIndex] << std::endl;
    }
    debug_info << "----------------------------------" << std::endl;
  }

  const GV& gv;
  const SubGV& elecGV;
  const SubGV& membGV;

  ElementSubdomainMapper& elemSubdomainMapper;
  Ax1ElementGroupMapper& elemGroupMapper;
  Electrolyte<T> electro;
  Membrane<T> membrane;
  Acme2CylParameters& params;
  const T TIME_SCALE;
  const T LENGTH_SCALE;

  ElementMapper elementMapper;  // Dune entity->index mapper
  SubDomainElementMapper elecElementMapper;
  SubDomainElementMapper membElementMapper;

  std::vector<int> nGridElementsPerAxis;
  std::vector<int> nOffsetPerAxis;
  std::vector<std::vector<int> > nSubdomainElementsPerAxis;

  NodeVolumesMap nodeVolumesMap;
  MembraneInterfaces* membraneInterfaces;
  MembraneInterfaceMap* membraneInterfaceMap_primary;
  MembraneInterfaceMap* membraneInterfaceMap_bidirectional;
  MembraneInterfaceMap* membraneInterfaceMap_self;

  std::vector<int> partition_x;
  std::vector<int> partition_y;
  std::vector<int> offset_x;
  std::vector<int> offset_y;

  std::vector<std::string> groupNames;
  T dt;

  std::vector<T> permittivities;
};

#endif /* DUNE_AX1_ACME2CYL_PHYSICS_HH */
