#ifndef DUNE_AX1_ACME1_PHYSICS_HH
#define DUNE_AX1_ACME1_PHYSICS_HH

#include <vector>
#include <iomanip>

#include <dune/ax1/common/electrolyte.hh>
#include <dune/ax1/common/membrane.hh>
#include <dune/ax1/common/element_subdomain_mapper.hh>
#include <dune/ax1/common/membrane_interface_mapper.hh>
#include <dune/ax1/common/ax1_subgrid_hostelement_mapper.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/acme1/common/acme1_parametertree.hh>

template<typename GV, typename T, typename ConfigTraits, typename SubGV>
class Physics {
public:

  // General typedefs
  typedef GV GridView;
  typedef T FieldType;
  typedef ConfigTraits Traits;
  typedef SubGV SubGridView;

  // Grid-specific typedefs
  typedef typename GV::template Codim<0>::Entity Element;
  typedef typename GV::template Codim<0>::EntityPointer ElementPointer;
  typedef typename GV::template Codim<0>::Iterator ElementIterator;
  typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

  typedef typename GV::Intersection ElementIntersection;
  typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 1> IntersectionMapper;
  typedef typename Element::LeafIntersectionIterator ElementIntersectionIterator;
  
  // Dune::SubGrid typedefs
  typedef typename SubGV::template Codim<0>::Entity SubGridElement;
  typedef typename SubGV::template Codim<0>::Iterator SubGridElementIterator;
  typedef Dune::SingleCodimSingleGeomTypeMapper<SubGV, 0> SubElementMapper;

  typedef typename SubGridElement::LeafIntersectionIterator SubElementIntersectionIterator;
  typedef typename SubElementIntersectionIterator::Intersection SubElementIntersection;
  typedef Dune::SingleCodimSingleGeomTypeMapper<SubGV, 1> SubElementIntersectionMapper;

  typedef SubGridHostElementMapper<ElementPointer,ElementMapper> HostElementMapper;

  typedef typename Membrane<T>::ChannelSet ChannelSet;


  Physics ( const GV& gv_, ElementSubdomainMapper& elementGroupMapper_,
           Electrolyte<T>& electro_, Membrane<T>& memb_, Acme1Parameters& params_,
           const T& timeScale_, const T& lengthScale_)
  : gv(gv_),
    elemSubdomainMapper(elementGroupMapper_),
    electro(electro_),
    membrane(memb_),
    params(params_),
    TIME_SCALE(timeScale_),
    LENGTH_SCALE(lengthScale_),
    elementMapper(gv),
    membraneInterfaceMapper(elementMapper),
    groups(params.membrane.getSubKeys())
  {
    // First two groups are always the (extracellular/intracellular) solutions;
    // Insert those at the beginning
    std::vector<std::string>::iterator git = groups.begin();
    git = groups.insert(git, "solution_in");
    groups.insert(git, "solution_ex");
    initVectors();
  }

  void initVectors()
  {
    // resize to number of vertices
    nodePositions.resize(gv.size(0)+1);
    // resize to number of elements
    cellCenterPositions.resize(gv.size(GV::dimension));
  }

  void initPosition(std::vector<T> coords)
  {
    assert(coords.size() == nodePositions.size());
    for(int i=0; i<nodePositions.size(); ++i)
    {
      nodePositions[i] = coords[i];
    }
    //TODO Works in 1D only!
    for(int i=0; i<cellCenterPositions.size(); ++i)
    {
      cellCenterPositions[i] = (coords[i+1] - coords[i])/2;
    }
  }

  /*! \brief This method tests if an element is a membrane element or not.
   */
  bool isMembrane(const Element& e) const
  {
    return isMembrane(getSubdomainIndex(getElementIndex(e)));
  }

  bool hasMembraneInterface(const Element& e)
  {
    int elemIndex = getElementIndex(e);

    // Check if we find an existing mapping in the membraneInterfaceMapper
    return membraneInterfaceMapper.contains(elemIndex);
  }

  bool hasMembraneInterface(const SubGridElement& e, const SubGV& subGV)
  {
    int hostElemIndex_Inside = getElementIndex(e, subGV);

    // Check if we find an existing mapping in the membraneInterfaceMapper
    return membraneInterfaceMapper.contains(hostElemIndex_Inside);
  }

  //! For host grid Intersections: Check if it is a membrane interface
  bool isMembraneInterface(const ElementIntersection& is)
  {
    // Exactly one of the neighboring elements is a membrane element
    return (isMembrane(*is.inside()) != isMembrane(*is.outside()));
  }

  //! For SubGridIntersections: Check if it is a membrane interface
  bool isMembraneInterface(const SubElementIntersection& is, const SubGV& subGV)
  {
    //TODO Check geometry! Are we really on the membrane?

    // According to our current assumptions, the membrane interface is always located at
    // the boundary of a subgrid!
    if(is.boundary())
    {
      int hostElemIndex_Inside = getElementIndex(*is.inside(), subGV);
      return membraneInterfaceMapper.contains(hostElemIndex_Inside);
    }
    return false;
  }

  int getLocalMembraneElementIndex(const int elemIndex)
  {
    const std::map<int,int> membElements = membrane.getChannelSet().getMembraneElements();
    std::map<int,int>::const_iterator it = membElements.find(elemIndex);
    if(it == membElements.end())
    {
      DUNE_THROW(Dune::Exception, "Element #" << elemIndex << " is not a membrane element!");
    }
    return it->second;
  }

  //! \brief Get membrane element index for a given membrane intersection
  int getMembraneElementIndex(const SubElementIntersection& is, const SubGV& subGV)
  {
    // This is slightly dumb: Get opposite membrane intersection first (living on the host grid)
    const SubElementIntersectionIterator iit = getOppositeMembraneIntersection(is, subGV);
    // Then get the element index of the 'outside' (membrane) element!
    const ElementPointer mep(iit->outside());
    return getElementIndex(*mep);
  }

  //! \brief get group index
  int getSubdomainIndex (const Element& e) const
  {
    return getSubdomainIndex(getElementIndex(e));
  }

  //! \brief get group name for a given group index
  std::string getGroupName(const Element& e) const
  {
    return groups[getSubdomainIndex(e)];
  }


  //! \brief get element index
  int getElementIndex (const Element& e) const
  {
    return elementMapper.map(e);
  }
  
  //! \brief get element index
  //! Compatibility method for subgrid functionality
  template<typename Element>
  int getElementIndex (const Element& e, const GV& gv) const
  {
    return elementMapper.map(e);
  }

#if USE_SUBGRID==1
  // ax1-subgrid stuff

  //! \brief get host element index of a subgrid entity
  int getElementIndex (const SubGridElement& e, const SubGV& gv) const
  {
    if(! gv.contains(e))
    {
      DUNE_THROW(Dune::Exception, "SubGridElement and SubGridView do not match!");
    }
    return getElementIndex(gv.getHostElement(e));
  }
#endif
#if USE_SUBGRID==2
  // dune-subgrid stuff

  //! \brief get element index of a subgrid entity
  int getElementIndex (const SubGridElement& e, const SubGV& gv) const
  {
    Element& hostEntity = *(gv.grid().template getHostEntity<0>(e));
    return elementMapper.map(hostEntity);
  }
#endif


  //!" \brief Find opposite interfaces located on each side of the membrane
  void addMembraneInterface(const ElementPointer& e)
  {
    typedef typename Element::LeafIntersectionIterator IntersectionIterator;

    ElementPointer ep(e);

    int countMembraneInterfaces = 0;
    for(IntersectionIterator iit = e->ileafbegin(); iit !=e->ileafend(); ++iit)
    {
      bool notMembrane = ! isMembrane(*iit->outside());
      if(notMembrane)
      {
        countMembraneInterfaces++;
        if(countMembraneInterfaces == 1)
        {
          ep = iit->outside();
        }
        if(countMembraneInterfaces == 2) {
          membraneInterfaceMapper.addInterfaceElements(ep, iit->outside());
        }
      }
    }
    if(countMembraneInterfaces != 2)
    {
      DUNE_THROW(Dune::Exception, "Wrong number of neighboring non-membrane elements!");
    }
  }
  
  // ================================ BEGIN SubGrid <-> HostGrid stuff =====================================

  const Element& getOppositeMembraneInterfaceElement(const Element& e)
  {
    int elemIndex = getElementIndex(e);
    if(isMembrane(e))
    {
      DUNE_THROW(Dune::Exception, "This method only works for non-membrane elements!");
    }
    if(not hasMembraneInterface(e))
    {
      DUNE_THROW(Dune::Exception, "This method only works for elements with a membrane interface!");
    }
    return *membraneInterfaceMapper.map(elemIndex);
  }

  const Element& getOppositeMembraneInterfaceElement(const SubGridElement& e, const SubGV& subGV)
  {
    int hostElemIndex = getElementIndex(e, subGV);
    if(isMembrane(hostElemIndex))
    {
      DUNE_THROW(Dune::Exception, "This method only works for non-membrane elements!");
    }
    if(not hasMembraneInterface(e, subGV))
    {
      DUNE_THROW(Dune::Exception, "This method only works for elements with a membrane interface!");
    }
    return *membraneInterfaceMapper.map(hostElemIndex);
  }

  const ElementIntersectionIterator
  getOppositeMembraneIntersection(const ElementIntersection& is)
  {
    if(! isMembraneInterface(is))
    {
      DUNE_THROW(Dune::Exception, "The given intersection is not a membrane interface!");
    }
    const Element& eOpposite = getOppositeMembraneInterfaceElement(*is.inside());

    SubElementIntersectionIterator iit=eOpposite.ileafbegin();
    for( ; iit != eOpposite.ileafend(); ++iit)
    {
      // Check if this (host grid!) intersection is on the membrane interface
      if(isMembraneInterface(*iit))
      {
        return iit;
      }
    }
    DUNE_THROW(Dune::Exception, "Opposite membrane interface could not be found!");
  }

  const SubElementIntersectionIterator
  getOppositeMembraneIntersection(const SubElementIntersection& is, const SubGV& subGV)
  {
    if(! isMembraneInterface(is, subGV))
    {
      DUNE_THROW(Dune::Exception, "The given intersection is not a membrane interface!");
    }
    const SubGridElement& eOpposite = getOppositeMembraneInterfaceElement(*is.inside(),subGV);

    SubElementIntersectionIterator iit=eOpposite.ileafbegin();
    for( ; iit != eOpposite.ileafend(); ++iit)
    {
      // Check if this (host grid!) intersection is on the membrane interface
      if(isMembraneInterface(*iit))
      {
        return iit;
      }
    }
    DUNE_THROW(Dune::Exception, "Opposite membrane interface could not be found!");
  }

  template<typename DGF_POT>
  void getMembranePotential(const Element& me, DGF_POT& dgfPot, typename DGF_POT::Traits::RangeType& membPot)
  {
    getMembranePotentialJump(me, dgfPot, membPot);
    membPot *= -1;
  }


#if USE_SUBGRID==0

  //! \brief Compatibility method for subgrid usage
  template<typename DGF_POT>
  typename DGF_POT::Traits::RangeType getMembranePotentialJump(const SubElementIntersection& is,
        const SubGV& subGV, DGF_POT& dgfPot, typename DGF_POT::Traits::RangeType& potJump)
  {
    potJump = 0.0;
  }
#else

  //! \brief Calculated potential jump between given intersection and opposite membrane interface
  //! intersection using a discrete grid function for the calculated potential
  template<typename DGF_POT>
  typename DGF_POT::Traits::RangeType getMembranePotentialJump(const SubElementIntersection& is,
      const SubGV& subGV, DGF_POT& dgfPot, typename DGF_POT::Traits::RangeType& potJump)
  {
    const ElementIntersectionIterator& iitOpposite = getOppositeMembraneIntersection(is, subGV);
    const ElementIntersection& isOpposite = *iitOpposite;

    // Return potential jump at neighboring membrane element
    return getMembranePotentialJump(*isOpposite.outside(), dgfPot, potJump);
  }

  //! \brief Calculated potential jump between given intersection and opposite membrane interface
  //! intersection using a discrete grid function for the calculated potential
  template<typename DGF_POT>
  typename DGF_POT::Traits::RangeType getMembranePotentialJump(const ElementIntersection& is, DGF_POT& dgfPot,
      typename DGF_POT::Traits::RangeType& potJump)
  {
    if(! isMembraneInterface(is))
    {
      DUNE_THROW(Dune::Exception, "This element intersection is not a membrane interface!");
    }
    if(isMembrane(*is.inside()))
    {
      return getMembranePotentialJump(*is.inside(), dgfPot, potJump);
    } else {
      return getMembranePotentialJump(*is.outside(), dgfPot, potJump);
    }
  }

  //! \brief Calculated potential jump between given intersection and opposite membrane interface
  //! intersection using a discrete grid function for the calculated potential
  template<typename DGF_POT>
  void getMembranePotentialJump(const Element& me, DGF_POT& dgfPot, typename DGF_POT::Traits::RangeType& potJump)
  {
    if(not isMembrane(me))
    {
      DUNE_THROW(Dune::Exception, "HUHN-HAHN!!!!!!!");
    }

    typename DGF_POT::Traits::RangeType pot_Inside = 0.0;
    typename DGF_POT::Traits::RangeType pot_Outside = 0.0;
    potJump = 0.0;

    int countMembraneInterfaces = 0;
    // Iterate intersections
    for(ElementIntersectionIterator iit = me.ileafbegin(); iit != me.ileafend(); ++iit)
    {
      if(isMembraneInterface(*iit))
      {
        countMembraneInterfaces++;
        switch(getSubdomainIndex(*iit->outside()))
        {
          case CYTOSOL:
          {
            dgfPot.evaluate(me, iit->geometryInInside().center(), pot_Inside);
            break;
          }
          case ES:
          {
            dgfPot.evaluate(me, iit->geometryInInside().center(), pot_Outside);
            break;
          }
          default:
            DUNE_THROW(Dune::Exception, "Outside element of membrane interface is neither CYTOSOL nor ES!");
        }
      }

    }
    // Assume we have only one layer of membrane elements
    // => exactly two membrane interfaces for each membrane element!
    assert(countMembraneInterfaces == 2);

    potJump = pot_Outside;
    potJump -= pot_Inside;
    debug_verb << "Potential jump: " << potJump
        << " [" << convertTo_mV(potJump) << " mV]"<< std::endl;
  }
#endif
#if USE_SUBGRID==0

  //! \brief Compatibility method for subgrid usage
  template<typename DGF_CON>
  void getMembraneConcentrationJump(const SubElementIntersection& is,
        const SubGV& subGV, DGF_CON& dgfCon, typename DGF_CON::Traits::RangeType& conJump)
  {
    conJump = 0.0;
  }
#else

  //! \brief Calculated potential jump between given intersection and opposite membrane interface
  //! intersection using a discrete grid function for the calculated potential
  template<typename DGF_CON>
  void getMembraneConcentrationJump(const SubElementIntersection& is,
      const SubGV& subGV, DGF_CON& dgfCon, typename DGF_CON::Traits::RangeType& conJump,
      typename DGF_CON::Traits::RangeType& conUp)
  {
    const ElementIntersectionIterator& iitOpposite = getOppositeMembraneIntersection(is, subGV);
    const ElementIntersection& isOpposite = *iitOpposite;

    // Return potential jump at neighboring membrane element
    return getMembraneConcentrationJump(*isOpposite.outside(), dgfCon, conJump, conUp);
  }

  //! \brief Calculated potential jump between given intersection and opposite membrane interface
  //! intersection using a discrete grid function for the calculated potential
  template<typename DGF_CON>
  void getMembraneConcentrationJump(const ElementIntersection& is, DGF_CON& dgfCon,
      typename DGF_CON::Traits::RangeType& conJump, typename DGF_CON::Traits::RangeType& conUp)
  {
    if(! isMembraneInterface(is))
    {
      DUNE_THROW(Dune::Exception, "This element intersection is not a membrane interface!");
    }
    if(isMembrane(*is.inside()))
    {
      return getMembraneConcentrationJump(*is.inside(), dgfCon, conJump, conUp);
    } else {
      return getMembraneConcentrationJump(*is.outside(), dgfCon, conJump, conUp);
    }
  }

  //! \brief Calculated potential jump between given intersection and opposite membrane interface
  //! intersection using a discrete grid function for the calculated potential
  template<typename DGF_CON>
  void getMembraneConcentrationJump(const Element& me, DGF_CON& dgfCon,
    typename DGF_CON::Traits::RangeType& conJump, typename DGF_CON::Traits::RangeType& conUp)
  {
    if(not isMembrane(me))
    {
      DUNE_THROW(Dune::Exception, "HUHN-HAHN!!!!!!!");
    }

    typename DGF_CON::Traits::RangeType con_Inside(0.0);
    typename DGF_CON::Traits::RangeType con_Outside(0.0);
    conJump = 0.0;
    conUp = 0.0;

    int countMembraneInterfaces = 0;
    // Iterate intersections
    for(ElementIntersectionIterator iit = me.ileafbegin(); iit != me.ileafend(); ++iit)
    {
      if(isMembraneInterface(*iit))
      {
        countMembraneInterfaces++;
        switch(getSubdomainIndex(*iit->outside()))
        {
          // We need to evaluate concentrations on the outside of the intersection, as they are 0 on the membrane
          case CYTOSOL:
          {
            // naming is a little confusing:
            // "iit->outside" is the Cytosol element on which "con_Inside" is evaluated
            dgfCon.evaluate(*iit->outside(), iit->geometryInOutside().center(), con_Inside);
            break;
          }
          case ES:
          {
            dgfCon.evaluate(*iit->outside(), iit->geometryInOutside().center(), con_Outside);
            break;
          }
          default:
            DUNE_THROW(Dune::Exception, "Outside element of membrane interface is neither CYTOSOL nor ES!");
        }
      }

    }
    // Assume we have only one layer of membrane elements
    // => exactly two membrane interfaces for each membrane element!
    assert(countMembraneInterfaces == 2);

    conJump = con_Outside;
    conJump -= con_Inside;

    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      conUp[j] = std::max(con_Inside[j], con_Outside[j]);
    }

    debug_verb << "Concentration inside/outside ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_verb << "[" << ION_NAMES[j] << "] = " << con_Inside[j] << "/" << con_Outside[j] << "  ";
    }
    debug_verb << std::endl;

    debug_verb << "Concentration jump: ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_verb << "[" << ION_NAMES[j] << "] = " << conJump[j] << "  ";
    }
    debug_verb << std::endl;

    /*
    debug_verb << "Upwind concentrations: ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_verb << ION_NAMES[j] << "_up = " << conUp[j] << "  ";
    }
    debug_verb << std::endl;
    */
  }
#endif

  // ================================ END SubGrid <-> HostGrid stuff =======================================

  //! \brief convert dimensionless potential to units of mV
  std::valarray<T> convertTo_mV (const std::valarray<T>& potential) const
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



  template<typename GF_IONIC_STRENGTH>
  double getDebyeLength(GF_IONIC_STRENGTH& gfIonicStrength)
  {
    typedef double Real;
    const int dim = GV::dimension;
    const int intorder = 2;

    typedef typename GF_IONIC_STRENGTH::Traits::DomainType DT;
    typedef typename GF_IONIC_STRENGTH::Traits::DomainFieldType DF;
    typedef typename GF_IONIC_STRENGTH::Traits::RangeType RT;
    typedef typename GF_IONIC_STRENGTH::Traits::RangeFieldType RF;

    RF debyeLength = 1e100; // aha!

    for (ElementIterator eit=gv.template begin<0>(); eit!=gv.template end<0>(); ++eit)
    {
      int elemIndex = elementMapper.map(*eit);
      int subdomainIndex = getSubdomainIndex(elemIndex);

      if ( not isMembrane(subdomainIndex) ) // debye length not evaluated on membrane
      {
        Real permittivity = getPermittivity(subdomainIndex);

        RT ionicStrength(0.0);
        gfIonicStrength.evaluate(*eit, eit->geometry().center(), ionicStrength);

        RF debyeLengthLocal = sqrt( con_eps0 * permittivity * con_k * electro.getTemperature() /
            ( 2.0 * con_e * con_e * electro.getStdCon() * ionicStrength) );

        debyeLength = std::min(debyeLength, debyeLengthLocal);
      }
    }
    return debyeLength;
  }

  void gridInfo(SubGV& subGV_Inside, SubGV& subGV_Outside)
  {
    debug_verb << "================ Grid info: =========================" << std::endl;
    debug_verb << "Subgrid INSIDE (" << subGV_Inside.size(0) << " elements)" << std::endl;
    for(SubGridElementIterator sit_Inside = subGV_Inside.template begin<0>();
        sit_Inside != subGV_Inside.template end<0>(); ++sit_Inside)
    {
      int subElementIndex = subGV_Inside.indexSet().index(*sit_Inside);
      debug_verb << "subGridElement[inside] #" << subElementIndex
         << " -> element #" << getElementIndex(*sit_Inside,subGV_Inside)
         << " [ " <<  sit_Inside->geometry().corner(0)
         << " - " << sit_Inside->geometry().center()
         << " - " <<  sit_Inside->geometry().corner(1) << "]"
         << std::endl;
    }
    debug_verb << "Subgrid OUTSIDE (" << subGV_Outside.size(0) << " elements)" << std::endl;
    for(SubGridElementIterator sit_Outside = subGV_Outside.template begin<0>();
            sit_Outside != subGV_Outside.template end<0>(); ++sit_Outside)
    {
      int subElementIndex = subGV_Outside.indexSet().index(*sit_Outside);
      debug_verb << "subGridElement[outside] #" << subElementIndex
         << " -> element #" << getElementIndex(*sit_Outside,subGV_Outside)
         << " [ " <<  sit_Outside->geometry().corner(0)
         << " - " << sit_Outside->geometry().center()
         << " - " <<  sit_Outside->geometry().corner(1) << "]"
         << std::endl;
    }

    debug_verb << std::endl;
    debug_verb << "subGV.contains() test:" << std::endl;
    // Test of member function "contains()"
    for(SubGridElementIterator sit_Inside = subGV_Inside.template begin<0>();
                sit_Inside != subGV_Inside.template end<0>(); ++sit_Inside)
    {
      debug_verb << subGV_Inside.contains(*sit_Inside) << " ";
    }
    for(SubGridElementIterator sit_Outside = subGV_Outside.template begin<0>();
                      sit_Outside != subGV_Outside.template end<0>(); ++sit_Outside)
    {
      debug_verb << subGV_Inside.contains(*sit_Outside) << " ";
    }
    debug_verb << std::endl;
    for(SubGridElementIterator sit_Inside = subGV_Inside.template begin<0>();
                sit_Inside != subGV_Inside.template end<0>(); ++sit_Inside)
    {
      debug_verb << subGV_Outside.contains(*sit_Inside) << " ";
    }
    for(SubGridElementIterator sit_Outside = subGV_Outside.template begin<0>();
                          sit_Outside != subGV_Outside.template end<0>(); ++sit_Outside)
    {
      debug_verb << subGV_Outside.contains(*sit_Outside) << " ";
    }
    debug_verb << std::endl;

    membraneInterfaceMapper.info();
    debug_verb << "=====================================================" << std::endl;
  }



  //! \brief debug output
  void info()
  {
    debug_info << "-------- Physics info -----------" << std::endl;
    debug_info << "Elements groups:" << std::endl;
    for(int i=0; i<groups.size(); ++i)
    {
      debug_info << "G[" << i << "] = " << groups[i] << std::endl;
    }
    debug_info << std::endl;
    debug_info << "Membrane elements:" << std::endl;
    const std::map<int,int>& membraneElements = membrane.getChannelSet().getMembraneElements();
    std::map<int, int>::const_iterator it;
    for(it = membraneElements.begin(); it != membraneElements.end(); ++it)
    {
      int subdomainIndex = getSubdomainIndex(it->first);
      debug_info << "E[" << it->first << "] --> ME[" << it->second << "]  -  G["
          << subdomainIndex << "] = " << getGroupName(subdomainIndex)
          << std::endl;
    }
    debug_info << "---------------------------------" << std::endl;
  }

   template<typename GridView>
   int numberSubDomainInterfaces(const GridView& gv_) const
   {
     return params.useMembrane() ? 2 : 0;
   }

   Membrane<T>& getMembrane()
   {
     return membrane;
   }

   Electrolyte<T>& getElectrolyte()
   {
     return electro;
   }

   //! \brief get number of ion species
   int numOfSpecies () const
   {
     return electro.numOfSpecies();
   }

   int nElements() const
   {
     return elemSubdomainMapper.size();
   }

   //! \brief get permittivity
   T getPermittivity (const int elemIndex) const
   {
     int subdomainIndex = getSubdomainIndex(elemIndex);
     if ( not isMembrane(subdomainIndex) )      // electrolyte
     {
       return electro.getPermittivity();
     } else {// membrane
       return membrane.getPermittivity();
     }
   }

  //! \brief get parameter tree (representing config file)
  const Acme1Parameters& getParams() const
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
    { return electro.getIonName( ionSpecies ); }

  //! \brief get squared length constant for poisson equation
  T getPoissonConstant() const
  {
    return electro.getPoissonConstant();
  }

  T getTimeStep()
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

  std::valarray<T>& getPosition()
  {
    return nodePositions;
  }

  std::valarray<T>& getCellCenterPositions()
  {
    return cellCenterPositions;
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


private:

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

  //! \brief get group index
  int getSubdomainIndex (const int elementIndex) const
  {
    return elemSubdomainMapper.map( elementIndex );
  }

  //! \brief get group name for a given group index
  std::string getGroupName(const int subdomainIndex) const
  {
   return groups[subdomainIndex];
  }

  //! \brief get diffusion coefficient
  T getDiffCoeff (int subdomainIndex, int ionSpecies, int elemIndex) const
  {
    // We don't calculate concentrations inside the membrane!
    //assert(isMembrane(subdomainIndex) == false);
    if(isMembrane(subdomainIndex)) return 0.0;

    return electro.getDiffConst( ionSpecies );
  }


  const GV & gv;
  ElementSubdomainMapper& elemSubdomainMapper;
  Electrolyte<T> electro;
  Membrane<T> membrane;
  Acme1Parameters& params;
  const T TIME_SCALE;
  const T LENGTH_SCALE;

  ElementMapper elementMapper;  // Dune entity->index mapper

  // Membrane interface neighbor element mapper
  MembraneInterfaceMapper<ElementPointer,ElementMapper> membraneInterfaceMapper;

  std::vector<std::string> groups;
  std::valarray<T> nodePositions;
  std::valarray<T> cellCenterPositions;
  T dt;
};

#endif /* DUNE_AX1_ACME1_PHYSICS_HH */
