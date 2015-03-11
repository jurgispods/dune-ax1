#ifndef DUNE_AX1_ACME1MD_PHYSICS_HH
#define DUNE_AX1_ACME1MD_PHYSICS_HH

#include <vector>
#include <iomanip>

#include <dune/ax1/common/electrolyte.hh>
#include <dune/ax1/common/membrane.hh>
#include <dune/ax1/common/element_subdomain_mapper.hh>
#include <dune/ax1/common/membrane_interface_mapper.hh>
#include <dune/ax1/common/ax1_subgrid_hostelement_mapper.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/acme1MD/common/acme1MD_parametertree.hh>

template<typename GV, typename T, typename ConfigTraits>
class Physics {
public:

  // General typedefs
  typedef GV GridView;
  typedef T FieldType;
  typedef ConfigTraits Traits;

  typedef typename GV::Grid::SubDomainGrid               SubGrid;
  typedef typename GV::Grid::SubDomainGrid::LeafGridView SubGV;

  // Grid-specific typedefs
  typedef typename GV::template Codim<0>::Entity Element;
  typedef typename GV::template Codim<0>::EntityPointer ElementPointer;
  typedef typename GV::template Codim<0>::Iterator ElementIterator;
  typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

  typedef typename GV::Intersection ElementIntersection;
  typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 1> IntersectionMapper;
  typedef typename Element::LeafIntersectionIterator ElementIntersectionIterator;
  
  // Dune::SubGrid typedefs
  typedef typename SubGV::template Codim<0>::Entity SubDomainElement;
  typedef typename SubGV::template Codim<0>::EntityPointer SubDomainElementPointer;
  typedef typename SubGV::template Codim<0>::Iterator SubDomainElementIterator;
  typedef Dune::SingleCodimSingleGeomTypeMapper<SubGV, 0> SubDomainElementMapper;

  typedef typename SubDomainElement::LeafIntersectionIterator SubDomainElementIntersectionIterator;
  typedef typename SubDomainElementIntersectionIterator::Intersection SubDomainElementIntersection;
  typedef Dune::SingleCodimSingleGeomTypeMapper<SubGV, 1> SubDomainElementIntersectionMapper;

  typedef SubGridHostElementMapper<ElementPointer,ElementMapper> HostElementMapper;

  typedef typename Membrane<T>::ChannelSet ChannelSet;


  Physics ( const GV& gv_, const SubGV& elecGV_, const SubGV& membGV_, ElementSubdomainMapper& elementGroupMapper_,
           Electrolyte<T>& electro_, Membrane<T>& memb_, Acme1MDParameters& params_,
           const T& timeScale_, const T& lengthScale_)
  : gv(gv_),
    elecGV(elecGV_),
    membGV(membGV_),
    elemSubdomainMapper(elementGroupMapper_),
    electro(electro_),
    membrane(memb_),
    params(params_),
    TIME_SCALE(timeScale_),
    LENGTH_SCALE(lengthScale_),
    elementMapper(gv),
    elecElementMapper(elecGV),
    membElementMapper(membGV),
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
  bool isMembraneInterface(const ElementIntersection& is)
  {
    if(is.boundary())
      return false;
    else
      // Exactly one of the neighboring elements is a membrane element
      return (isMembrane(*is.inside()) != isMembrane(*is.outside()));
  }

  //! For SubGridIntersections: Check if it is a membrane interface
  bool isMembraneInterface(const SubDomainElementIntersection& is)
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

  int getSubDomainMembraneElementIndex(const ElementIntersection& is)
  {
    assert(isMembraneInterface(is));

    if(isMembrane(*is.inside()))
    {
      return getSubDomainElementIndex(*membGV.grid().subDomainEntityPointer(*is.inside()));
    } else {
      // We are on the membrane subdomain, inside element is membrane element
      return getSubDomainElementIndex(*membGV.grid().subDomainEntityPointer(*is.outside()));
    }
  }

  int getSubDomainMembraneElementIndex(const SubDomainElementIntersection& is)
  {
    if(elecGV.contains(*is.inside()))
    {
      // We are on the electrolyte subdomain, outside element is membrane element
      Element& mde_out = *(gv.grid().multiDomainIntersection(is).outside());
      return getSubDomainElementIndex(*membGV.grid().subDomainEntityPointer(mde_out));
    } else {
      // We are on the membrane subdomain, inside element is membrane element
      return getSubDomainElementIndex(*is.inside());
    }
  }

  //! \brief Get membrane element index for a given membrane intersection
  int getMembraneElementIndex(const ElementIntersection& is)
  {
    assert(isMembraneInterface(is));

    if(isMembrane(*is.inside()))
    {
      // We are on the electrolyte subdomain, outside element is membrane element
      return getElementIndex(*is.inside());
    } else {
      // We are on the membrane subdomain, inside element is membrane element
      return getElementIndex(*is.outside());
    }
  }

  //! \brief Get membrane element index for a given membrane intersection
  int getMembraneElementIndex(const SubDomainElementIntersection& is)
  {
    assert(isMembraneInterface(is));

    if(elecGV.contains(*is.inside()))
    {
      // We are on the electrolyte subdomain, outside element is membrane element
      return getElementIndex(*gv.grid().multiDomainIntersection(is).outside());
    } else {
      // We are on the membrane subdomain, inside element is membrane element
      return getElementIndex(*gv.grid().multiDomainIntersection(is).inside());
    }
  }

  //! \brief get group index
  int getSubdomainIndex (const Element& e) const
  {
    return getSubdomainIndex(getElementIndex(e));
  }

  //! \brief get group index
  int getSubdomainIndex (const SubDomainElement& e) const
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
  

  //! \brief get host element index of a subgrid entity
  int getElementIndex (const SubDomainElement& e) const
  {
    return getElementIndex(gv.grid().multiDomainEntity(e));
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


  // ================================ BEGIN SubGrid <-> HostGrid stuff =====================================

  const ElementIntersectionIterator
  getOppositeMembraneIntersection(const ElementIntersection& is)
  {

    DUNE_THROW(Dune::Exception, "Opposite membrane interface could not be found!");
  }

  const SubDomainElementIntersectionIterator
  getOppositeMembraneIntersection(const SubDomainElementIntersection& is)
  {

    DUNE_THROW(Dune::Exception, "Opposite membrane interface could not be found!");
  }

  template<typename DGF_POT>
  void getMembranePotential(const Element& me, DGF_POT& dgfPot, typename DGF_POT::Traits::RangeType& membPot)
  {
    getMembranePotentialJump(me, dgfPot, membPot);
    membPot *= -1;
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


  /*
  //! \brief Calculated potential jump between given intersection and opposite membrane interface
  //! intersection using a discrete grid function for the calculated potential
  template<typename DGF_CON>
  void getMembraneConcentrationJump(const Element& me, DGF_CON& dgfCon,
    typename DGF_CON::Traits::RangeType& conCytosol, typename DGF_CON::Traits::RangeType& conExtracellular)
  {
    if(not isMembrane(me))
    {
      DUNE_THROW(Dune::Exception, "HUHN-HAHN!!!!!!!");
    }

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
            // "iit->outside" is the Cytosol element
            dgfCon.evaluate(*iit->outside(), iit->geometryInOutside().center(), conCytosol);
            break;
          }
          case ES:
          {
            // "iit->outside" is the extracellular element
            dgfCon.evaluate(*iit->outside(), iit->geometryInOutside().center(), conExtracellular);
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

    debug_jochen << "Concentration inside/outside ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_jochen << "[" << ION_NAMES[j] << "] = " << conCytosol[j] << "/" << conExtracellular[j] << "  ";
    }
    debug_jochen << std::endl;

    debug_verb << "Concentration jump: ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_verb << "[" << ION_NAMES[j] << "] = " << (conCytosol[j] - conExtracellular[j]) << "  ";
    }
    debug_verb << std::endl;
  }
  */

  //! \brief Calculated concentration jump across this membrane element. It calculates the concentration on both
  //! ends (i.e., neighboring electrolyte boundary) and saves the oriented difference in the parameter conJump.
  //! The sign of potJump gives the orientation of the jump in positive y-axis direction
  template<typename DGF_CON>
  void getMembraneConcentrationRatio(const Element& me, DGF_CON& dgfCon,
    typename DGF_CON::Traits::RangeType& conRatio)
  {
    if(not isMembrane(me))
    {
      DUNE_THROW(Dune::Exception, "HUHN-HAHN!!!!!!!");
    }

    int countMembraneInterfaces = 0;
    typename DGF_CON::Traits::RangeType con1(0.0), con2(0.0);
    typename DGF_CON::Traits::DomainType x1(0.0), x2(0.0);
    conRatio = 0.0;

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
          case ES:
          {
            if(countMembraneInterfaces == 1)
            {
              x1 = iit->outside()->geometry().global(iit->geometryInOutside().center());
              dgfCon.evaluate(*iit->outside(), iit->geometryInOutside().center(), con1);
            } else {
              x2 = iit->outside()->geometry().global(iit->geometryInOutside().center());
              dgfCon.evaluate(*iit->outside(), iit->geometryInOutside().center(), con2);
            }
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

    // Sign of y difference (= y-axis direction)
    int ySign = Dune::sign(x1[1] - x2[1]);
    debug_verb << "x1 = " << x1 << ", x2 = " << x2
      << " => sign = "<< ySign  << (ySign < 0 ? "↓" : "↑") << std::endl;

    if(ySign > 0)
    {
      conRatio = con2;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        conRatio[j] /= con1[j];
      }
    }
    if(ySign < 0)
    {
      conRatio = con1;
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        conRatio[j] /= con2[j];
      }
    }

    debug_verb << "Concentration bottom/top ";
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_verb << "[" << ION_NAMES[j] << "] = " << con1[j] << "/" << con2[j] << "  ";
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

  //! \brief convert dimensionless potential to units of mV
  T convertFrom_mV (const T& potential) const
  {
    return potential * con_e / (1.0e3 * con_k * electro.getTemperature());
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
        gfIonicStrength.evaluate(*eit, eit->geometry().local(eit->geometry().center()), ionicStrength);

        //debug_jochen << "ionic strength @" << eit->geometry().center() << " = " << ionicStrength << std::endl;

        RF debyeLengthLocal = std::sqrt( con_eps0 * permittivity * con_k * electro.getTemperature() /
            ( 2.0 * con_e * con_e * electro.getStdCon() * ionicStrength) );

        debyeLength = std::min(debyeLength, debyeLengthLocal);
      }
    }
    return debyeLength;
  }

  void gridInfo()
  {
    const SubGV& elecGV = gv.grid().subDomain(0).leafView();
    const SubGV& membGV = gv.grid().subDomain(1).leafView();

    debug_verb << "================ Grid info: =========================" << std::endl;
    debug_verb << "Subgrid INSIDE (" << elecGV.size(0) << " elements)" << std::endl;
    for(SubDomainElementIterator elec_it = elecGV.template begin<0>();
        elec_it != elecGV.template end<0>(); ++elec_it)
    {
      int subElementIndex = elecGV.indexSet().index(*elec_it);
      debug_verb << "subGridElement[elec] #" << subElementIndex
         << " -> element #" << getElementIndex(*elec_it)
         << " [ " <<  elec_it->geometry().corner(0)
         << " - " << elec_it->geometry().center()
         << " - " <<  elec_it->geometry().corner(1) << "]"
         << std::endl;
    }
    debug_verb << "Subgrid OUTSIDE (" << membGV.size(0) << " elements)" << std::endl;
    for(SubDomainElementIterator memb_it = membGV.template begin<0>();
            memb_it != membGV.template end<0>(); ++memb_it)
    {
      int subElementIndex = membGV.indexSet().index(*memb_it);
      debug_verb << "subGridElement[memb] #" << subElementIndex
         << " -> element #" << getElementIndex(*memb_it)
         << " [ " <<  memb_it->geometry().corner(0)
         << " - " << memb_it->geometry().center()
         << " - " <<  memb_it->geometry().corner(1) << "]"
         << std::endl;
    }

    debug_verb << std::endl;
    debug_verb << "subGV.contains() test:" << std::endl;
    // Test of member function "contains()"
    for(SubDomainElementIterator elec_it = elecGV.template begin<0>();
                elec_it != elecGV.template end<0>(); ++elec_it)
    {
      debug_verb << elecGV.contains(*elec_it) << " ";
    }
    for(SubDomainElementIterator memb_it = membGV.template begin<0>();
                      memb_it != membGV.template end<0>(); ++memb_it)
    {
      debug_verb << elecGV.contains(*memb_it) << " ";
    }
    debug_verb << std::endl;
    for(SubDomainElementIterator elec_it = elecGV.template begin<0>();
                elec_it != elecGV.template end<0>(); ++elec_it)
    {
      debug_verb << membGV.contains(*elec_it) << " ";
    }
    for(SubDomainElementIterator memb_it = membGV.template begin<0>();
                          memb_it != membGV.template end<0>(); ++memb_it)
    {
      debug_verb << membGV.contains(*memb_it) << " ";
    }
    debug_verb << std::endl;

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
  const Acme1MDParameters& getParams() const
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


  const GV& gv;
  const SubGV& elecGV;
  const SubGV& membGV;

  ElementSubdomainMapper& elemSubdomainMapper;
  Electrolyte<T> electro;
  Membrane<T> membrane;
  Acme1MDParameters& params;
  const T TIME_SCALE;
  const T LENGTH_SCALE;

  ElementMapper elementMapper;  // Dune entity->index mapper
  SubDomainElementMapper elecElementMapper;
  SubDomainElementMapper membElementMapper;

  std::vector<std::string> groups;
  std::valarray<T> nodePositions;
  std::valarray<T> cellCenterPositions;
  T dt;
};

#endif /* DUNE_AX1_ACME1MD_PHYSICS_HH */
