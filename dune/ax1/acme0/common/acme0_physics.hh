#ifndef DUNE_AX1_ACME0_PHYSICS_HH
#define DUNE_AX1_ACME0_PHYSICS_HH

#include <vector>
#include <iomanip>

#include <dune/ax1/common/electrolyte.hh>
#include <dune/ax1/common/membrane.hh>
#include <dune/ax1/common/element_subdomain_mapper.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/acme0/common/acme0_parametertree.hh>

template<typename GV, typename T, typename ConfigTraits>
class Physics {
public:

  typedef GV GridView;
  typedef T FieldType;

  typedef typename GV::template Codim<0>::Entity Element;
  typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;
  typedef typename GV::template Codim<0>::Iterator ElementIterator;

  typedef ConfigTraits Traits;
  
  Physics ( const GV& gv_, ElementMapper elementMapper_,
           Electrolyte<T>& electro_, Membrane<T>& memb_, Acme0Parameters& params_,
           const T& timeScale_, const T& lengthScale_)
  : gv(gv_),
    elementMapper(elementMapper_),
    electro(electro_),
    membrane(memb_),
    params(params_),
    groups(params.membrane.getSubKeys()),
    TIME_SCALE(timeScale_),
    LENGTH_SCALE(lengthScale_)
  {
    // First two groups are always the (extracellular/intracellular) solutions;
    // Insert those at the beginning
    std::vector<std::string>::iterator git = groups.begin();
    git = groups.insert(git, "solution_in");
    groups.insert(git, "solution_ex");
    initVectors();
  }
  
  ~Physics () {}

  Membrane<T>& getMembrane()
  {
    return membrane;
  }

  Electrolyte<T>& getElectrolyte()
  {
    return electro;
  }

  void setElementSubdomainMapper(ElementSubdomainMapper elemSubdomainMapper_)
  {
    elemSubdomainMapper = elemSubdomainMapper_;
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

  /*! \brief This method tests if an element is a membrane element or not.
   */
  template<typename Element>
  bool isMembrane(const Element& e) const
  {
    return isMembrane(getSubdomainIndex(getElementIndex(e)));
  }

  //! \brief get group index
  int getSubdomainIndex (const int elementIndex) const
  {
    return elemSubdomainMapper.map( elementIndex );
  }
  
  //! \brief get group index
  template<typename Element>
  int getSubdomainIndex (const Element& e) const
  {
    return getSubdomainIndex(getElementIndex(e));
  }

  //! \brief get group name for a given group index
  std::string getGroupName(const Element& e) const
  {
    return groups[getSubdomainIndex(e)];
  }

  //! \brief get group name for a given group index
  std::string getGroupName(const int subdomainIndex) const
  {
    return groups[subdomainIndex];
  }

  //! \brief get element index
  template<typename Element>
  int getElementIndex (const Element& e) const
  {
    return elementMapper.map(e);
  }
  
  //! \brief get permittivity
  T getPermittivity (const int subdomainIndex) const
  {
    if ( not isMembrane(subdomainIndex) )      // electrolyte
    {
      return electro.getPermittivity();
    } else {// membrane
      return membrane.getPermittivity();
    }
  }
  
  /*
  //! \brief set permittivity vector
  void setPermittivityVector()
  {
    permittivity.resize(gv.size(1));
    debug_verb << "permittivity.size = " << permittivity.size() << std::endl;
    
    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    for ( ElementIterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
    {
      int elemIndex = elementMapper.map(*it);
      permittivity[elemIndex]=getPermittivity(getSubdomainIndex(elemIndex));
      //std::cout << elemIndex << " " << permittivity[elemIndex] << std::endl;
      //debug_verb << "element #" << elemIndex << ", permittivity = " << permittivity[elemIndex] << std::endl;
    }
    // Handle last node
    permittivity[gv.size(1)-1] = permittivity[gv.size(1)-2];
    //debug_verb << "element #" << gv.size(1)-1 << ", permittivity = " << permittivity[gv.size(1)-1] << std::endl;
  }
  */
  
  //! \brief get parameter tree (representing config file)
  const Acme0Parameters& getParams() const
  {
    return params;
  }

  //! \brief get diffusion coefficient
  T getDiffCoeff (int ionSpecies, int elemIndex) const
  {
    int subdomainIndex = getSubdomainIndex(elemIndex);
    return getDiffCoeff(subdomainIndex, ionSpecies, elemIndex);
  }

  //! \brief get diffusion coefficient
  T getDiffCoeff (int subdomainIndex, int ionSpecies, int elemIndex) const
  {
    
    if ( not isMembrane(subdomainIndex) )
    {
      return electro.getDiffConst( ionSpecies );
    } else {

      T g = membrane.getDiffCoeff( ionSpecies, elemIndex);
      
      // calculate diffusion coeffcient D from conductance g
      // (given in mS/cm2 - my guess) --> correct! [jürgen]
      
      g *= 10.0; // convert g to S/m2 (SI-units)
      T membraneThickness = 4.0 * LENGTH_SCALE;
      T D = con_k*electro.getTemperature()*g*membraneThickness/
        (con_e*con_e*electro.getValence(ionSpecies)*
         electro.getValence(ionSpecies)*electro.getStdCon());
      
      //debug_verb << "## " << std::setprecision(16) << ionSpecies << "  " << D << std::endl;

      // using old diffusion constants for testing
      return membrane.getDiffConst( ionSpecies );
      //return D;
    }
    //return 0.0; // freeze everything
  }
  
  //! \brief Update membrane channels
  void updateDiffCoeffs(int subdomainIndex, int elemIndex, T dt, T deltaPot)
  {
    if ( isMembrane(subdomainIndex) )
    {
      membrane.updateDiffCoeff(elemIndex, dt, deltaPot);
    }
  }

  //! \brief get valence
  T getValence (const int ionSpecies) const
  {
    return electro.getValence( ionSpecies );
  }
  
  //! \brief get name of ion species
  std::string getIonName ( const unsigned int ionSpecies ) const
    { return electro.getIonName( ionSpecies ); }
  
  //! \brief get squared length constant for poisson equation
  T getPoissonConstant() const
  {
    return electro.getPoissonConstant();
  }
  
  void initVectors()
  {
    nodePositions.resize(gv.size(0)+1);
    cellCenterPositions.resize(gv.size(GV::dimension));
    //chargeDensity.resize(gv.size(0)+1);
    //permittivity.resize(gv.size(0)+1);
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

  template<typename GF_IONIC_STRENGTH>
  double getDebyeLength(GF_IONIC_STRENGTH& gfIonicStrength)
  {
    typedef double Real;
    const int dim = GV::dimension;
    const int intorder = 2;

    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef typename Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;
    typedef typename GF_IONIC_STRENGTH::Traits::DomainType DT;
    typedef typename GF_IONIC_STRENGTH::Traits::DomainFieldType DF;
    typedef typename GF_IONIC_STRENGTH::Traits::RangeType RT;
    typedef typename GF_IONIC_STRENGTH::Traits::RangeFieldType RF;

    ElementMapper elementMapper(gv);

    RF debyeLength = 1e100; // aha!

    for(int j=0; j<NUMBER_OF_SPECIES; ++j)
    {

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
    }

    return debyeLength;
  }

  void gridInfo()
  {
      /*
    typedef typename GV::Codim<0>::Iterator Iterator;
    typedef typename GV::Traits::template Codim<0>::EntityPointer EntityPointer;

    for(Iterator it = gv.begin(); it<gv.end(); ++it)
    {
      EntityPointer e = it;
      int elemIndex = getElementIndex(*ep);

      chargeDensity[elemIndex] = initialChargeDensity.evaluate(e, );
    }
        */

    //typedef typename GF::Traits::GridViewType GV;
    typedef typename GV::template Codim<0>::Entity Entity;
    typedef typename GV::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GV::template Codim<0>::Iterator ElementLeafIterator;

    //const GV& gv = initalChargeDensity.getGridView();
    int solSize = gv.size(0);

    //typename GV::Traits::DomainType x;
    //typename GV::Traits::RangeType y;
    double x, y;

    int i=0;
    for (ElementLeafIterator nit=gv.template begin<0>(); nit!=gv.template end<0>(); ++nit) {

      Entity& e = *nit;
      int elemIndex = getElementIndex(*nit);
      int subdomainIndex = getSubdomainIndex(elemIndex);
      x = e.geometry().local(e.geometry().center());

      //initalChargeDensity.evaluate(e, x, y);
      //chargeDensity[elemIndex] = y;
      y =  getDiffCoeff (subdomainIndex, 0, elemIndex);
      debug_verb << "initial chargeDensity [" << elemIndex  << "] @ x = "
          << e.geometry().center() << " : "<< y << std::endl;
    }

  }

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

  /*
  std::valarray<T>& getChargeDensity()
  {
    return chargeDensity;
  }

  std::valarray<T>& getPermittivity()
  {
    return permittivity;
  }
  */

  const T getTimeScale() const
  {
    return TIME_SCALE;
  }

  const T getLengthScale() const
  {
    return LENGTH_SCALE;
  }

  /*
  template<typename IntersectionGeometry>
  IntersectionGeometry getOppositeMembraneInterface(const IntersectionGeometry& face)
  {
    typedef typename IntersectionGeometry::EntityPointer EntityPointer;
    typedef typename IntersectionGeometry::Entity        Entity;
    typedef typename Entity::IntersectionIterator        IntersectionIterator;

    EntityPointer pElemInside = face.inside();


    for(IntersectionIterator iit = pElemInside->ileafbegin(); iit !=pElemInside->ileafend(); ++iit)
    {
      // TODO Find opposite intersection
      // TODO Get index of intersection
      int index = ?
      IntersectionGeometry<Intersection> ig(*iit,index);

      return ig;
    }
    DUNE_THROW(Dune::Exception, "No corresponding interface found!");

  }
  */

  const GV& gridView() const
  {
    return gv;
  }



private:
  const GV & gv;
  ElementMapper elementMapper;
  Electrolyte<T> electro;
  Membrane<T> membrane;
  Acme0Parameters& params;
  std::vector<std::string> groups;
  std::valarray<T> nodePositions;
  std::valarray<T> cellCenterPositions;
  //std::valarray<T> chargeDensity, permittivity;
  ElementSubdomainMapper elemSubdomainMapper;
  T dt;
  const T TIME_SCALE;
  const T LENGTH_SCALE;
};

#endif /* DUNE_AX1_ACME0_PHYSICS_HH */
