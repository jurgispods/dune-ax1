
#ifndef DUNE_AX1_MEMBRANE_PHYSICS_HH
#define DUNE_AX1_MEMBRANE_PHYSICS_HH

#include <vector>

#include<dune/ax1/common/electrolyte.hh>
#include<dune/ax1/common/membrane.hh>
#include<dune/ax1/common/element_properties.hh>
#include<dune/ax1/common/element_subdomain_mapper.hh>

template<class T>
class Physics {
public:
  Electrolyte<T> elecl, elecr;
  Membrane<T> membrane;
  
  Physics ( Electrolyte<T>& elecl_, Membrane<T> memb_, Electrolyte<T>& elecr_ )
  : elecl( elecl_ ), membrane( memb_ ), elecr( elecr_ ), groups(2), elementGroupProperties(2)
  {
    groups[0] = "solution";
    groups[1] = "membrane";

    initElementGroupProperties();
  }
  
  ~Physics () {}
  
  // potential boundary values
  
  void setBoundaries( T boundl_, T boundr_ )
  {
    boundl = boundl_;
    boundr = boundr_;
  }
  
  T getBoundaryl() { return boundl; }
  T getBoundaryr() { return boundr; }
  
  // mebrane start/end position
  
  void setMembranePos( T posL, T posR )
  {
    membranePosL = posL;
    membranePosR = posR;
  }
  
  T getMembranePosL() { return membranePosL; }
  T getMembranePosR() { return membranePosR; }
  
  // membrane maximum thickness in limit case for given dPhi
  
  T getMaxThicknessForDPhi ( T dl, int z, T dPhi )
  {
    return dl * membrane.getPermittivity() / sqrt(elecl.getPermittivity()) * dPhi * exp( 0.25 * z * dPhi );
  }
  
  // analytical solutions ###############################################
  
  // membrane thickness zero
  T potential ( T x, T dl, int z )
  {
    return 2.0 / z * log( 1.0 + 0.5 * z * x / dl );
  }
  
  // membrane thickness d
  T potential2 ( T x, T dl, int z, T d, T dPhi )
  {
    return 2.0 / z * log( exp( 0.25 * z * dPhi ) +
                         0.5 * z * ( x - 0.5 * d ) / ( dl * sqrt(elecl.getPermittivity()) ) );
  }
  


  // element properties ################################################
  
  void initElementGroupProperties()
  {
    for(int i=0; i<elementGroupProperties.size(); ++i)
    {
      ElementProperties elemProps;

      elemProps.setName(groups[i]);
      if(groups[i] == "membrane")
      {
        elemProps.setMembrane(true);
      }

      elementGroupProperties[i] = elemProps;
    }
  }

  const ElementProperties& getElementProperties(int elementIndex)
  {
    return elementGroupProperties[elemSubdomainMapper.map(elementIndex)];
  }

  void setElementSubdomainMapper(ElementSubdomainMapper elemSubdomainMapper_)
  {
    elemSubdomainMapper = elemSubdomainMapper_;
  }


private:

  T boundl, boundr;
  T membranePosL, membranePosR;

  std::vector<std::string> groups;
  std::vector<ElementProperties> elementGroupProperties;

  ElementSubdomainMapper elemSubdomainMapper;
};

#endif /* DUNE_AX1_MEMBRANE_PHYSICS_HH */
