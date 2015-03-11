/*
 * membrane_interface_mapper.hh
 *
 *  Created on: Feb 13, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_MEMBRANE_INTERFACE_MAPPER_HH
#define DUNE_AX1_MEMBRANE_INTERFACE_MAPPER_HH

#include <map>

template<typename EntityPointer, typename EntityMapper>
class MembraneInterfaceMapper
{
  public:

    MembraneInterfaceMapper(EntityMapper eMapper_)
    : eMapper(eMapper_)
    {}

    typedef typename EntityPointer::Entity Entity;
    typedef std::map<int, EntityPointer> InterfaceMap;

    void addInterfaceElements(const EntityPointer& ep_Inside, const EntityPointer& ep_Outside)
    {
      int index_Inside = eMapper.map(*ep_Inside);
      int index_Outside = eMapper.map(*ep_Outside);
      membraneInterfaceMap.insert(std::make_pair(index_Inside, ep_Outside));
      membraneInterfaceMap.insert(std::make_pair(index_Outside, ep_Inside));
    }

    /*
    void addInterfaceElements(const Entity& e_Inside, const Entity& e_Outside)
    {
      EntityPointer ep_Inside(e_Inside);
      EntityPointer ep_Outside(e_Outside);
      membraneInterfaceMap.insert(std::make_pair(ep_Inside, ep_Outside));
      membraneInterfaceMap.insert(std::make_pair(ep_Outside, ep_Inside));
    }*/

    bool contains(const int elemIndex) const
    {
      //typename InterfaceMap::const_iterator it = membraneInterfaceMap.find(elemIndex);
      //return ! (it == InterfaceMap::end);
      return membraneInterfaceMap.count(elemIndex);
    }

    const EntityPointer& map(const int elemIndex) const
    {
      // TODO Find if this is also working on non-GCC compilers
      return membraneInterfaceMap.at(elemIndex);
    }

    const EntityPointer& map(const Entity& e) const
    {
      return map(eMapper.map(e));
    }

    int size() const
    {
      return membraneInterfaceMap.size();
    }

    void info() const
    {
      debug_verb << "----------- MembraneInterfaceMapper ----------------" << std::endl;
      for(typename InterfaceMap::const_iterator it = membraneInterfaceMap.begin();
          it != membraneInterfaceMap.end(); ++it)
      {
        debug_verb << it->first << " <=> "
            << eMapper.map(*it->second)
            << "[" << it->second->geometry().center() <<"]" << std::endl;
      }
      debug_verb << "----------------------------------------------------" << std::endl;
    }


  private:
    std::map<int, EntityPointer> membraneInterfaceMap;
    EntityMapper& eMapper;


};

#endif /* DUNE_AX1_MEMBRANE_INTERFACE_MAPPER_HH */
