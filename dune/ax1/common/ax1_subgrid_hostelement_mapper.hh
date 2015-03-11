/*
 * subgroup_element_mapper.hh
 *
 *  Created on: Feb 2, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_SUBGRID_ELEMENT_INDEX_MAPPER_HH
#define DUNE_AX1_SUBGRID_ELEMENT_INDEX_MAPPER_HH

#include <map>

template<typename EntityPointer, typename EntityMapper>
class SubGridHostElementMapper
{
  public:

    typedef typename EntityPointer::Entity Entity;

    SubGridHostElementMapper(EntityMapper& eMapper_)
    : eMapper(eMapper_)
    {}

    void addElement(int subElementIndex, const EntityPointer& ep)
    {
      subIndexHostElementMap.insert(std::make_pair(subElementIndex, ep));
    }

    bool contains(const int subElemIndex) const
    {
      //typename InterfaceMap::const_iterator it = membraneInterfaceMap.find(elemIndex);
      //return ! (it == InterfaceMap::end);
      return subIndexHostElementMap.count(subElemIndex);
    }

    const EntityPointer& map(const int subElementIndex) const
    {
      // TODO Find if this is also working on non-GCC compilers
      return subIndexHostElementMap.at(subElementIndex);
    }

    int size() const
    {
      return subIndexHostElementMap.size();
    }




  private:
    std::map<int, EntityPointer> subIndexHostElementMap;
    EntityMapper& eMapper;


};

#endif /* DUNE_AX1_SUBGRID_ELEMENT_INDEX_MAPPER_HH */
