/*
 * element_to_subdomain_mapper.hh
 *
 *  Created on: Jun 1, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ELEMENT_SUBDOMAIN_MAPPER_HH
#define DUNE_AX1_ELEMENT_SUBDOMAIN_MAPPER_HH

#include <map>

class ElementSubdomainMapper
{
  public:

    ElementSubdomainMapper()
    {
    }

    void setSubdomain(int elementIndex, int subdomainIndex)
    {
      elementIndexMap[elementIndex] = subdomainIndex;
    }

    int map(const int elementIndex) const
    {
      // TODO Find out if this is also working on non-GCC compilers
      return elementIndexMap.at(elementIndex);
    }

    int size() const
    {
      return elementIndexMap.size();
    }


  private:
    std::map<int, int> elementIndexMap;


};


#endif /* DUNE_AX1_ELEMENT_SUBDOMAIN_MAPPER_HH */
