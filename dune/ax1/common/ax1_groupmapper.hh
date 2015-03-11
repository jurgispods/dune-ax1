/*
 * ax1_groupmapper.hh
 *
 *  Created on: Jun 4, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_GROUPMAPPER_HH
#define DUNE_AX1_GROUPMAPPER_HH

#include <map>

class Ax1ElementGroupMapper
{
  public:

    Ax1ElementGroupMapper()
    {
    }

    void setGroup(int elementIndex, int groupIndex)
    {
      elementIndexMap[elementIndex] = groupIndex;
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

#endif /* DUNE_AX1_GROUPMAPPER_HH */
