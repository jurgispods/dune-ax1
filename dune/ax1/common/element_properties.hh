/*
 * element_properties.hh
 *
 *  Created on: Jun 1, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ELEMENT_PROPERTIES_HH
#define DUNE_AX1_ELEMENT_PROPERTIES_HH

#include <string>
#include <iostream>

class ElementProperties
{
  public:

    ElementProperties() :
      is_membrane(false)
    {
    }

    bool isMembrane() const
    {
      return is_membrane;
    }

    bool setMembrane(bool is_membrane_)
    {
      is_membrane = is_membrane_;
    }

    std::string getName() const
    {
      return name;
    }

    void setName(std::string name_)
    {
      name = name_;
    }

  private:
    bool is_membrane;
    std::string name;

};


#endif /* DUNE_AX1_ELEMENT_PROPERTIES_HH */
