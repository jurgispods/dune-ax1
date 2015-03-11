/*
 * ax1_gridtools.hh
 *
 *  Created on: Aug 27, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_GRIDTOOLS_HH
#define DUNE_AX1_GRIDTOOLS_HH

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/ax1_gridvector.hh>

class Ax1GridTools
{

  public:

    template<typename GridType>
    static void prepareGrid(Dune::GridFactory<GridType>& factory,
        GridVector<typename GridType::ctype>& x,
        GridVector<typename GridType::ctype>& y,
        const Dune::GeometryType& gt)
    {
      //insert vertices
      typename GridVector<typename GridType::ctype>::iterator itxend = x.end();
      typename GridVector<typename GridType::ctype>::iterator ityend = y.end();

      int vertexCount = 0;
      for (typename GridVector<typename GridType::ctype>::iterator ity = y.begin(); ity!=ityend; ++ity)
      {
        for (typename GridVector<typename GridType::ctype>::iterator itx = x.begin(); itx!=itxend; ++itx)
        {
          //set vertex position
          Dune::FieldVector<typename GridType::ctype,2> pos;
          pos[0] = *itx;
          pos[1] = *ity;
          //debug_jochen << "Inserting vertex " << (vertexCount++) << ": "  << pos << " into UG grid" << std::endl;
          factory.insertVertex(pos);
        }
      }

      //insert elements
      for (unsigned int i=0; i<y.size()-1;i++)
      {
        for (unsigned int j=0; j<x.size()-1;j++)
        {
          //set node positions
          std::vector<unsigned int> pos;
          pos.push_back(j+i*x.size());
          pos.push_back(j+1+i*x.size());
          pos.push_back(j+(i+1)*x.size());
          pos.push_back(j+1+(i+1)*x.size());

          //debug_jochen << "Inserting element with vertices ";
          //Output::printVector(pos);

          factory.insertElement(gt,pos);
        }
      }

    }

};

#endif /* DUNE_AX1_GRIDTOOLS_HH */
