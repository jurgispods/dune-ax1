/*
 * ax1_tensorgrid_transformation.hh
 *
 *  Created on: Apr 30, 2013
 *      Author: jpods
 */

#include <dune/grid/geometrygrid/coordfunction.hh>

#ifndef DUNE_AX1_TENSORGRID_TRANSFORMATION_HH
#define DUNE_AX1_TENSORGRID_TRANSFORMATION_HH

#define CACHE_VERTICES 1

template<typename ctype>
class Ax1TensorGridTransformation
  : public Dune::AnalyticalCoordFunction<ctype,2,2,Ax1TensorGridTransformation<ctype> >
{
  typedef Ax1TensorGridTransformation<ctype> This;
  typedef Dune::AnalyticalCoordFunction<ctype,2,2,This> Base;

public:
  typedef typename Base::DomainVector DomainVector;
  typedef typename Base::RangeVector RangeVector;

  static constexpr ctype EPS = 1e-6;

  Ax1TensorGridTransformation(ctype xmax_, ctype ymax_,
      const std::vector<ctype>& x_, const std::vector<ctype>& y_)
  : xmax(xmax_),
    ymax(ymax_),
    xTensor(x_),
    yTensor(y_),
    h_x(xmax_/(xTensor.size()-1)),
    h_y(ymax_/(yTensor.size()-1))
  {
    //std::cout << "h_x = " << h_x << std::endl;
    //std::cout << "h_y = " << h_y << std::endl;

    for(int i=0; i<xTensor.size(); i++)
    {
      ctype x_yasp = i*h_x;
      xMap[x_yasp] = xTensor[i];

      //std::cout << "Map x: " << x_yasp << " -> " << xTensor[i] << std::endl;
    }
    for(int i=0; i<yTensor.size(); i++)
    {
      ctype y_yasp = i*h_y;
      yMap[y_yasp] = yTensor[i];

      //std::cout << "Map y: " << y_yasp << " -> " << yTensor[i] << std::endl;
    }
  }

  struct fuzzy_less_than
  {

    bool operator() (const ctype& x1, const ctype& x2) const
    {
      if ( x1 < x2-EPS ) return true;   // is less than
      //if ( x1 > x2+EPS ) return false;  // is greater than
      return false; // is equal
    }
  };

  void evaluate ( const DomainVector &x, RangeVector &y ) const
  {
#if CACHE_VERTICES
//    if(xMap.count(x[0]) == 0)
//      std::cerr << "ERROR, x = " << x[0] << " not found!" << std::endl;
//    if(yMap.count(x[1]) == 0)
//      std::cerr << "ERROR, y = " << x[1] << " not found!" << std::endl;
    y[0] = xMap.at(x[0]);
    y[1] = yMap.at(x[1]);
#else
    int xi = x[0]/h_x;
    int yi = x[1]/h_y;

    assert(xi < xTensor.size());
    assert(yi < yTensor.size());

    y[0] = xTensor[xi];
    y[1] = yTensor[yi];
#endif
    //std::cout << "x = " << x << " -> " << y << std::endl;
  }

private:
  const ctype xmax;
  const ctype ymax;
  const std::vector<ctype>& xTensor;
  const std::vector<ctype>& yTensor;
  const ctype h_x;
  const ctype h_y;
  std::map<ctype,ctype,fuzzy_less_than> xMap;
  std::map<ctype,ctype,fuzzy_less_than> yMap;
};

#endif /* DUNE_AX1_TENSORGRID_TRANSFORMATION_HH */
