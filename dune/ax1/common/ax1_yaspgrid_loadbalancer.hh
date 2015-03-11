/*
 * ax1_yaspgrid_loadbalancer.hh
 *
 *  Created on: May 2, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_YASPGRID_LOADBALANCER_HH
#define DUNE_AX1_YASPGRID_LOADBALANCER_HH

/*
  With this class you can specify how to distribute the total number of
  processes to the YASP grid by passing a vector of type
  Dune::FieldVector<int,dim> to the constructor.
*/
template<int dim, class iTupel>
class Ax1YaspPartition : public Dune::YLoadBalance<dim>
{

public:
  //constructor:
  Ax1YaspPartition( const iTupel& yasppartitions_ )
    : yasppartitions( yasppartitions_ )
  {
  }

  void loadbalance (const iTupel& size, int P, iTupel& dims) const
  {
    dims = yasppartitions;
    debug_info << "Partitioning grid, size = " << size << ", dims = " << dims << std::endl;
  }

private:
  const iTupel& yasppartitions;

};


#endif /* DUNE_AX1_YASPGRID_LOADBALANCER_HH */
