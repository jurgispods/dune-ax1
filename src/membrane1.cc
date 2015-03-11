// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file

    \brief Solve elliptic problem in unconstrained spaces with conforming finite elements
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>

#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#include <dune/grid/onedgrid.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_ALUGRID
#include<dune/grid/alugrid.hh>
#include<dune/grid/io/file/dgfparser/dgfalu.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif

#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/p12dfem.hh>
#include<dune/pdelab/finiteelementmap/pk2dfem.hh>
#include<dune/pdelab/finiteelementmap/pk3dfem.hh>
#include<dune/pdelab/finiteelementmap/q12dfem.hh>
#include<dune/pdelab/finiteelementmap/q22dfem.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/finiteelementmap/p1fem.hh>
#include<dune/pdelab/finiteelementmap/rannacher_turek2dfem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/constraints.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>


#include<dune/ax1/membrane1/membrane1_Pk.hh>

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
	std::cout<< "This is Membrane1. Geil." << std::endl;
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
	  {
		if(helper.rank()==0)
		  std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
	  }

	if (argc!=2)
	  {
		if(helper.rank()==0)
		  std::cout << "usage: ./membrane1 <level>" << std::endl;
		return 1;
	  }


	int level = 0;
	//sscanf(argv[1],"%d",&level);
  
  
	const int dim = 1;

	// sequential version
    if (helper.size()==1)
    {
      // grid generation #####################################################################
      
      /*
      Dune::FieldVector<double,dim> L(1.0);
      Dune::FieldVector<int,dim> N(1);
      Dune::FieldVector<bool,dim> periodic(false);
      int overlap=0;
      Dune::YaspGrid<dim> grid(L,N,periodic,overlap);
      */
      
      //std::vector<double> coords = {-5., -4., -3., -2., -1., 0., 1., 2., 3., 4., 5.};
      
      std::vector<double> coords;
      for (int i=0; i<11; i++)
      {
        coords.push_back((double) i-5);
        //std::cout << coords[i] << std::endl;
      }
      
      Dune::OneDGrid grid(coords);

      grid.globalRefine(level);
      //typedef Dune::YaspGrid<dim>::LeafGridView GV;
      typedef Dune::OneDGrid::LeafGridView GV;
      const GV& gv=grid.leafView();
      membrane1_Pk(gv);
    }
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
	return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
	return 1;
  }
}
