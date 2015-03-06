#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <typeinfo>

#include <dune/common/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfwriter.hh>
#include <dune/grid/geometrygrid/grid.hh>

#include <dune/ax1/common/constants.hh>

#include <dune/ax1/common/ax1_yaspgrid_loadbalancer.hh>
#include <dune/ax1/common/ax1_subgridview.hh>
#include <dune/ax1/common/ax1_subgrid_hostelement_mapper.hh>
#include <dune/ax1/common/ax1_parallelhelper.hh>
//#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/hdf5_tools.hh>
#include <dune/ax1/channels/channel_builder.hh>




//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is acme1MD. Extrem cremig." << std::endl;
    else
	  {
      if(helper.rank()==0)
        std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
	  }
    
    if (argc<4)
	  {
      if(helper.rank()==0)
        std::cout << "usage: ./acme1MD <level> <dtstart> <tend> [config-file]" << std::endl;
      return 1;
	  }
    
    int level;
    sscanf(argv[1],"%d",&level);
    
    double dtstart;
    sscanf(argv[2],"%lg",&dtstart);
    
    double tend;
    sscanf(argv[3],"%lg",&tend);
    
    std::string configFileName;
    if (argc>=5)
    {
      configFileName = argv[4];
    } else {
      configFileName = "acme1MD.config";
    }
    debug_verb << "Using config file " << configFileName << std::endl;

    Dune::Timer timer;

    //GridType hostGrid(coords);

    const int dim = 2;
    Dune::FieldVector<double,dim> L;
    L[0] = 100;
    L[1] = 50;
    Dune::FieldVector<int,dim> N;
    N[0] = 10;
    N[1] = 5;
    Dune::FieldVector<bool,dim> periodic (false);
    int overlap = 1;

    // Set up custom load balancing
    int px = helper.size();
    int py = 1; // Do not partition y-direction!

    typedef Ax1YaspPartition<2,Dune::FieldVector<int,2> > YP;
    if (px*py == 0 || px*py !=helper.size())
    {
      DUNE_THROW(Dune::Exception, "px*py = 0 or != np");
    }
    if(px > N[0])
    {
      DUNE_THROW(Dune::Exception, "Cannot partition grid onto " << px << " processes, as there are only "
          << N[0] << " element stripes in y-direction!");
    }

    Dune::FieldVector<int,2> yasppartitions;
    yasppartitions[0] = px;
    yasppartitions[1] = py;
    YP* yp = new YP(yasppartitions);
    if( helper.rank() == 0 )
    {
      debug_info << "Partitioning of YASP: " << yasppartitions << std::endl;
    }
    Dune::YaspGrid<2> yaspGrid(helper.getCommunicator(),L,N,periodic,overlap,yp);

    typedef Dune::YaspGrid<2>::LeafGridView GV;
    GV gv = yaspGrid.leafGridView();

    // MPI test stuff
    typedef std::tuple<std::string,int,double,float> MyTuple;
    std::stringstream str;
    str << "Hello world, #" << helper.rank() << " speaking!";
    MyTuple tuple(str.str(), 12, 3.14, helper.rank());

    std::cout << "int: " << sizeof(int) << ", double: " << sizeof(double) << ", float: " << sizeof(float) << std::endl;
    std::cout << "string: " << sizeof(std::string) << std::endl;
    std::cout << sizeof(str.str()) << std::endl;

    MyTuple recv_tuple;

    if(helper.rank() == 0)
    {
      int val = MPI_Send(&tuple, 1, Dune::MPITraits<MyTuple>::getType(), 1, 12, gv.comm());

      std::cout << "Sent tuple: " << std::get<2>(tuple) << " | " << val << std::endl;
    } else {
      MPI_Status status;
      int val = MPI_Recv(&recv_tuple, 1, Dune::MPITraits<MyTuple>::getType(), 0, 12, gv.comm(), &status);

      std::cout << "Received tuple: " << std::get<2>(recv_tuple) << " | " << val  << std::endl;
    }


      
    double elapsed = timer.stop();
    debug_info << "Time elapsed: " << elapsed << " s" << std::endl;

  } catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 0;
  }

}
