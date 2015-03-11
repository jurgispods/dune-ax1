 #ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

// Preprocessor defines
// YaspGrid=1, UG=2
#define USE_GRID 1
// Sequential=0, Parallel=1
#ifndef AX1_PARALLEL
#define USE_PARALLEL 0
#else
#define USE_PARALLEL 1
#define USE_OVERLAP 1 // Does not converge without overlap, why?
#endif

// Use LocalBasisCache in local operator
#define USE_CACHE 1

//#define DUNE_ISTL_WITH_CHECKING

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
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG 
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
#if HAVE_DUNE_MULTIDOMAINGRID
#include <dune/grid/multidomaingrid/multidomaingrid.hh>
#endif

#include <dune/pdelab/common/logtag.hh>

// Include constants.hh here in order to make the DOF-ordering-hack in dune-multidomain work
// TODO Use PDELab native ordering once it was ported to the new branch
#include <dune/ax1/common/constants.hh>

#if HAVE_DUNE_MULTIDOMAIN
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinggridfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinglocalfunctionspace.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#endif

#include <dune/ax1/common/ax1_gridgenerator.hh>
#include <dune/ax1/common/ax1_gridtools.hh>
#include <dune/ax1/common/ax1_gridvector.hh>
#include <dune/ax1/common/ax1_tensorgrid_transformation.hh>
#include <dune/ax1/common/ax1_yaspgrid_loadbalancer.hh>
#include <dune/ax1/common/ax1_subgridview.hh>
#include <dune/ax1/common/ax1_subgrid_hostelement_mapper.hh>
#include <dune/ax1/common/gnuplot_tools.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_factory.hh>
//#include <dune/ax1/acme2_cyl/common/acme2_cyl_geometrycheck.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_physics.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_setup.hh>
#include <dune/ax1/channels/channel_builder.hh>
#include <dune/ax1/common/ax1_lfs_tools.hh>
#include <dune/ax1/common/ax1_parallelhelper.hh>


template<typename Element>
void test(const Element& e0)
{
  debug_jochen << "testHallo1 @" << e0.geometry().center() << std::endl;
  if(e0.ileafbegin()->neighbor())
  {
    auto ep1 = e0.ileafbegin()->outside();
    //const Element& e1 = e0;
    debug_jochen << "testHallo2 @" << ep1->geometry().center() << std::endl;
    {
      debug_jochen << "call1" << std::endl;
      auto iit_test = ep1->ileafbegin();
      debug_jochen << "OK" << std::endl;
      debug_jochen << "call2" << std::endl;
      auto iit_test2 = ep1->ileafbegin();
      debug_jochen << "OK" << std::endl;
//      debug_jochen << "call3" << std::endl;
//      auto iit_test3 = e1.ileafbegin();
//      debug_jochen << "OK" << std::endl;
//      debug_jochen << "call4" << std::endl;
//      auto iit_test4 = e1.ileafbegin();
//      debug_jochen << "OK" << std::endl;
    }
  }
}

template<typename Grid, typename Config, typename GV, typename SubGV>
void run(Grid& grid, Acme2CylParameters& params, Config& config, GV& gv, double dtstart, double tend,
    ElementSubdomainMapper& elemSubdomainMapper, Ax1ElementGroupMapper& elemGroupMapper, SubGV& elecGV, SubGV& membGV)
{
  debug_verb << "Configuration '" << config.getConfigName() << "' loaded." << std::endl;

  typedef typename Acme2CylFactory<GV,double,Config,SubGV>::PHYSICS PHYSICS;
  PHYSICS physics = Acme2CylFactory<GV,double,Config,SubGV>::setup(config, gv, elecGV, membGV, params, dtstart,
      elemSubdomainMapper, elemGroupMapper);
  //debug_info << Tools::getTypeName(physics) << std::endl;

  //physics.gridInfo(subGV_Inside, subGV_Outside);

  // Check parameters used
  debug_verb << "==============================================" << std::endl;
  debug_verb << "LENGTH_SCALE = " << physics.getLengthScale() << std::endl;
  debug_verb << "TIME_SCALE   = " << physics.getTimeScale() << std::endl;
  for(int i=0; i<NUMBER_OF_SPECIES; ++i)
  {
    debug_verb << "con_diffWater" << ION_NAMES[i] << " = " << physics.getElectrolyte().getDiffConst(i) << std::endl;
  }
  debug_verb << "electrolyte temperature: " << physics.getElectrolyte().getTemperature() << std::endl;
  debug_verb << "electrolyte stdCon: " << physics.getElectrolyte().getStdCon() << std::endl;
  debug_verb << "poissonConstant: " << physics.getPoissonConstant() << std::endl;
  debug_verb << "==============================================" << std::endl;

  // Happy testing starts here!

  int i = 0;
  //for(typename PHYSICS::ElementIterator_All eit = gv.template begin<0>();
  //      eit != gv.template end<0>(); ++eit)
  for(typename PHYSICS::SubDomainElementIterator_All meit = membGV.template begin<0>();
        meit != membGV.template end<0>(); ++meit)
  {
    i++;

    if(physics.isMembrane(*meit))
    {
      const typename PHYSICS::Element& e0 = membGV.grid().multiDomainEntity(*meit);

      // WARUM IST DAS HIER :

      debug_jochen << "Hallo1 @" << e0.geometry().center() << std::endl;
      if(e0.ileafbegin()->neighbor())
      {
        typename PHYSICS::ElementPointer ep1 = e0.ileafbegin()->outside();
        debug_jochen << "Hallo2 @" << ep1->geometry().center() << std::endl;
        {
          debug_jochen << "call1" << std::endl;
          auto iit_test = ep1->ileafbegin();
          debug_jochen << "OK" << std::endl;
          debug_jochen << "call2" << std::endl;
          auto iit_test2 = ep1->ileafbegin();
          debug_jochen << "OK" << std::endl;
//          debug_jochen << "call3" << std::endl;
//          auto iit_test3 = e1.ileafbegin();
//          debug_jochen << "OK" << std::endl;
//          debug_jochen << "call4" << std::endl;
//          auto iit_test4 = e1.ileafbegin();
//          debug_jochen << "OK" << std::endl;
        }
        //sleep(0.1);
      }

      debug_jochen << "-----------" << std::endl;

      // ANDERS ALS DAS HIER !=!=!??!??!?!?!?!?!??)(/&()%&6&&6&&(&%%&E$§
      //test(membGV.grid().multiDomainEntity(*meit));
      test(e0);

      debug_jochen << "===========" << std::endl;
    }

  }



  int count = 0;
  // Loop over membrane elements and find primary intersection
  typename PHYSICS::ElementIntersectionIterator iit_primary = gv.template begin<0>()->ileafbegin();
  for(typename PHYSICS::SubDomainElementIterator_All meit = membGV.template begin<0>();
      meit != membGV.template end<0>(); ++meit)
  {
    const typename PHYSICS::ElementIntersectionIterator& iit = physics.getNextMembraneInterface(*meit);

    count++;
    if(count == 25)
    {
      //iit_primary = iit;
      //iit_primary = membGV.grid().multiDomainEntity(*meit).ileafbegin();
    }
  }

  debug_jochen << "Found intersection @" << iit_primary->geometry().center() << std::endl;
  debug_jochen << " with inside @" << iit_primary->inside()->geometry().center() << std::endl;
  if(iit_primary->neighbor())
    debug_jochen << " and outside @" << iit_primary->outside()->geometry().center() << std::endl;

  // =====================================================================================================

}

//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    std::stringstream logtag;

    if(Dune::MPIHelper::isFake)
      debug_info<< "This is acme2_cyl. Extrem cremig." << std::endl;
    else
	{
      logtag << "[p" << helper.rank() << "] ";

      debug_verb.setLogTag(logtag.str());
      debug_info.setLogTag(logtag.str());
      debug_warn.setLogTag(logtag.str());
      debug_minimal.setLogTag(logtag.str());
      debug_jochen.setLogTag(logtag.str());

      if(helper.rank()==0)
      {
        debug_info << "parallel run on " << helper.size() << " process(es)" << std::endl;
      } else {
        // Decativate debug output on non-root processes
        debug_verb.push(false);
        debug_info.push(false);
        debug_warn.push(false);
        debug_minimal.push(false);
        debug_jochen.push(false);
      }
	}
    
    if (argc<4)
	  {
      if(helper.rank()==0)
        debug_info << "usage: ./acme2_cyl <level> <dtstart> <tend> [config-file]" << std::endl;
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
      if(USE_CYLINDER_COORDINATES)
        configFileName = "acme2_cyl.config";
      else
        configFileName = "acme2_2d.config";
    }
    debug_verb << "Using config file " << configFileName << std::endl;

    Dune::Timer timer;

    // Read config file
    Acme2CylParameters params;
    Dune::ParameterTreeParser::readINITree(configFileName, params, true);
    params.init(configFileName, argc, argv);

    if(! params.doOutputRootNodeOnly() && helper.rank() > 0)
    {
      // Activate debug output on non-root processes
      debug_verb.pop();
      debug_info.pop();
      debug_warn.pop();
      debug_minimal.pop();
      debug_jochen.pop();
    }

    // ========================= Grid generation ==============================
    Ax1GridGenerator gridGenerator(params, level);
    GridVector<double> x;
    CylinderGridVector<double> y(params.dX());
    gridGenerator.generateTensorGrid(x, y);

    // Now copy x,y gridvectors in to a std::vector which can be passed to the parameter class
    std::vector<double> x_std(x);
    std::vector<double> y_std(y);
    params.setX(x_std);
    params.setY(y_std);

    debug_jochen << "===== x: ===================" << std::endl;
    Output::printVector(x_std);
    debug_jochen << "============================" << std::endl;
    debug_jochen << "===== y: ===================" << std::endl;
    Output::printVector(y_std);
    debug_jochen << "============================" << std::endl;

    debug_jochen << "x.size() = " << x.size() << std::endl;
    debug_jochen << "y.size() = " << y.size() << std::endl;
    debug_jochen << "=> Will give grid of " << ((x.size()-1) * (y.size()-1)) << " cells!" << std::endl;

    if(params.doEquilibration() && x.size() > 2 && !params.doForceEquilibration())
      DUNE_THROW(Dune::Exception, std::string("The new data transfer only works with grids having exactly 1 element ")
         + std::string("in x-direction! Please consider adjusting your grid parameters."));


    // YaspGrid + GeometryGrid stuff
#if USE_GRID==1
    //GridType hostGrid(coords);

    const int dim = 2;
    Dune::FieldVector<double,dim> L;
    L[0] = params.xMax();
    L[1] = params.yMax();
    Dune::FieldVector<int,dim> N;
    N[0] = x.size()-1;
    N[1] = y.size()-1;
    Dune::FieldVector<bool,dim> periodic (false);
#if USE_OVERLAP == 1
    int overlap = 1;
#else
    int overlap = 0;
#endif

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
    //yaspGrid.globalRefine(level);

    typedef Ax1TensorGridTransformation<double> CoordFunction;
    CoordFunction coordFunction(params.xMax(),params.yMax(),x,y);
    typedef Dune::GeometryGrid<BaseGrid,CoordFunction> HostGrid;
    HostGrid hostGrid(yaspGrid,coordFunction);

#elif USE_GRID==2 // UG case
    //BaseGrid::setDefaultHeapSize(2000);
    typedef BaseGrid HostGrid;

    //get GridFactory instance
    Dune::GridFactory<BaseGrid> factory;

    //get geometry type
    const Dune::GeometryType gt(Dune::GeometryType::cube,2);

    Ax1GridTools::prepareGrid(factory, x, y, gt);

    //make grid
    Dune::shared_ptr<BaseGrid> pBaseGrid(factory.createGrid());
    HostGrid& hostGrid = *pBaseGrid;

#else
#warning "Selected grid could not be found!"
#endif
    debug_info << "Grid creation finished, time elapsed: " << timer.elapsed() << "s" << std::endl;
    // ############## END grid generation ###############

    // ################### Test grid decomposition ###############################
    debug_jochen << "Leaf GV size before grid.loadBalance(): " << hostGrid.leafGridView().size(0) << std::endl;
    hostGrid.loadBalance();
    debug_jochen << "Leaf GV size after grid.loadBalance(): " << hostGrid.leafGridView().size(0) << std::endl;

    typename HostGrid::Partition<Dune::PartitionIteratorType::Interior_Partition>::LeafGridView interiorGV
     = hostGrid.leafGridView<Dune::PartitionIteratorType::Interior_Partition>();
    debug_jochen << "Leaf interior GV size after grid.loadBalance(): " << interiorGV.size(0) << std::endl;

    typedef HostGrid::LeafGridView BaseGV;
    BaseGV basegv = hostGrid.leafGridView();

    debug_jochen << "Processor #" << helper.rank() << " has " << basegv.size(0)
        << "elements." << std::endl;
    debug_jochen << "Processor #" << helper.rank() << " has " << basegv.overlapSize(0)
        << " overlap element layers" << std::endl;
    debug_jochen << "Processor #" << helper.rank() << " has " << basegv.ghostSize(0)
        << " ghost element layers" << std::endl;

    typedef BaseGV::Codim<0>::EntityPointer BElementPointer;
    typedef Dune::SingleCodimSingleGeomTypeMapper<BaseGV, 0> BElementMapper;

    // Custom mapper to map element indices to subdomain indices
    BElementMapper belementMapper(basegv);

//      typedef BaseGV::Codim<0>::Partition<
//          Dune::PartitionIteratorType::Interior_Partition>::Iterator IElementIterator;
//      // print info about interior elements
//      for(IElementIterator ep = basegv.begin<0,Dune::PartitionIteratorType::Interior_Partition>();
//          ep != basegv.end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++ep)
//      {
//        debug_jochen << "Processor #" << helper.rank() << " has element #"
//            << belementMapper.map(*ep) << " @"
//            << ep->geometry().center() << " partition: "
//            << ep->partitionType() << std::endl;
//      }
//
//      typedef BaseGV::Codim<0>::Partition<
//          Dune::PartitionIteratorType::Overlap_Partition>::Iterator OElementIterator;
//      // print info about overlap elements
//      for(OElementIterator ep = basegv.begin<0,Dune::PartitionIteratorType::Overlap_Partition>();
//          ep != basegv.end<0,Dune::PartitionIteratorType::Overlap_Partition>(); ++ep)
//      {
//        debug_jochen << "Processor #" << helper.rank() << " has element #"
//            << belementMapper.map(*ep) << " @"
//            << ep->geometry().center() << " partition: "
//            << ep->partitionType() << std::endl;
//      }
//
//      typedef BaseGV::Codim<0>::Partition<
//          Dune::PartitionIteratorType::Ghost_Partition>::Iterator GElementIterator;
//      // print info about ghost elements
//      for(GElementIterator ep = basegv.begin<0,Dune::PartitionIteratorType::Ghost_Partition>();
//          ep != basegv.end<0,Dune::PartitionIteratorType::Ghost_Partition>(); ++ep)
//      {
//        debug_jochen << "Processor #" << helper.rank() << " has element #"
//            << belementMapper.map(*ep) << " @"
//            << ep->geometry().center() << " partition: "
//            << ep->partitionType() << std::endl;
//      }

    // ============== Multidomaingrid stuff ================================================================
    //typedef Dune::MultiDomainGrid<HostGrid, Dune::mdgrid::FewSubDomainsTraits<HostGrid::dimension,2> > Grid;
    typedef Dune::MultiDomainGrid<HostGrid,
        Dune::mdgrid::FewSubDomainsTraits<HostGrid::dimension,2,Dune::mdgrid::CellAndVertexCodims> > Grid;
    Grid grid(hostGrid,false);

    typedef Grid::LeafGridView GV;
    GV gv = grid.leafGridView();

    typedef GV::Codim<0>::Iterator ElementIterator;
    typedef GV::Codim<0>::EntityPointer ElementPointer;
    typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

    // Custom mapper to map element indices to subdomain indices
    Ax1ElementGroupMapper elemGroupMapper;
    ElementSubdomainMapper elemSubdomainMapper;
    ElementMapper elementMapper(gv);

    // Grid parameters
    double yMemb = params.yMemb();
    double dMemb = params.dMemb();
    int numMembranes = params.nMembranes();
    double cellWidth = params.getCellWidth();

    double xNode = params.xNode();
    double nodeWidth = params.dNode();

    int nElementsThisProcess = 0;
    // Tag elements by type
    for(ElementIterator ep = gv.begin<0>(); ep != gv.end<0>(); ++ep)
    {
      // Only count interior elements
      if(ep->partitionType() == Dune::PartitionType::InteriorEntity)
        nElementsThisProcess++;

      int elemIndex = elementMapper.map(*ep);

      // default: intracellular
      int subdomainIndex = CYTOSOL;

      typename ElementIterator::Entity::Geometry::GlobalCoordinate center = ep->geometry().center();

      typename ElementIterator::Entity::Geometry::GlobalCoordinate diagonal
        = ep->geometry().corner(ep->geometry().corners()-1);
      diagonal -= ep->geometry().corner(0);

      double h_y = std::abs(diagonal[1]);

      //assert(std::abs(h_y - d) < 1e-8);

      // Assume there is a single layer of elements for each membrane in the following
      if(params.useMembrane())
      {
        // We are on the membrane{
        if(center[1] > yMemb and center[1] < yMemb+dMemb)
        {
          subdomainIndex = MEMBRANE;
        }

        if(numMembranes == 1)
        {
          // There is only one membrane => the rest is outside of the cell
          if(center[1] > yMemb+dMemb)
          {
            subdomainIndex = ES;
          }
        }
        // In case of a second membrane layer
        if(numMembranes == 2)
        {
          // Tag second membrane
          if(center[1] > yMemb+dMemb+cellWidth and center[1] < yMemb+dMemb+cellWidth+dMemb)
          {
            subdomainIndex = MEMBRANE;
          }
          // Tag this cell as extracellular if it is outside the two membrane layers
          if(center[1] < yMemb or center[1] > yMemb+dMemb+cellWidth+dMemb)
          {
            subdomainIndex = ES;
          }
        }
      } else {
        // If no membrane is present, everything is extracellular space
        subdomainIndex = ES;
      }
      //debug_jochen << "Tagging element @ " << center[1] << " as " << subdomainIndex << std::endl;

      // Store subdomain of this element [M->G]
      elemSubdomainMapper.setSubdomain(elemIndex, subdomainIndex);

      if(subdomainIndex != MEMBRANE)
      {
        // not membrane: subdomainIndex == groupIndex
        elemGroupMapper.setGroup(elemIndex, subdomainIndex);
      } else {
        // membrane: check if node of Ranvier or not
        if(center[0] > xNode and center[0] < xNode+nodeWidth)
        {
          elemGroupMapper.setGroup(elemIndex, 2);
        } else {
          elemGroupMapper.setGroup(elemIndex, 3);
        }

      }
    }
    params.setNElementsThisProcess(nElementsThisProcess);

    typedef Grid::SubDomainGrid SubDomainGrid;
    SubDomainGrid& elecGrid = grid.subDomain(0);
    SubDomainGrid& membGrid = grid.subDomain(1);
    typedef Grid::ctype ctype;

    typedef SubDomainGrid::LeafGridView SDGV;
    SDGV elecGV = elecGrid.leafGridView();
    SDGV membGV = membGrid.leafGridView();

    grid.startSubDomainMarking();
    // Only add interior elements, multidomaingrid handles overlap automatically!
    typedef GV::Codim<0>::Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator ElementIterator_Interior;
    for (ElementIterator_Interior eit = gv.begin<0,Dune::PartitionIteratorType::Interior_Partition>();
        eit != gv.end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++eit)
    {
      int elemIndex = elementMapper.map(*eit);
      int subdomainIndex = elemSubdomainMapper.map(elemIndex);

      // Init subgrids
      switch(subdomainIndex)
      {
        case CYTOSOL:
          grid.addToSubDomain(0,*eit);
          break;
        case ES:
          grid.addToSubDomain(0,*eit);
          break;
        case MEMBRANE:
          grid.addToSubDomain(1,*eit);
          break;
        default:
          DUNE_THROW(Dune::Exception, "Element could not be assigned to a pre-defined subdomain!");
      }

    }
    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();
    // ===============================================================================================


    // ================ Cylindric coordinates check ==================================================
//      Acme2CylGeometryTools cylTools(params);
//
//      // Check volume/surface
//      double vol = 0.0;
//      cylTools.checkCylinderVolume(gv, vol, 4);
//      cylTools.checkCylinderSurface(gv, vol, 4);
    // ===============================================================================================

    bool configFound = false;
    std::string configName = params.getConfigName();

    if(configName == "default")
    {
      DefaultConfiguration<double> config;
      run(grid,params,config,gv,dtstart,tend,
          elemSubdomainMapper,elemGroupMapper,elecGV,membGV);
      configFound = true;
    }
    if(configName == "ES")
    {
      ESConfiguration<double> config;
      run(grid,params,config,gv,dtstart,tend,
          elemSubdomainMapper,elemGroupMapper,elecGV,membGV);
      configFound = true;
    }
    // Other configurations disabled to speed up compilation
    /*
    if(configName == "step")
    {
      StepConfiguration<double> config;
      run(grid,params,config,gv,dtstart,tend,
          elemSubdomainMapper,elemGroupMapper,elecGV,membGV);
      configFound = true;
    }
    if(configName == "test_scales")
    {
      TestScalesConfiguration<double> config;
      run(grid,params,config,gv,dtstart,tend,
          elemSubdomainMapper,elemGroupMapper,elecGV,membGV);
      configFound = true;
    }
    */


    if(! configFound)
    {
      DUNE_THROW(Dune::Exception, "No configuration named '" << configName << "' could be found!");
    }

    double elapsed = timer.stop();
    debug_info << "Time elapsed: " << elapsed << " s" << std::endl;

  } catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    throw e; // Throw exception in order to get a full stacktrace for debugging
  }

}
