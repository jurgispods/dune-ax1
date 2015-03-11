#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

// Preprocessor defines
// YaspGrid=1, UG=2
#define USE_GRID 2
// Sequential=0, Parallel=1
#define USE_PARALLEL 0

// Use LocalBasisCache in local operator
#define USE_CACHE 1

#define DUNE_ISTL_WITH_CHECKING

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

#include <dune/grid/io/file/dgfparser/dgfwriter.hh>
#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>
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

  // Set up and start simulation
  Acme2CylSetup<Grid,GV,PHYSICS,SubGV> acme2_cylSetup(grid, gv, physics, elecGV, membGV);
  acme2_cylSetup.setup(dtstart, tend);
}

//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is acme2_cyl. Extrem cremig." << std::endl;
    else
	  {
      if(helper.rank()==0)
        std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
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

    // sequential version
    if (1 && helper.size()==1)
    {
      // Read config file
      Acme2CylParameters params;
      Dune::ParameterTreeParser::readINITree(configFileName, params, true);
      params.init(configFileName, argc, argv);
      
#if USE_GRID==1
      typedef BaseGrid GridType;
      //GridType baseGrid(coords);
      
      const int dim = 2;
      Dune::FieldVector<double,2*1> L;
      L[0] = params.xMax();
      L[1] = params.yMax();
      Dune::FieldVector<int,dim> N (1);
      Dune::FieldVector<bool,dim> periodic (false);
      int overlap = 1;
      GridType baseGrid (L,N,periodic,overlap);
      baseGrid.globalRefine(level);
#elif USE_GRID==2 // UG case
      typedef BaseGrid GridType;

      //BaseGrid::setDefaultHeapSize(2000);

      //get GridFactory instance
      Dune::GridFactory<BaseGrid> factory;

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

      //get geometry type
      const Dune::GeometryType gt(Dune::GeometryType::cube,2);

      Ax1GridTools::prepareGrid(factory, x, y, gt);

      //make grid
      debug_jochen << "Maken GridPointer." << std::endl;
      Dune::shared_ptr<BaseGrid> pBaseGrid(factory.createGrid());
      debug_jochen << "Maken Grid." << std::endl;
      BaseGrid& baseGrid = *pBaseGrid;

#else
#warning "Selected grid could not be found!"
#endif
      debug_info << "Grid creation finished, time elapsed: " << timer.elapsed() << "s" << std::endl;
      // ############## END grid generation ###############
      
      
      // ============== Multidomaingrid stuff ================================================================
      typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,2> > Grid;
      Grid grid(baseGrid,false);

      typedef Grid::LeafGridView GV;
      GV gv = grid.leafView();

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
          // membrane: check which membrane group this element belongs to
          if(params.membrane.getSubKeys().size() < 1)
            DUNE_THROW(Dune::Exception, "Coulrd not find a matching membrane group in config file [membrane] section for element @"
                << center << "!");
          // There is only one membrane group (besides solution_in, solution_out), assign it to element
          if(params.membrane.getSubKeys().size() == 1)
          {
            elemGroupMapper.setGroup(elemIndex, 2);
          } else {
            // TODO Implement group-wise definition of location along the membrane where the group is defined!
            // 1) Add config file parameters 'start', 'width', 'stride' to every group
            //    Also: Add group-specific permittivity value there and hand it over to membrane later on!
            // 2) Check here into which interval the element belongs
            // 3) Check if interval definition is consistent, i.e. if it is a partition of the whole membrane subdomain!
            // 4) Determine group index as group index in config file +2 (solution_in, solution_ex are always there)

            assert(params.membrane.getSubKeys().size() == 2);
            // For now: Use hardcoded method with only two groups and global parameters xNode and nodeWidth:
            if(center[0] > xNode and center[0] < xNode+nodeWidth)
            {
              // Node of ranvier
              elemGroupMapper.setGroup(elemIndex, 2);
            } else {
              // Myelin
              elemGroupMapper.setGroup(elemIndex, 3);
            }
          }


        }
      }
      params.setNElementsThisProcess(nElementsThisProcess);

      typedef Grid::SubDomainGrid SubDomainGrid;
      SubDomainGrid& elecGrid = grid.subDomain(0);
      SubDomainGrid& membGrid = grid.subDomain(1);
      typedef Grid::ctype ctype;

      typedef SubDomainGrid::LeafGridView SDGV;
      SDGV elecGV = elecGrid.leafView();
      SDGV membGV = membGrid.leafView();

      grid.startSubDomainMarking();
      for (GV::Codim<0>::Iterator eit = gv.begin<0>(); eit != gv.end<0>(); ++eit)
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
      /*
      Acme2CylGeometryTools cylTools(params);

      // Check volume/surface
      double vol = 0.0;
      cylTools.checkCylinderVolume(gv, vol, 4);
      cylTools.checkCylinderSurface(gv, vol, 4);
      */

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
      // Other configurations disabled to speed up compilation
      /*
      if(configName == "ES")
      {
        ESConfiguration<double> config;
        run(grid,params,config,gv,dtstart,tend,
            elemSubdomainMapper,elemGroupMapper,elecGV,membGV);
        configFound = true;
      }
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
    }
  } catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    throw e; // Throw exception in order to get a full stacktrace for debugging
  }

}
