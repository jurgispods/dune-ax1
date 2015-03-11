 #ifdef HAVE_CONFIG_H
#include "config.h"     
#endif


// Preprocessor defines
// YaspGrid=1, UG=2
#define USE_GRID 1
// Sequential=0, Parallel=1
#ifndef AX1_PARALLEL
#define USE_PARALLEL (0)
#else
#define USE_PARALLEL (1)
#define USE_OVERLAP (1) // Does not converge without overlap, why?
#define OVERLAP_SIZE (1)
#endif

// Use LocalBasisCache in local operator
//#define USE_CACHE 1

//#define DUNE_ISTL_WITH_CHECKING

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <typeinfo>

#include <dune/common/tuples.hh>
#include <dune/common/parallel/mpihelper.hh>
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
#include <dune/ax1/acme2_cyl/common/laplace_setup.hh>
#include <dune/ax1/acme2_cyl/common/eval_grids_setup.hh>
#include <dune/ax1/channels/channel_builder.hh>
#include <dune/ax1/common/ax1_lfs_tools.hh>
#include <dune/ax1/common/ax1_parallelhelper.hh>

template<typename Grid, typename Config, typename GV, typename SubGV>
void run(std::vector<Dune::shared_ptr<Grid> >& grids, Acme2CylParameters& params, Config& config, GV& gv, double dtstart, double tend,
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
  for(int i=0; i<Config::numberOfSpecies; ++i)
  {
    debug_verb << "con_diffWater" << ION_NAMES[i] << " = " << physics.getElectrolyte().getDiffConst(i) << std::endl;
  }
  debug_verb << "electrolyte temperature: " << physics.getElectrolyte().getTemperature() << std::endl;
  debug_verb << "electrolyte stdCon: " << physics.getElectrolyte().getStdCon() << std::endl;
  debug_verb << "poissonConstant: " << physics.getPoissonConstant() << std::endl;
  debug_verb << "==============================================" << std::endl;

  // Set grid evaluation
  EvalGridsSetup<Grid,GV,PHYSICS,SubGV> evalGridsSetup(grids, gv, physics, elecGV, membGV);
  evalGridsSetup.setup(dtstart, tend);
}

//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  Dune::Timer timer;

  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    std::stringstream logtag;

    if(helper.rank() == 0)
    {
      debug_info<< "This is acme2_cyl. Extrem cremig." << std::endl;
      debug_info << "Compiled with flags:" << std::endl;
      debug_info << " USE_GRID: " << USE_GRID << std::endl;
      debug_info << " USE_PARALLEL: " << USE_PARALLEL << std::endl;
      debug_info << " USE_OVERLAP: " << USE_OVERLAP << std::endl;
      debug_info << " USE_CYLINDER_COORDINATES: " << USE_CYLINDER_COORDINATES << std::endl;
      debug_info << " DO_MULTIDOMAIN_BLOCKING: " << DO_MULTIDOMAIN_BLOCKING << std::endl;
      debug_info << " AX1_BLOCKSIZE: " << AX1_BLOCKSIZE << std::endl;
    }

    if(Dune::MPIHelper::isFake)
    {
      debug_info << "This is sequential run!" << std::endl;
    } else {
      logtag << "[p" << helper.rank() << "] ";

      debug_verb.setLogTag(logtag.str());
      debug_info.setLogTag(logtag.str());
      debug_warn.setLogTag(logtag.str());
      debug_minimal.setLogTag(logtag.str());
      debug_jochen.setLogTag(logtag.str());

      if(helper.rank() == 0)
      {
        debug_info << "parallel run on " << helper.size() << " process(es)" << std::endl;
      } else {
        // Deactivate debug output on non-root processes
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

    // Switch from output rank 0 to another rank, if desired
    if(params.doOutputRootNodeOnly() && helper.size() > 1)
    {
      int rootOutputNode = params.general.get("rootOutputNode", 0);

      // Deactivate rank 0
      if(helper.rank() == 0 && helper.rank() != rootOutputNode)
      {
        debug_verb.push(false);
        debug_info.push(false);
        debug_warn.push(false);
        debug_minimal.push(false);
        debug_jochen.push(false);
      }
      // Activate new root output node
      if(helper.rank() != 0 && helper.rank() == rootOutputNode)
      {
        debug_verb.pop();
        debug_info.pop();
        debug_warn.pop();
        debug_minimal.pop();
        debug_jochen.pop();
      }
    }

    // ========================= Grid generation ==============================
    Ax1GridGenerator gridGenerator(params, level);
    GridVector<double> x_orig;
    CylinderGridVector<double> y_orig(params.dX());
    std::vector<typename Ax1GridGenerator::MembraneGroupTuple> membGroups;
    gridGenerator.generateTensorGrid(x_orig, y_orig, membGroups);
    params.setMembraneGroups(membGroups);

    if(params.doEquilibration() && x_orig.size() > 2 && !params.doForceEquilibration())
      DUNE_THROW(Dune::Exception, std::string("The new data transfer only works with grids having exactly 1 element ")
         + std::string("in x-direction! Please consider adjusting your grid parameters."));

    typedef Ax1GridGenerator::Grid Grid;

    const int nGrids = params.general.get("nGrids",1);

    std::vector<int> defaultLevels(nGrids,0);
    std::vector<int> levels = params.general.get("gridLevels",defaultLevels);

    typedef Dune::shared_ptr<Grid> GridP;
    std::vector<GridP> grids(nGrids);

    Ax1ElementGroupMapper elemGroupMapper_last;
    ElementSubdomainMapper elemSubdomainMapper_last;

    GridVector<double> x(x_orig);
    CylinderGridVector<double> y(y_orig);
    for(int i=0; i<nGrids; i++)
    {
      // Copy original vector, otherwise we add up global refines of successive iterations
      x = x_orig;
      y = y_orig;

      // New method globalRefine which basically does bisection, but excludes membrane elements in y-direction
      if(params.general.get("doAdditionalGlobalRefine",false))
      {
        gridGenerator.globalRefine(x,y,levels[i]);
      }

      // Now copy x,y gridvectors in to a std::vector which can be passed to the parameter class
      std::vector<double> x_std(x);
      std::vector<double> y_std(y);
      params.setX(x_std);
      params.setY(y_std);

      debug_jochen << "===== x: ===================" << std::endl;
      Output::printVector(x_std);
      debug_jochen << "============================" << std::endl;
      gridGenerator.checkForDoublettes(x);
      debug_jochen << "===== y: ===================" << std::endl;
      Output::printVector(y_std);
      debug_jochen << "============================" << std::endl;
      gridGenerator.checkForDoublettes(y);

      debug_jochen << "x.size() = " << x.size() << std::endl;
      debug_jochen << "y.size() = " << y.size() << std::endl;
      debug_jochen << "=> Will give grid of " << ((x.size()-1) * (y.size()-1)) << " cells!" << std::endl;


      Ax1ElementGroupMapper elemGroupMapper;
      ElementSubdomainMapper elemSubdomainMapper;
      gridGenerator.makeMultidomainGrid(x,y,helper,membGroups,
          elemGroupMapper,elemSubdomainMapper);

      //debug_jochen << "Grid size: " << grids[i]->size(0) << std::endl;

      elemGroupMapper_last = elemGroupMapper;
      elemSubdomainMapper_last = elemSubdomainMapper;
     }

    grids = gridGenerator.getGrids();

    // Check grids
    for(int i=0; i<grids.size(); i++)
    {
      debug_jochen << "Grid size: " << grids[i]->size(0) << std::endl;
    }

    Grid& grid_last = *(grids.back());

    typedef Grid::LeafGridView GV;
    GV gv = grid_last.leafGridView();

    typedef Grid::SubDomainGrid SubDomainGrid;
    SubDomainGrid& elecGrid = grid_last.subDomain(0);
    SubDomainGrid& membGrid = grid_last.subDomain(1);

    typedef SubDomainGrid::LeafGridView SDGV;
    SDGV elecGV = elecGrid.leafGridView();
    SDGV membGV = membGrid.leafGridView();

    bool configFound = false;
    std::string configName = params.getConfigName();

    if(configName == "default")
    {
      DefaultConfiguration<double> config(params);
      run(grids,params,config,gv,dtstart,tend,
          elemSubdomainMapper_last,elemGroupMapper_last,elecGV,membGV);
      configFound = true;
    }

    if(! configFound)
    {
      DUNE_THROW(Dune::Exception, "No configuration named '" << configName << "' could be found!");
    }

    double elapsed = timer.stop();
    debug_info << "Time elapsed: " << elapsed << " s" << std::endl;

  } catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    double elapsed = timer.stop();
    debug_info << "Time elapsed: " << elapsed << " s" << std::endl;

    throw e; // Throw exception in order to get a full stacktrace for debugging

  }

}
