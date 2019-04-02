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
#define USE_CACHE 1

//#define DUNE_ISTL_WITH_CHECKING

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <typeinfo>
#include <execinfo.h>
#include <signal.h>
//#include <ucontext.h>
#include <unistd.h>

// Steffen's global timer for new C++11 low-overhead timings
#include <dune/ax1/common/timers.hh>

#include <dune/common/array.hh>
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
//#include <dune/ax1/acme2_cyl/common/acme2_cyl_mori_setup.hh>
//#include <dune/ax1/acme2_cyl/common/laplace_setup.hh>
#include <dune/ax1/channels/channel_builder.hh>
#include <dune/ax1/common/ax1_lfs_tools.hh>
#include <dune/ax1/common/ax1_parallelhelper.hh>

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
  for(int i=0; i<Config::numberOfSpecies; ++i)
  {
    debug_verb << "con_diffWater" << ION_NAMES[i] << " = " << physics.getElectrolyte().getDiffConst(i) << std::endl;
  }
  debug_verb << "electrolyte temperature: " << physics.getElectrolyte().getTemperature() << std::endl;
  debug_verb << "electrolyte stdCon: " << physics.getElectrolyte().getStdCon() << std::endl;
  debug_verb << "poissonConstant: " << physics.getPoissonConstant() << std::endl;
  debug_verb << "==============================================" << std::endl;

//#ifdef USE_LAPLACE_OPERATOR
//  // Set up and start simulation
//  LaplaceSetup<Grid,GV,PHYSICS,SubGV> laplaceSetup(grid, gv, physics, elecGV, membGV);
//  laplaceSetup.setup(dtstart, tend);
//#else
//#ifdef USE_MORI_OPERATOR_SPLIT
//  if(! params.general.get("useMoriOperatorSplit",false))
//  {
//    DUNE_THROW(Dune::Exception, "You are using the operator-split executable, please set flag "
//        << "'useMoriOperatorSplit = yes' in config file!");
//  }
//
//  // Set up and start simulation
//  Acme2CylMoriSetup<Grid,GV,PHYSICS,SubGV> acme2_cylMoriSetup(grid, gv, physics, elecGV, membGV);
//  acme2_cylMoriSetup.setup(dtstart, tend);
//#else
    // Set up and start simulation
    Acme2CylSetup<Grid,GV,PHYSICS,SubGV> acme2_cylSetup(grid, gv, physics, elecGV, membGV);
    acme2_cylSetup.setup(dtstart, tend);
//#endif
//#endif
}

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "[dune-ax1 backtrace handler] Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

///* This structure mirrors the one found in /usr/include/asm/ucontext.h */
//typedef struct _sig_ucontext {
// unsigned long     uc_flags;
// struct ucontext   *uc_link;
// stack_t           uc_stack;
// struct sigcontext uc_mcontext;
// sigset_t          uc_sigmask;
//} sig_ucontext_t;
//
//void crit_err_hdlr(int sig_num, siginfo_t * info, void * ucontext)
//{
//  //std::ostream& err_str = std::cerr;
//  std::ofstream err_str("ax1_stacktrace.log", std::ios_base::trunc);
//
//
//  err_str << "[Ax1 crit_err_hdlr] Generating backtrace..." << std::endl;
//    sig_ucontext_t * uc = (sig_ucontext_t *)ucontext;
//
////    void * caller_address = (void *) uc->uc_mcontext.eip; // x86 specific
////
////    std::cerr << "signal " << sig_num
////              << " (" << strsignal(sig_num) << "), address is "
////              << info->si_addr << " from " << caller_address
////              << std::endl << std::endl;
//
//    void * array[50];
//    int size = backtrace(array, 50);
//
////    array[1] = caller_address;
//
//    char ** messages = backtrace_symbols(array, size);
//
//    // skip first stack frame (points here)
//    for (int i = 1; i < size && messages != NULL; ++i)
//    {
//        char *mangled_name = 0, *offset_begin = 0, *offset_end = 0;
//
//        // find parantheses and +address offset surrounding mangled name
//        for (char *p = messages[i]; *p; ++p)
//        {
//            if (*p == '(')
//            {
//                mangled_name = p;
//            }
//            else if (*p == '+')
//            {
//                offset_begin = p;
//            }
//            else if (*p == ')')
//            {
//                offset_end = p;
//                break;
//            }
//        }
//
//        // if the line could be processed, attempt to demangle the symbol
//        if (mangled_name && offset_begin && offset_end &&
//            mangled_name < offset_begin)
//        {
//            *mangled_name++ = '\0';
//            *offset_begin++ = '\0';
//            *offset_end++ = '\0';
//
//            int status;
//            char * real_name = abi::__cxa_demangle(mangled_name, 0, 0, &status);
//
//            // if demangling is successful, output the demangled function name
//            if (status == 0)
//            {
//                err_str << "[bt]: (" << i << ") " << messages[i] << " : "
//                          << real_name << "+" << offset_begin << offset_end
//                          << std::endl;
//
//            }
//            // otherwise, output the mangled function name
//            else
//            {
//                err_str << "[bt]: (" << i << ") " << messages[i] << " : "
//                          << mangled_name << "+" << offset_begin << offset_end
//                          << std::endl;
//            }
//            free(real_name);
//        }
//        // otherwise, print the whole line
//        else
//        {
//          err_str << "[bt]: (" << i << ") " << messages[i] << std::endl;
//        }
//    }
//    err_str << std::endl;
//
//    err_str.close();
//
//    free(messages);
//
//    exit(EXIT_FAILURE);
//}

//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  signal(SIGSEGV, handler);
//  signal(SIGABRT, handler);

//  struct sigaction sigact;
//
//   sigact.sa_sigaction = crit_err_hdlr;
//   sigact.sa_flags = SA_RESTART | SA_SIGINFO;
//
//   if (sigaction(SIGSEGV, &sigact, (struct sigaction *)NULL) != 0)
//   {
//    fprintf(stderr, "error setting signal handler for %d (%s)\n",
//      SIGSEGV, strsignal(SIGSEGV));
//
//    exit(EXIT_FAILURE);
//   }
//   if (sigaction(SIGABRT, &sigact, (struct sigaction *)NULL) != 0)
//   {
//    fprintf(stderr, "error setting signal handler for %d (%s)\n",
//      SIGABRT, strsignal(SIGSEGV));
//
//    exit(EXIT_FAILURE);
//   }

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
    GridVector<double> x;
    CylinderGridVector<double> y(params.dX());
    std::vector<typename Ax1GridGenerator::MembraneGroupTuple> membGroups;
    gridGenerator.generateTensorGrid(x, y, membGroups);

    // Hack for Mori model to remove grid vertices close to the membrane
    if(params.general.get("removeDebyeLayerVertices",false))
    {
      debug_jochen << "===== y, original: ===================" << std::endl;
      Output::printVector(y);
      debug_jochen << "============================" << std::endl;

      typename GridVector<double>::iterator start = y.end();
      typename GridVector<double>::iterator end = y.end();
      bool foundStart = false;
      double hmin = -1;

      double yMemb = params.yMemb()[0];
      double dMemb = params.dMemb();

      // If not using membrane, yMemb[0] is set to 0.0; workaround this by forcing usage of
      // the (otherwise unused) value if yMemb in config file; this is hack anyway, so it doesn't matter.
      if(! params.useMembrane())
      {
        yMemb = params.general.get("y_memb",std::vector<double>(1,0.0))[0];
        dMemb = params.general.get("d_memb",dMemb);
        debug_jochen << "Using yMemb = " << yMemb << ", dMemb = " << dMemb << std::endl;
      }

      double dmin = params.general.get("minimalMembraneDistance",params.dYMin());
      for(typename GridVector<double>::iterator it = y.begin(); it != y.end(); ++it)
      {
        if(!foundStart && *it > yMemb + dMemb)
        {
          start = it;
          foundStart = true;
        }
        if(foundStart && (*it-(yMemb + dMemb) > dmin
            || std::abs(*it-(yMemb + dMemb) - dmin) < 1e-6))
        {
          // This is the spacing between the last non-removed vertex and its following neighbor
          hmin = -*it + *(it+1);
          end = it;
          break;
        }
      }
      if(start != end)
      {
        // Delete all vertices with membrane distance smaller than desired value
        debug_info << "Removing interval [" << *start  << ", " << *end << ")" << std::endl;
        end = y.erase(start,end);

        debug_jochen << "===== y, removed: ===================" << std::endl;
        Output::printVector(y);
        debug_jochen << "============================" << std::endl;

        // Now fill the resulting gap with equidistant-spaced vertices to avoid one huge element
        if(params.general.get("fillRemovedRange",false))
        {
          start = end-1;
          int nEquidistant = (int) std::floor((*end - *start) / hmin);
          double hEquidistant = (*end - *start) / nEquidistant;

          debug_info << "Filling interval (" << *start << ", " << *end << ") with equidistant elements of size "
              << hEquidistant << std::endl;

          debug_info << "hmin = " << hmin << std::endl;
          debug_info << "nEquidistant = " << nEquidistant << std::endl;
          debug_info << "hEquidistant = " << hEquidistant << std::endl;

          for(int i=nEquidistant-1; i>0; --i)
          {
            debug_verb << "Inserting vertex " << (*start + i*hEquidistant) << std::endl;
            end = y.insert(end, *start + i*hEquidistant);
          }
        }
      }
    }

    params.setMembraneGroups(membGroups);

    // New method globalRefine which basically does bisection, but excludes membrane elements in y-direction
    if(params.general.get("doAdditionalGlobalRefine",false))
    {
      gridGenerator.globalRefine(x,y,level);
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

    if(params.doEquilibration() && x.size() > 2 && !params.doForceEquilibration())
      DUNE_THROW(Dune::Exception, std::string("The new data transfer only works with grids having exactly 1 element ")
         + std::string("in x-direction! Please consider adjusting your grid parameters."));


    typedef Ax1GridGenerator::Grid Grid;

    Ax1ElementGroupMapper elemGroupMapper;
    ElementSubdomainMapper elemSubdomainMapper;

    gridGenerator.makeMultidomainGrid(x,y,helper,membGroups,elemGroupMapper,elemSubdomainMapper);

    Grid& grid = *gridGenerator.getGrids()[0];

    typedef Grid::LeafGridView GV;
    GV gv = grid.leafGridView();

    typedef Grid::SubDomainGrid SubDomainGrid;
    SubDomainGrid& elecGrid = grid.subDomain(0);
    SubDomainGrid& membGrid = grid.subDomain(1);

    typedef SubDomainGrid::LeafGridView SDGV;
    SDGV elecGV = elecGrid.leafGridView();
    SDGV membGV = membGrid.leafGridView();


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

//#ifdef USE_LAPLACE_OPERATOR
//    if(configName == "laplace")
//        {
//          LaplaceConfiguration<double> config;
//          run(grid,params,config,gv,dtstart,tend,
//              elemSubdomainMapper,elemGroupMapper,elecGV,membGV);
//          configFound = true;
//        }
//#else
//    if(configName == "mori")
//    {
//      MoriConfiguration<double> config(params);
//      run(grid,params,config,gv,dtstart,tend,
//          elemSubdomainMapper,elemGroupMapper,elecGV,membGV);
//      configFound = true;
//    }
//  #ifndef USE_MORI_OPERATOR_SPLIT
    if(configName == "default")
    {
      DefaultConfiguration<double> config(params);
      run(grid,params,config,gv,dtstart,tend,
          elemSubdomainMapper,elemGroupMapper,elecGV,membGV);
      configFound = true;
    }
//    if(configName == "ES")
//    {
//      ESConfiguration<double> config;
//      run(grid,params,config,gv,dtstart,tend,
//          elemSubdomainMapper,elemGroupMapper,elecGV,membGV);
//      configFound = true;
//    }
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
//  #endif
//#endif

    if(! configFound)
    {
      DUNE_THROW(Dune::Exception, "No configuration named '" << configName << "' could be found!");
    }

    if(helper.rank() == 0)
      timers.print_timers_per_call(std::cout);

    double elapsed = timer.stop();
    debug_info << "Time elapsed: " << elapsed << " s" << std::endl;

  } catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;

    timers.print_timers_per_call(std::cout);

    double elapsed = timer.stop();
    debug_info << "Time elapsed: " << elapsed << " s" << std::endl;

    throw e; // Throw exception in order to get a full stacktrace for debugging
  }
}
