#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#define USE_SUBGRID 0

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <typeinfo>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/onedgrid.hh>
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG 
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#endif

#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/acme1/common/acme1_factory.hh>
#include <dune/ax1/acme1/common/acme1_parametertree.hh>
#include <dune/ax1/acme1/common/acme1_setup_nomembrane.hh>
#include <dune/ax1/channels/channel_builder.hh>

template<typename Config, typename GV>
void run(Acme1Parameters& params, Config& config, GV& gv, std::vector<double>& coords,
    double dtstart, double tend, ElementSubdomainMapper& elemSubdomainMapper)
{
  debug_verb << "Configuration '" << config.getConfigName() << "' loaded." << std::endl;

  typedef typename Acme1Factory<GV,double,Config>::PHYSICS PHYSICS;
  PHYSICS physics = Acme1Factory<GV,double,Config>::setup(config, gv, params, coords, dtstart,elemSubdomainMapper);
  //std::cout << Tools::getTypeName(physics) << std::endl;


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

  //############# operator-split ######################
  Acme1Setup<GV,PHYSICS> acme1SetupOperatorSplit(gv, physics);
  acme1SetupOperatorSplit.setup(dtstart, tend);
  //###################################################
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
      std::cout<< "This is a sequential program." << std::endl;
    else
	  {
      if(helper.rank()==0)
        std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
	  }
    
    if (argc<4)
	  {
      if(helper.rank()==0)
        std::cout << "usage: ./nemo1 <level> <dtstart> <tend> [config-file]" << std::endl;
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
      configFileName = "nemo1.config";
    }
    debug_verb << "Using config file " << configFileName << std::endl;


    // sequential version
    if (1 && helper.size()==1)
    {


      // Read config file
      Acme1Parameters params;
      Dune::ParameterTreeParser::readINITree(configFileName, params, true);
      params.init();
      
      //assert(params.useMembrane() == false);

      // ############## grid generation ##################
      // Seed 1d grid coordinates
      double d = params.dMemb();
      std::vector<double> coords = { 0.5*d, params.xMax() };

      Tools::globalRefineVector( coords, level, params.refineMembrane() );
      //Tools::globalRefineVectorLogarithmic( coords, level );
      Output::printVector(coords);

      // mirror grid
      for (int i=coords.size()-1; i>=0; i--)
      {
        // Do not mirror coordinate x=0
        if(not std::abs(coords[i] < 1e-12))
        {
          coords.push_back( -coords[i] );
        }
      }
      sort(coords.begin(), coords.end());
      Output::printVector(coords);
      
      typedef Dune::OneDGrid GridType;
      GridType grid(coords);
      
      //typedef Dune::YaspGrid<dim>::LeafGridView GV;
      typedef Dune::OneDGrid::LeafGridView GV;
      // ############## END grid generation ###############
      
      // ############# Tag elements by type and initialize mapper #############
      const GV& gv=grid.leafView();

      typedef GV::Codim<0>::Iterator ElementIterator;
      typedef GV::Codim<0>::EntityPointer ElementPointer;
      typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

      // Custom mapper to map element indices to subdomain indices
      ElementSubdomainMapper elemSubdomainMapper;
      ElementMapper elementMapper(gv);

      // Tag elements by type
      for(ElementIterator ep = gv.begin<0>(); ep != gv.end<0>(); ++ep)
      {
        int elemIndex = elementMapper.map(*ep);
        int subdomainIndex = CYTOSOL;

        double center = ep->geometry().center();

        // We are on the membrane
        if(params.useMembrane() and std::abs(center) < 0.5*d)
        {
          subdomainIndex = MEMBRANE;
        }
        // There is a membrane and we are on the outside of the cell
        if(params.useMembrane() and center > 0.5*d)
        {
          subdomainIndex = ES;
        }

        // Store subdomain of this element [M->G]
        elemSubdomainMapper.setGroup(elemIndex, subdomainIndex);
      }
      // ############# END Tag elements by type and initialize mapper #############


      bool configFound = false;
      std::string configName = params.getConfigName();
      if(configName == "default")
      {
        DefaultConfiguration<double> config;
        run(params,config,gv,coords,dtstart,tend,elemSubdomainMapper);
        configFound = true;
      }
      if(configName == "hamburger")
      {
        HamburgerConfiguration<double> config;
        run(params,config,gv,coords,dtstart,tend,elemSubdomainMapper);
        configFound = true;
      }
      if(configName == "cheeseburger")
			{
				CheeseburgerConfiguration<double> config;
				run(params,config,gv,coords,dtstart,tend,elemSubdomainMapper);
				configFound = true;
			}
      if(configName == "bigmac")
			{
				BigmacConfiguration<double> config;
				run(params,config,gv,coords,dtstart,tend,elemSubdomainMapper);
				configFound = true;
			}
      if(configName == "step")
			{
				StepConfiguration<double> config;
				run(params,config,gv,coords,dtstart,tend,elemSubdomainMapper);
				configFound = true;
			}

      if(! configFound)
      {
        DUNE_THROW(Dune::Exception, "No configuration named '" << configName << "' could be found!");
      }
    }
  } catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 0;
  }
}
