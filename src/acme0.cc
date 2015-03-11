#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
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
#include <dune/ax1/acme0/common/acme0_factory.hh>
#include <dune/ax1/acme0/common/acme0_physics.hh>
#include <dune/ax1/acme0/common/acme0_parametertree.hh>
//#ifdef ACME0_FULLY_IMPLICIT
//#include <dune/ax1/acme0/fully-implicit/acme0_Pk_fully_implicit.hh>
//#else
#include <dune/ax1/acme0/common/acme0_setup.hh>
//#endif

#include <dune/ax1/channels/channel_builder.hh>

template<typename Config, typename GV>
void run(Acme0Parameters& params, Config& config, GV& gv, std::vector<double>& coords, double dtstart, double tend)
{
  debug_verb << "Configuration '" << config.getConfigName() << "' loaded." << std::endl;

  typedef typename Acme0Factory<GV,double,Config>::PHYSICS PHYSICS;
  PHYSICS physics = Acme0Factory<GV,double,Config>::setup(config, gv, params, coords, dtstart);
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

#ifdef ACME0_FULLY_IMPLICIT
  //############# fully-implicit ######################
  //Acme0PkFullyImplicit<GV,PHYSICS> acme0PkFullyImplicit(gv, physics);
  //acme0PkFullyImplicit.acme0_Pk(dtstart, tend);
  Acme0Setup<GV,PHYSICS> acme0SetupFullyImplicit(gv, physics, false);
  acme0SetupFullyImplicit.setup(dtstart, tend);
  //###################################################
#else
  //############# operator-split ######################
  Acme0Setup<GV,PHYSICS> acme0SetupOperatorSplit(gv, physics, true);
  acme0SetupOperatorSplit.setup(dtstart, tend);
  //###################################################
#endif
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
        std::cout << "usage: ./acme0 <level> <dtstart> <tend> [config-file]" << std::endl;
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
      configFileName = "acme0.config";
    }
    debug_verb << "Using config file " << configFileName << std::endl;


    // sequential version
    if (1 && helper.size()==1)
    {


      // Read config file
      Acme0Parameters params;
      Dune::ParameterTreeParser::readINITree(configFileName, params, true);
      params.init();

      /*
      if(params.useLogScaling())
      {
        debug_info << "=== USING LOGARITHMIC SCALING! ====================" << std::endl;
      }

      Membrane<double>::ChannelSet channelSet = ChannelBuilder<double>::buildChannelSet(params);

      // electrolyte definition ###########################
      Ion<double>    sodium(  1.0, "na" );
      Ion<double> potassium(  1.0, "k"  );
      Ion<double>  chloride( -1.0, "cl" );
      
      Solvent<double> water( 80.0 );
      //Solvent<double> water( 1.0 ); // test hack
      Electrolyte<double> electro(water, 279.45, con_mol); // concentrations now given in "1 mM"
      //Electrolyte<double> electro(water, 1.0, 1.0); // test hack
      
      electro.addIon(sodium);
      if(NUMBER_OF_SPECIES > 1) electro.addIon(potassium);
      if(NUMBER_OF_SPECIES > 2) electro.addIon(chloride);
      
      Membrane<double> memb(2.0, channelSet);
      */
      
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
      
      Dune::OneDGrid grid(coords);
      
      //typedef Dune::YaspGrid<dim>::LeafGridView GV;
      typedef Dune::OneDGrid::LeafGridView GV;
      // ############## END grid generation ###############
      
      const GV& gv=grid.leafView();
      
      /*
      //typedef GV::Codim<0>::EntityPointer ElementPointer;
      typedef GV::Codim<0>::Iterator ElementIterator;
      typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

      // Dune mapper to map an entity to its index
      ElementMapper elementMapper(gv);
      


      // physics instantiation ###############################################
      typedef Physics<GV,ElementMapper,double> PHYSICS;
      PHYSICS physics(gv,elementMapper,electro,memb,params);
      physics.initPosition(coords);
      physics.setTimeStep(dtstart);
      //######################################################################
      
      // Custom mapper to map element indices to subdomain indices
      ElementSubdomainMapper elemSubdomainMapper;

      // Channelset containing all ion channels
      Membrane<double>::ChannelSet& channels = physics.getMembrane().getChannelSet();

      std::cout << "Initializing mappers" << std::endl;
      // Initialize Element->Group mapper and Element->MembraneElement mapper
      for(ElementIterator ep = gv.begin<0>(); ep != gv.end<0>(); ++ep)
      {
        int elemIndex = elementMapper.map(*ep);
        int subdomainIndex = 0;

        double center = ep->geometry().center();

        if(params.useMembrane() and std::abs(center) < 0.5*d)
        {
          subdomainIndex = 2;
          // Add this membrane element to channel set [M->ME]
          channels.addMembraneElement(elemIndex);
        }
        // Store subdomain of this element [M->G]
        elemSubdomainMapper.setGroup(elemIndex, subdomainIndex);

      }
      physics.setElementSubdomainMapper(elemSubdomainMapper);
      channels.resize();

      std::cout << "Channels.init()" << std::endl;
      channels.init(physics);
      physics.info();
      channels.info();
      */

      bool configFound = false;
      std::string configName = params.getConfigName();
      if(configName == "default")
      {
        DefaultConfiguration<double> config;
        run(params,config,gv,coords,dtstart,tend);
        configFound = true;
      }
      if(configName == "hamburger")
      {
        HamburgerConfiguration<double> config;
        run(params,config,gv,coords,dtstart,tend);
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
