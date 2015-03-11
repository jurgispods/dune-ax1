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

#if HAVE_DUNE_SUBGRID && USE_SUBGRID==2
#include <dune/subgrid/subgrid.hh>
#endif

#if HAVE_DUNE_MULTIDOMAINGRID
#include <dune/grid/multidomaingrid/multidomaingrid.hh>
#endif

#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinggridfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinglocalfunctionspace.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/multidomain/constraints.hh>

#include <dune/ax1/common/ax1_subgridview.hh>
#include <dune/ax1/common/ax1_subgrid_hostelement_mapper.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/acme1MD/common/acme1MD_factory.hh>
#include <dune/ax1/acme1MD/common/acme1MD_physics.hh>
#include <dune/ax1/acme1MD/common/acme1MD_parametertree.hh>
#include <dune/ax1/acme1MD/common/acme1MD_setup.hh>
#include <dune/ax1/channels/channel_builder.hh>


template<typename Config, typename GV, typename SubGV>
void run(Acme1MDParameters& params, Config& config, GV& gv, std::vector<double>& coords, double dtstart, double tend,
    ElementSubdomainMapper& elemSubdomainMapper, SubGV& elecGV, SubGV& membGV)
{
  debug_verb << "Configuration '" << config.getConfigName() << "' loaded." << std::endl;

  typedef typename Acme1MDFactory<GV,double,Config,SubGV>::PHYSICS PHYSICS;
  PHYSICS physics = Acme1MDFactory<GV,double,Config,SubGV>::setup(config, gv, elecGV, membGV, params, coords, dtstart,
      elemSubdomainMapper);


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

  double v_rest = -65.0;
  double dt = dtstart;

  for(double time = 0.0; time < tend; time +=dt)
  {
    std::cout << "time = " << time << std::endl;
    typename PHYSICS::ChannelSet& channels = physics.getMembrane().getChannelSet();
    std::map<int,int> membraneElements = channels.getMembraneElements();

    for(typename PHYSICS::SubDomainElementIterator sdeit = membGV.template begin<0>();
        sdeit != membGV.template end<0>(); ++sdeit)
    {
      v_rest += 1;


      int eIndex = physics.getElementIndex(*sdeit);
      if(time == 0.0)
      {
        // Initialize channels to steady state with respect to resting potential
        channels.initChannels(eIndex, v_rest);
      } else {
        // Carry channels one step further; when v_rest does not change, the channels' conductances should not change either!
        double real_dt = 1e-9 * dt;
        channels.timeStep(eIndex, real_dt, v_rest);

      }


      std::cout << "---------------------------------------" << std::endl;
      std::cout << "Membrane element #" <<  eIndex
          << " @" << sdeit->geometry().center() << std::endl;

      int mIndex = membraneElements[eIndex];
      for(int k=0; k<channels.size(); k++)
      {
        int nGatingParticles = channels.getChannel(k).numGatingParticles();
        for(int j=0; j<nGatingParticles; j++)
        {
          std::cout << "p_" << j << " = " << channels.getGatingParticle(k,j,eIndex)
              << std::endl;
        }
        std::cout << "g_" << mIndex << " = " << channels.getEffConductance(k,eIndex)
            << " --> g_" << mIndex << " = " << channels.getNewEffConductance(k,eIndex)
            << std::endl;
      }
      std::cout << std::endl;
    }
    channels.updateState();

  }
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

    // sequential version
    if (1 && helper.size()==1)
    {
      // Read config file
      Acme1MDParameters params;
      Dune::ParameterTreeParser::readINITree(configFileName, params, true);
      params.init();
      
      assert(params.useMembrane() == true);

      // ############## grid generation ##################
      // Seed 1d grid coordinates
      double d = params.dMemb();
      std::vector<double> coords = { 0.5*d, params.xMax() };

      Tools::globalRefineVector( coords, level, params.refineMembrane() );
      //Tools::globalRefineVectorLogarithmic( coords, level );
      //Output::printVector(coords);

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
      //Output::printVector(coords);
      
      typedef Dune::OneDGrid BaseGrid;
      typedef BaseGrid GridType;
      GridType baseGrid(coords);
      // ############## END grid generation ###############
      
      
      // ============== Multidomaingrid stuff ================================================================
      typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,2> > Grid;
      Grid grid(baseGrid,false);

      typedef Grid::LeafGridView GV;
      GV gv = grid.leafView();

      typedef GV::Codim<0>::Iterator ElementIterator;
      typedef GV::Codim<0>::EntityPointer ElementPointer;
      typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

      // Custom mapper to map element indices to group indices
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

        // Store group of this element [M->G]
        elemSubdomainMapper.setGroup(elemIndex, subdomainIndex);
      }

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
              DUNE_THROW(Dune::Exception, "Element could not be assigned to a pre-defined group!");
          }

        }
      grid.preUpdateSubDomains();
      grid.updateSubDomains();
      grid.postUpdateSubDomains();
      // ===============================================================================================

      bool configFound = false;
      std::string configName = params.getConfigName();
      if(configName == "default")
      {
        DefaultConfiguration<double> config;
        run(params,config,gv,coords,dtstart,tend,
            elemSubdomainMapper,elecGV,membGV);
        configFound = true;
      }
      if(configName == "hamburger")
      {
        HamburgerConfiguration<double> config;
        run(params,config,gv,coords,dtstart,tend,
            elemSubdomainMapper,elecGV,membGV);
        configFound = true;
      }

      if(! configFound)
      {
        DUNE_THROW(Dune::Exception, "No configuration named '" << configName << "' could be found!");
      }

      double elapsed = timer.stop();
      debug_info << "Time elapsed: " << elapsed << " s" << std::endl;
    }
  } catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 0;
  }

}
