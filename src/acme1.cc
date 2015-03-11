#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

// 0 = no subgrid, 1 = ax1-subgrid 2 = dune-subgrid
#define USE_SUBGRID 1

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

#if HAVE_DUNE_SUBGRID && USE_SUBGRID==2
#include <dune/subgrid/subgrid.hh>
#endif

#include <dune/ax1/common/ax1_subgridview.hh>
#include <dune/ax1/common/ax1_subgrid_hostelement_mapper.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/acme1/common/acme1_factory.hh>
#include <dune/ax1/acme1/common/acme1_physics.hh>
#include <dune/ax1/acme1/common/acme1_parametertree.hh>
#include <dune/ax1/acme1/common/acme1_setup.hh>
#include <dune/ax1/channels/channel_builder.hh>

template<typename Config, typename GV, typename SubGV>
void run(Acme1Parameters& params, Config& config, GV& gv, std::vector<double>& coords, double dtstart, double tend,
    ElementSubdomainMapper& elemSubdomainMapper, SubGV& subGV_Inside, SubGV& subGV_Outside)
{
  debug_verb << "Configuration '" << config.getConfigName() << "' loaded." << std::endl;

  typedef typename Acme1Factory<GV,double,Config,SubGV>::PHYSICS PHYSICS;
  PHYSICS physics = Acme1Factory<GV,double,Config,SubGV>::setup(config, gv, params, coords, dtstart,
      elemSubdomainMapper);
  //std::cout << Tools::getTypeName(physics) << std::endl;

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

  //############# operator-split ######################
  Acme1Setup<GV,PHYSICS,SubGV> acme1SetupOperatorSplit(gv, physics, subGV_Inside, subGV_Outside);
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
        std::cout << "usage: ./acme1 <level> <dtstart> <tend> [config-file]" << std::endl;
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
      configFileName = "acme1.config";
    }
    debug_verb << "Using config file " << configFileName << std::endl;

    Dune::Timer timer;

    // sequential version
    if (1 && helper.size()==1)
    {
      // Read config file
      Acme1Parameters params;
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
      
      typedef Dune::OneDGrid GridType;
      GridType grid(coords);
      
      //typedef Dune::YaspGrid<dim>::LeafGridView GV;
      typedef Dune::OneDGrid::LeafGridView GV;
      // ############## END grid generation ###############
      
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

#if USE_SUBGRID==0
      std::cout << "Please use nemo1.cc for simulations without membrane!" << std::endl;
      exit(0);
#endif
#if USE_SUBGRID==2
      typedef Dune::SubGrid<1,GridType,true> SubGrid;
      SubGrid subGrid_Inside(grid);
      SubGrid subGrid_Outside(grid);

      // ========================== SUBGRID CREATION ========================================
      // Init subgrids
      subGrid_Inside.createBegin();
      subGrid_Outside.createBegin();
      for(ElementIterator ep = gv.begin<0>(); ep != gv.end<0>(); ++ep)
      {
        int elemIndex = elementMapper.map(*ep);
        int subdomainIndex = elemSubdomainMapper.map(elemIndex);

        // Init subgrids
        switch(subdomainIndex)
        {
          case CYTOSOL:
            subGrid_Inside.insert(*ep);
            break;
          case ES:
            subGrid_Outside.insert(*ep);
            break;
          case MEMBRANE:
            // This element does not belong to any subgrid
            break;
          default:
            DUNE_THROW(Dune::Exception, "Element could not be assigned to a pre-defined subdomain!");
        }


      }
      subGrid_Inside.createEnd();
      subGrid_Outside.createEnd();

      subGrid_Inside.report();
      subGrid_Outside.report();

      /*
      debug_verb << "=== Intracellular elements: " << subGrid_Inside.size(0) << std::endl;
      for(SubGridElementIterator ep = subGrid_Inside.template leafbegin<0>();
          ep != subGrid_Inside.template leafend<0>(); ++ep)
      {
        debug_verb << "element @ " << ep->geometry().center()
                   << "  with host entity #"
                   << elementMapper.map(*(subGrid_Inside.template getHostEntity<0>(*ep)))
                   << std::endl;
        for(SubGridIntersectionIterator iit = ep->ileafbegin(); iit != ep->ileafend(); ++iit)
        {
          debug_verb << "  intersection: neighbor()=" << iit->neighbor()
              << " , boundary()=" << iit->boundary() << std::endl;
        }

      }

      debug_verb << "=== Extracellular elements: " << subGrid_Outside.size(0) << std::endl;
      for(SubGridElementIterator ep = subGrid_Outside.template leafbegin<0>();
          ep != subGrid_Outside.template leafend<0>(); ++ep)
      {
        debug_verb << "element @ " << ep->geometry().center() << std::endl;
      }
      */


      typedef SubGrid::LeafGridView SubGV;
      const SubGV subGV_Inside = subGrid_Inside.leafView();
      const SubGV subGV_Outside = subGrid_Outside.leafView();
      // ===============================================================================================
#endif
#if USE_SUBGRID==1
      typedef Dune::OneDGrid SubGrid;

      std::vector<double> coords_inside;
      std::vector<double> coords_outside;
      for(int i=0; i<coords.size(); i++)
      {
        if(coords[i] < 0.5*d || std::abs(coords[i]+0.5*d) < 1e-12)
        {
          coords_inside.push_back(coords[i]);
          //debug_verb << "[INSIDE] Inserting node " << coords[i] << std::endl;
        }
        if(coords[i] > 0.5*d || std::abs(coords[i]-0.5*d) < 1e-12)
        {
          coords_outside.push_back(coords[i]);
          //debug_verb << "[OUTSIDE] Inserting node " << coords[i] << std::endl;
        }
      }

      SubGrid subGrid_Inside(coords_inside);
      SubGrid subGrid_Outside(coords_outside);

      typedef SubGrid::LeafGridView SubGV;
      const SubGV subGV_Inside = subGrid_Inside.leafView();
      const SubGV subGV_Outside = subGrid_Outside.leafView();

      // Dune element mappers
      typedef Dune::SingleCodimSingleGeomTypeMapper<SubGV, 0> SubElementMapper;
      SubElementMapper subElementMapper_Inside(subGV_Inside);
      SubElementMapper subElementMapper_Outside(subGV_Outside);

      typedef SubGridHostElementMapper<ElementPointer,ElementMapper> HostElementMapper;
      HostElementMapper hostElementMapper_Inside(elementMapper);
      HostElementMapper hostElementMapper_Outside(elementMapper);

      typedef SubGV::Codim<0>::Iterator SubGridElementIterator;
      SubGridElementIterator sep_Inside = subGV_Inside.begin<0>();
      SubGridElementIterator sep_Outside = subGV_Outside.begin<0>();
      for(ElementIterator eit = gv.begin<0>(); eit != gv.end<0>(); ++eit)
      {
        int elemIndex = elementMapper.map(*eit);
        int subdomainIndex = elemSubdomainMapper.map(elemIndex);
        double center = eit->geometry().center();

        switch(subdomainIndex)
        {
          case CYTOSOL:
          {
            // Make sure we map the right subentity: Positions must match!
            assert(center == sep_Inside->geometry().center());
            assert(sep_Inside != subGV_Inside.end<0>());
            ElementPointer ep(eit);
            hostElementMapper_Inside.addElement(subElementMapper_Inside.map(*sep_Inside),ep);
            //debug_verb << subElementMapper_Inside.map(*sep_Inside) << " -> " << elemIndex << std::endl;
            ++sep_Inside;
            break;
          }
          case ES:
          {
            // Make sure we map the right subentity: Positions must match!
            assert(center == sep_Outside->geometry().center());
            assert(sep_Outside != subGV_Outside.end<0>());
            ElementPointer ep(eit);
            hostElementMapper_Outside.addElement(subElementMapper_Outside.map(*sep_Outside),ep);
            //debug_verb << subElementMapper_Outside.map(*sep_Outside) << " -> " << elemIndex << std::endl;
            ++sep_Outside;
            break;
          }
          case MEMBRANE:
            break;
          default:
            DUNE_THROW(Dune::Exception, "Element has an unknown subdomain index!");
        }
      }

      Ax1SubGridView<SubGV,HostElementMapper> ax1SubGV_Inside(subGV_Inside, hostElementMapper_Inside);
      Ax1SubGridView<SubGV,HostElementMapper> ax1SubGV_Outside(subGV_Outside, hostElementMapper_Outside);
#endif

      bool configFound = false;
      std::string configName = params.getConfigName();
      if(configName == "default")
      {
        DefaultConfiguration<double> config;
        run(params,config,gv,coords,dtstart,tend,
            elemSubdomainMapper,ax1SubGV_Inside,ax1SubGV_Outside);
        configFound = true;
      }
      if(configName == "hamburger")
      {
        HamburgerConfiguration<double> config;
        run(params,config,gv,coords,dtstart,tend,
            elemSubdomainMapper,ax1SubGV_Inside,ax1SubGV_Outside);
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
