/*
 * acme2_cyl_factory.hh
 *
 *  Created on: Oct 13, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_FACTORY_HH
#define DUNE_AX1_ACME2CYL_FACTORY_HH

#include <dune/ax1/common/ax1_subgridview.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_physics.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>
#include <dune/ax1/acme2_cyl/configurations/all_configurations.hh>
#include <dune/ax1/channels/channel_builder.hh>

template<typename GV, typename T, typename Config, typename SubGV=GV>
class Acme2CylFactory
{
  public:

    typedef typename GV::Grid Grid;

    typedef typename Config::template Traits<GV,T,SubGV> ConfigTraits;

    typedef Physics<GV,T,ConfigTraits> PHYSICS;
    //static const std::vector<std::string> configIndex = {"default", "hamburger"};

    static PHYSICS setup(Config& config, const GV& gv, const SubGV& elecGV, const SubGV& membGV,
        Acme2CylParameters& params, double dt,
        ElementSubdomainMapper& elemSubdomainMapper,
        Ax1ElementGroupMapper& elemGroupMapper)
    {
      std::string configName = params.getConfigName();

      //typedef Dune::BlockVector<Dune::FieldVector<T,1> > Vector;
      //typedef DynamicChannelSet<T,Vector> ChannelSet;
      typedef typename Membrane<T>::ChannelSet ChannelSet;
      ChannelSet channelSet = ChannelBuilder<T,ChannelSet>::buildChannelSet(params);

      Membrane<T> memb(2.0, channelSet); // 2.0
      T d = params.dMemb();

      T timeScale;
      T lengthScale;
      bool configFound = false;

      // Extract all necessary stuff from config object
      Electrolyte<T> electro(config.getElectrolyte());

      // Set user-defined permittivity, if there is one in the config file
      // TODO Enable having two different electrolytes to allow different permittivities! For now, consider only
      // the value defined in the extracellular solution.
      if(params.solution_ex.hasKey("permittivity"))
      {
        electro.setPermittivity(params.solution_ex.get("permittivity", electro.getPermittivity()));
        debug_info << "= Electrolyte permittivity: " << electro.getPermittivity() << std::endl;
      }

      timeScale = config.TIME_SCALE;
      lengthScale = config.LENGTH_SCALE;

      params.setUseLogScaling(config.useLogScaling);
      params.setHasAnalyticalSolution(config.hasAnalyticalSolution);

      // Set error tolerances for main iteration
      params.setReduction(config.reduction);
      params.setAbsLimit(config.absLimit);

      PHYSICS physics(gv,elecGV,membGV,elemSubdomainMapper,elemGroupMapper,electro,memb,params,timeScale,lengthScale
          /*,subGV_Inside,subGV_Outside*/);

      physics.setTimeStep(dt);
      //######################################################################

      // Print some debug info
      physics.info();
      //physics.gridInfo();
      //physics.channelInfo();

      // Some optional tests for console spamming
      //physics.testMembraneIntersectionMapping();

      return physics;
    }
};

#endif /* DUNE_AX1_ACME2CYL_FACTORY_HH */
