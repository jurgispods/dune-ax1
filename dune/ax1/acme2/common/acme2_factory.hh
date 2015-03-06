/*
 * acme2_factory.hh
 *
 *  Created on: Oct 13, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2_FACTORY_HH
#define DUNE_AX1_ACME2_FACTORY_HH

#include <dune/ax1/common/ax1_subgridview.hh>
#include <dune/ax1/acme2/common/acme2_physics.hh>
#include <dune/ax1/acme2/common/acme2_parametertree.hh>
#include <dune/ax1/acme2/configurations/all_configurations.hh>
#include <dune/ax1/channels/channel_builder.hh>

template<typename GV, typename T, typename Config, typename SubGV=GV>
class Acme2Factory
{
  public:

    typedef typename GV::Grid Grid;

    typedef typename Config::template Traits<GV,T,SubGV> ConfigTraits;

    typedef Physics<GV,T,ConfigTraits> PHYSICS;
    //static const std::vector<std::string> configIndex = {"default", "hamburger"};

    static PHYSICS setup(Config& config, const GV& gv, const SubGV& elecGV, const SubGV& membGV,
        Acme2Parameters& params, double dt,
        ElementSubdomainMapper& elemSubdomainMapper)
    {
      std::string configName = params.getConfigName();

      typename Membrane<T>::ChannelSet channelSet = ChannelBuilder<T>::buildChannelSet(params);
      Membrane<T> memb(2.0, channelSet); // 2.0
      T d = params.dMemb();

      T timeScale;
      T lengthScale;
      bool configFound = false;

      Electrolyte<T> electro(config.getElectrolyte());
      timeScale = config.TIME_SCALE;
      lengthScale = config.LENGTH_SCALE;

      params.setUseLogScaling(config.useLogScaling);
      params.setHasAnalyticalSolution(config.hasAnalyticalSolution);

      // Set error tolerances for main iteration
      params.setReduction(config.reduction);
      params.setAbsLimit(config.absLimit);

      PHYSICS physics(gv,elecGV,membGV,elemSubdomainMapper,electro,memb,params,timeScale,lengthScale
          /*,subGV_Inside,subGV_Outside*/);

      physics.setTimeStep(dt);
      //######################################################################

      // Channelset containing all ion channels
      typename Membrane<T>::ChannelSet& channels = physics.getMembrane().getChannelSet();

      for(typename PHYSICS::SubDomainElementIterator eit = membGV.template begin<0>(); eit != membGV.template end<0>(); ++eit)
      {
        int elemIndex = physics.getElementIndex(*eit);
        int subdomainIndex = physics.getSubdomainIndex(*eit);

        switch(subdomainIndex)
        {
          case CYTOSOL:
            /*
            // Make sure we map the right subentity: Positions must match!
            assert(center == sep_Inside->geometry().center());
            assert(sep_Inside != subGV_Inside.template end<0>());
            subIndexMapper_Inside.addElement(subElementMapper_Inside.map(*sep_Inside),ep);
            //debug_verb << subElementMapper_Inside.map(*sep_Inside) << " -> " << elemIndex << std::endl;
            ++sep_Inside;
            */
            DUNE_THROW(Dune::Exception, "Found a non-membrane entity inside the membrane subdomain grid!");
            break;
          case ES:
            /*
            // Make sure we map the right subentity: Positions must match!
            assert(center == sep_Outside->geometry().center());
            assert(sep_Outside != subGV_Outside.template end<0>());
            subIndexMapper_Outside.addElement(subElementMapper_Outside.map(*sep_Outside),ep);
            //debug_verb << subElementMapper_Outside.map(*sep_Outside) << " -> " << elemIndex << std::endl;
            ++sep_Outside;
            */
            DUNE_THROW(Dune::Exception, "Found a non-membrane entity inside the membrane subdomain grid!");
            break;
          case MEMBRANE:
          {
            channels.addMembraneElement(elemIndex);
            break;
          }
          default:
            DUNE_THROW(Dune::Exception, "Element has an unknown group index!");
        }
      }
      channels.resize();
      channels.init(physics);

      // Print some debug info
      //physics.info();
      //physics.gridInfo();
      //channels.info();

      return physics;
    }
};

#endif /* DUNE_AX1_ACME2_FACTORY_HH */
