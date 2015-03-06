/*
 * acme0_factory.hh
 *
 *  Created on: Oct 13, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME0_FACTORY_HH
#define DUNE_AX1_ACME0_FACTORY_HH

#include <dune/ax1/acme0/common/acme0_physics.hh>
#include <dune/ax1/acme0/common/acme0_parametertree.hh>
#include <dune/ax1/channels/channel_builder.hh>

#include <dune/ax1/acme0/configurations/default/default_config.hh>
#include <dune/ax1/acme0/configurations/hamburger/hamburger_config.hh>

template<typename GV, typename T, typename Config>
class Acme0Factory
{
  public:
    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;
    typedef typename Config::template Traits<GV,T> ConfigTraits;
    typedef Physics<GV,T,ConfigTraits> PHYSICS;
    //static const std::vector<std::string> configIndex = {"default", "hamburger"};

    //typedef CONFIG::INITIAL_CON  INITIAL_CON;
    //typedef CONFIG::SOLUTION_CON SOLUTION_CON;
    //typedef CONFIG::SOLUTION_POT SOLUTION_POT;



    static PHYSICS setup(Config& config, const GV& gv, Acme0Parameters& params, std::vector<T>& coords, double dt)
    {
      std::string configName = params.getConfigName();

      typename Membrane<T>::ChannelSet channelSet = ChannelBuilder<T>::buildChannelSet(params);
      Membrane<T> memb(2.0, channelSet);
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

       // Dune mapper to map an entity to its index
      ElementMapper elementMapper(gv);

      PHYSICS physics(gv,elementMapper,electro,memb,params,timeScale,lengthScale);
      physics.initPosition(coords);
      physics.setTimeStep(dt);
      //######################################################################

      // Custom mapper to map element indices to group indices
      ElementSubdomainMapper elemSubdomainMapper;

      // Channelset containing all ion channels
      typename Membrane<T>::ChannelSet& channels = physics.getMembrane().getChannelSet();

      std::cout << "Initializing mappers" << std::endl;
      // Initialize Element->Group mapper and Element->MembraneElement mapper
      for(ElementIterator ep = gv.template begin<0>(); ep != gv.template end<0>(); ++ep)
      {
        int elemIndex = elementMapper.map(*ep);
        int subdomainIndex = CYTOSOL;

        double center = ep->geometry().center();

        // We are on the membrane
        if(params.useMembrane() and std::abs(center) < 0.5*d)
        {
          subdomainIndex = MEMBRANE;
          // Add this membrane element to channel set [M->ME]
          channels.addMembraneElement(elemIndex);
        }
        // There is a membrane and we are on the outside of the cell
        if(params.useMembrane() and center > 0.5*d)
        {
          subdomainIndex = ES;
        }
        // Store group of this element [M->G]
        elemSubdomainMapper.setGroup(elemIndex, subdomainIndex);

      }
      physics.setElementSubdomainMapper(elemSubdomainMapper);
      channels.resize();

      std::cout << "Channels.init()" << std::endl;
      channels.init(physics);
      physics.info();
      channels.info();

      return physics;
    }
};

#endif /* DUNE_AX1_ACME0_FACTORY_HH */
