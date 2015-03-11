/*
 * acme1_factory.hh
 *
 *  Created on: Oct 13, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME1_FACTORY_HH
#define DUNE_AX1_ACME1_FACTORY_HH

#include <dune/ax1/common/ax1_subgridview.hh>
#include <dune/ax1/acme1/common/acme1_physics.hh>
#include <dune/ax1/acme1/common/acme1_parametertree.hh>
#include <dune/ax1/acme1/configurations/all_configurations.hh>
#include <dune/ax1/channels/channel_builder.hh>

template<typename GV, typename T, typename Config, typename SubGV=GV>
class Acme1Factory
{
  public:

    typedef typename GV::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

    typedef typename Config::template Traits<GV,T,SubGV> ConfigTraits;

    typedef Physics<GV,T,ConfigTraits,SubGV> PHYSICS;
    //static const std::vector<std::string> configIndex = {"default", "hamburger"};

    static PHYSICS setup(Config& config, const GV& gv, Acme1Parameters& params, std::vector<T>& coords, double dt,
        ElementSubdomainMapper& elemSubdomainMapper)
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

      PHYSICS physics(gv,elemSubdomainMapper,electro,memb,params,timeScale,lengthScale
          /*,subGV_Inside,subGV_Outside*/);

      physics.initPosition(coords);
      physics.setTimeStep(dt);
      //######################################################################

      // Channelset containing all ion channels
      typename Membrane<T>::ChannelSet& channels = physics.getMembrane().getChannelSet();

      /*
      // Subgrid element index -> element index mappers
      SubGridHostElementMapper<ElementPointer,ElementMapper>& subIndexMapper_Inside
        = physics.getSubGridElementIndexMapper_Inside();
      SubGridHostElementMapper<ElementPointer,ElementMapper>& subIndexMapper_Outside
        = physics.getSubGridElementIndexMapper_Outside();

      // Dune element mappers
      SubElementMapper subElementMapper_Inside(subGV_Inside);
      SubElementMapper subElementMapper_Outside(subGV_Outside);

      SubGridElementIterator sep_Inside = subGV_Inside.template begin<0>();
      SubGridElementIterator sep_Outside = subGV_Outside.template begin<0>();
      */

      for(ElementIterator eit = gv.template begin<0>(); eit != gv.template end<0>(); ++eit)
      {
        int elemIndex = physics.getElementIndex(*eit);
        int subdomainIndex = physics.getSubdomainIndex(*eit);
        double center = eit->geometry().center();

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
            break;
          case MEMBRANE:
          {
            channels.addMembraneElement(elemIndex);
            ElementPointer ep(eit);
            physics.addMembraneInterface(ep);
            break;
          }
          default:
            DUNE_THROW(Dune::Exception, "Element has an unknown group index!");
        }
      }
      channels.resize();

      channels.init(physics);
      physics.info();
      channels.info();

      return physics;
    }
};

#endif /* DUNE_AX1_ACME1_FACTORY_HH */
