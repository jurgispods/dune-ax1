
#ifndef DUNE_AX1_MEMBRANE_HH
#define DUNE_AX1_MEMBRANE_HH

#include <dune/istl/bvector.hh>
#include <dune/ax1/channels/channel.hh>
#include <dune/ax1/channels/channelset.hh>

template<class T>
class Membrane
{
  public:
    typedef Dune::BlockVector<Dune::FieldVector<T,1> > Vector;
    typedef DynamicChannelSet<T,Vector> ChannelSet;

    Membrane (T default_permittivity_ ) :
      default_permittivity(default_permittivity_)
    {}

    Membrane (T default_permittivity_, ChannelSet channels_) :
      default_permittivity(default_permittivity_), channels(channels_)
    {
      debug_jochen << "Init membrane, membrane indices size: " << channels.getMembraneIndices().size() <<
          std::endl;
    }

    void setThickness (T thickness_)
    {
      thickness = thickness_;
    }

    T getPermittivity() const
    {
      return default_permittivity;
    }
    
    T getDiffConst ( unsigned int ionSpecies ) const
    {
      if ( ionSpecies == Na )
      {
        return 0.0;
      }
      if ( ionSpecies == K  )
      {
        return 0.0;
      }
      if ( ionSpecies == Cl )
      {
        return 0.0;
      }

      // default value (should never be used)
      return 0.;
    }

    //! Update channel states on element i
    void timeStep(int elemIndex, const T dt, const T pot)
    {
      // Solve ODEs for all channels of element i
      channels.timeStep(elemIndex, dt, pot);
    }

    //! \brief get effective conductance for element i
    const T getDiffCoeff (int ionSpecies, int elemIndex) const
    {
      T diffCoeff = 0.;

      for(int k=0; k<channels.size(); ++k)
      {
        if(channels.getChannel(k).getIonSpecies() == ionSpecies)
        {
          diffCoeff += channels.getNewEffConductance(k,elemIndex);
        }
      }

      return diffCoeff;
    }

    //! \brief get the channel set containing all ion channels
    ChannelSet& getChannelSet()
    {
      return channels;
    }

    //! \brief get the channel set containing all ion channels
    const ChannelSet& getChannelSet() const
    {
      return channels;
    }

    void updateState()
    {
      channels.updateState();
      //TODO update concentration of each ion species for ion-dependent channels, e.g.:
      //channels.updateConcentration(calciumConcentrations);
    }

    //! \brief preliminary method to set conductances for a default HH channel set
    void setChannelDensities()
    {
      std::map<int, int>& membElements = channels.getMembraneElements();

      // Loop over all membrane elements
      for(std::map<int,int>::const_iterator it = membElements.begin();
          it!= membElements.end(); ++it)
      {
        // Set Na conductance to 120 for all membrane elements
        channels.setConductance(0, it->first, 120.);
        // Set K conductance to 36 for all membrane elements
        channels.setConductance(1, it->first, 36.);
      }
    }

  private:
    T thickness;
    T default_permittivity;
    ChannelSet channels;
};

#endif // DUNE_AX1_MEMBRANE_HH
