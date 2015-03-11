// class ChannelSet
// Simulate sets of nonlinear ion channels
// author: Stefan Lang

#ifndef CHANNELSET_HH
#define CHANNELSET_HH

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

#include <dune/istl/bvector.hh>

#include <dune/ax1/channels/channel.hh>


/**
 * \brief This class provides a consistent channel set containing pointers to
 * different channel classes (with a common interface).
 * It contains a map of a subset of element indices associated with membrane
 * element indices.
 * All the methods operating on grid elements need to map the element index
 * to the membrane index before calling the channel methods, as the channels
 * exclusively use membrane indices.
 */
template <class T, class V>
class DynamicChannelSet
{

public:
  typedef Dune::FieldVector<T, NUMBER_OF_SPECIES> VCON;

  typedef Channel<T,V> ChannelType;
  typedef LeakChannel<T,V> LeakChannelType;
  typedef VoltageGatedChannel<T,V> VoltageGatedChannelType;

  static constexpr double TIME_SCALE = 1e-3;

  DynamicChannelSet()
  {}

  inline void addChannel(Dune::shared_ptr<Channel<T,V> > channel)
  {
    channels.push_back(channel);
  }

  //! \brief Add an element index to the map of membrane indices
  inline void addMembraneElement(int iIndex)
  {
    // Successively add membrane elements, assume size()-1 is last index present in map!
    int membIndex = membraneIndices.size();

    //debug_verb << "Adding intersection #" << iIndex  << " with new membrane index "
    //    << membIndex << " to channelset!" << std::endl;

    membraneIndices[iIndex] = membIndex;
  }


  /**
   * \brief Check for a channel with the given name
   * \param[in] name
   * \return True if a channel with the given name is present in
   * this channel set
   */
  inline bool hasChannel(std::string name) const
  {
    for(int k=0; k<channels.size(); ++k)
    {
      if(channels[k]->getName() == name) {
        return true;
      }
    }
    return false;
  }

  /**
   * \brief Get Channel with the given name. Use hasChannel(string) to check if it exists first!
   * \param[in] name
   * \return The channel for a given name if it exists, a useless preinitialized shared_ptr else
   */
  inline Dune::shared_ptr<Channel<T,V> > getChannelPtr(std::string name) const
  {
    for(int k=0; k<channels.size(); ++k)
    {
      if(channels[k]->getName() == name) return channels[k];
    }

    DUNE_THROW(Dune::Exception, "No channel named '" << name << "' present in this channel set!");
  }

  inline Dune::shared_ptr<Channel<T,V> > getChannelPtr(int k)
  {
    assert(k < channels.size());
    return channels[k];
  }

  inline const Channel<T,V>&  getChannel(int k) const
  {
    assert(k < channels.size());
    return *(channels[k]);
  }


  inline int size() const
  {
    return channels.size();
  }


  virtual inline void resize ()
  {
    int nMembElements = membraneIndices.size();
    for(int k=0; k<channels.size(); ++k)
    {
      channels[k]->resize(nMembElements);
    }
  }

  virtual inline void init()
  {
    DUNE_THROW(Dune::NotImplemented, "Use init(Physics&) instead!");
  }

  //! \brief This is the first init() method that is called by Physics. It resizes the conductances vectors
  // to match the number of elements and initializes it with values from the config file
  template<class PHYSICS>
  inline void init (const PHYSICS& physics)
  {
    if(physics.getParams().useMembrane())
    {
      T sumLeakRatios = 0.0;
      for(int k=0; k<channels.size(); ++k)
      {
        // init channel states
        // Note: This must be called before setting conductances, otherwise they will be set to zero!
        channels[k]->init();

        debug_jochen << "Initializing channel " << channels[k]->getName() << std::endl;

        const typename PHYSICS::GridView& gv = physics.gridView();
        for(typename PHYSICS::MIterator mit = physics.mBegin(); mit != physics.mEnd(); ++mit)
        {
          //debug_jochen << "Membrane element @" << mit->geometry().center() << std::endl;

          int i = physics.getIntersectionIndex(*mit);
          int m = membraneIndices[i];
          std::string groupName = physics.getGroupName(*mit->inside());
          std::string chName = channels[k]->getName();

          //debug_jochen << "Membrane interface #" << m << " (" << i << "), group: " << groupName << std::endl;

          T conductance = 0.0;
          if(channels[k]->isLeakChannel())
          {
            // This channel is a leak channel, get total leak conductance and ratio contributed by this channel
            conductance = physics.getParams().membrane.sub(groupName).get("leak", 0.);
            T ratio =  physics.getParams().equilibration.get(chName, 0.);

            sumLeakRatios += ratio;

            setRatio(k, i, ratio);
          } else {
            conductance = physics.getParams().membrane.sub(groupName).get(chName, 0.);
          }
          channels[k]->setConductance(m, conductance);
        }

      }
      if(std::abs(sumLeakRatios - membraneIndices.size()) > 1e-6)
        DUNE_THROW(Dune::Exception,
          "The sum of leak conductance ratios (" << sumLeakRatios
          << ") in section [equilibration] of the config file is not equal to "
          << "the number of membrane elements (" << membraneIndices.size() << ")!");
    }
  }


  //! \brief This is the second init() method that is called by MembraneFluxGridFunction. It initialized the
  // channel states of a given element i with respect to a given membrane potential
  inline void initChannels(const int i, const T& v, const VCON& conCytosol, const VCON& conExtra)
  {
    checkIndex(i);
    int m = membraneIndices[i];


    std::vector<int> secondaryChannels;
    for(int k=0; k<channels.size(); ++k)
    {
      channels[k]->setV_rest(v);

      // Secondary channels depend on the state of other channels; handle those last
      if(channels[k]->isSecondary())
      {
        secondaryChannels.push_back(k);
      } else {
        channels[k]->init(m,v,conCytosol,conExtra);
      }
    }

    // Now adjust leak conductance ratios
    recalculateLeakConductances(i);

    // Handle secondary channels last
    if(secondaryChannels.size() > 0)
    {
      // Fill map with effective conductances of primary channels
      std::map<std::string,T> gOtherChannels;
      for(int k=0; k<channels.size(); ++k)
      {
        if(! channels[k]->isSecondary())
        {
          gOtherChannels[channels[k]->getName()] = channels[k]->getEffConductance(m);
        }
      }

      // Provide secondary channels with the conductances of the other (primary channels) and do init
      for(int k=0; k<secondaryChannels.size(); ++k)
      {
        channels[secondaryChannels[k]]->init(m,v,conCytosol,conExtra,gOtherChannels);
      }
    }
  }

  /**
   * Calculate new leak channel conductances after voltage-gated channel have been activated
   * @param i element index (not membrane element index!)
   */
  inline void recalculateLeakConductances(const int i)
  {
    debug_verb << "----------------------------------------" << std::endl;
    // Loop over all possible ion species
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      // Loop over all channels existing on this membrane element
      for(int k=0; k<channels.size(); k++)
      {
        // Is this channel 'k' permeable for the current ion species 'j' at all?
        if(channels[k]->getIonSpecies() == j)
        {

          debug_verb << channels[k]->getName() << ": ";
          if(not channels[k]->isLeakChannel())
          {
            debug_verb << "gbar = " << getConductance(k, i)
                << ", g_eff = " << getNewEffConductance(k, i) << std::endl;
          } else {

            // The total leak conductance is available via the maximum conductance of this leak channel
            T totalLeakConductance = getConductance(k, i);

            // Only update leak conductances ratios if total leak conductance is nonzero
            if(totalLeakConductance > 0.0)
            {
              // Separately calculate voltage-gated conductances for this ion species and other ion species
              T sumVoltageGatedConductancesThis = 0.0;
              T sumVoltageGatedConductancesOther = 0.0;
              for(int l=0; l<channels.size(); l++)
              {
                if(channels[l]->isVoltageGated())
                {
                 if(channels[l]->getIonSpecies() != j)
                 {
                   sumVoltageGatedConductancesOther += getNewEffConductance(l, i);
                 } else {
                   sumVoltageGatedConductancesThis += getNewEffConductance(l, i);
                 }
                }
              }
              // Get this channel's ratio of the total leak conductance
              T ratio = getGatingParticle(k, 0, i);

              T leakConductance = 0.0;
              leakConductance = (ratio/(1.-ratio)) * (totalLeakConductance + sumVoltageGatedConductancesOther)
                  - sumVoltageGatedConductancesThis;
              leakConductance /= 1. + (ratio/(1.-ratio));

              // Now calculate the ratio of this channel's leak conductance from the total leak conductance
              T newRatio = leakConductance / totalLeakConductance;

              /*
              debug_jochen << "totalLeakConductance = " << totalLeakConductance << std::endl;
              debug_jochen << "sumVoltageGatedConductancesThis = " << sumVoltageGatedConductancesThis << std::endl;
              debug_jochen << "sumVoltageGatedConductancesOther = " << sumVoltageGatedConductancesOther << std::endl;
              debug_jochen << "ratio = " << ratio << std::endl;
              debug_jochen << "leakConductance = " << leakConductance << std::endl;
              debug_jochen << "newRatio = " << newRatio << std::endl;
              */

              debug_verb << "Calculated [" << ION_NAMES[j] << "] leak conductance: " << leakConductance
                           << " (ratio " << newRatio << ", old ratio: " << ratio << ")" << std::endl;

              if(newRatio < 0.0 || newRatio > 1.0)
                DUNE_THROW(Dune::Exception, "Calculated leak channel ratio is " << newRatio << ", must be between 0 and 1!");

              setRatio(k, i, newRatio);
            } else {
              debug_verb << "Total leak conductance is 0, did not recalculate leak conductance ratios" << std::endl;
            }

          }
        }
      }
    }
    debug_verb << "----------------------------------------" << std::endl;
  }

  //! \brief Do one time step for membrane element i
  //! \param dt Real dt [s]
  //! \param v Potential [mV]
  virtual inline void timeStep (int i, const T dt, const T v,
      const VCON& conCytosol, const VCON& conExtra)
  {
    checkIndex(i);

    T dt_ms = dt / TIME_SCALE;
    debug_verb << "  Channels.timeStep(" << i << ", " << dt_ms << "[ms], "
        << v << "[mV])" << std::endl;
    debug_verb << "   real time step: dt = " << dt << "[s]" << std::endl;
    int m = membraneIndices[i];
    for(int k=0; k<channels.size(); ++k)
    {
      channels[k]->timeStep(m,dt_ms,v,conCytosol,conExtra);
    }
  }

  virtual inline void acceptState()
  {
    for(int k=0; k<channels.size(); ++k)
    {
      channels[k]->updateState();
    }
  }

  virtual inline T getConductance(int k, int i) const
  {
    checkIndex(i);
    int m = membraneIndices.at(i);
    assert(k < channels.size());
    return channels[k]->getConductance(m);
  }

  virtual inline void setConductance(int k, int i, const T g_)
  {
    checkIndex(i);
    int m = membraneIndices.at(i);
    assert(k < channels.size());
    channels[k]->setConductance(m,g_);
  }

  virtual inline void setRatio(int k, int i, const T ratio_)
  {
    checkIndex(i);
    int m = membraneIndices.at(i);
    assert(k < channels.size());
    assert(channels[k]->isLeakChannel());
    // Now this is ugly. I f***ing hate C++.
    dynamic_cast<LeakChannelType *>(channels[k].get())->setRatio(ratio_, m);
  }

  virtual inline T getEffConductance(int k, int i) const
  {
    checkIndex(i);
    int m = membraneIndices.at(i);
    assert(k < channels.size());
    return channels[k]->getEffConductance(m);
  }

  virtual inline T getNewEffConductance(int k, int i) const
  {
    checkIndex(i);
    int m = membraneIndices.at(i);
    assert(k < channels.size());
    return channels[k]->getNewEffConductance(m);
  }


  virtual inline T getGatingParticle (int k, int j, int i) const
  {
    checkIndex(i);
    int m = membraneIndices.at(i);
    assert(k < channels.size());
    return channels[k]->getGatingParticle(j,m);
  }

  virtual inline T getNewGatingParticle (int k, int j, int i) const
  {
    checkIndex(i);
    int m = membraneIndices.at(i);
    assert(k < channels.size());
    return channels[k]->getNewGatingParticle(j,m);
  }


  virtual inline void info() const
  {
    debug_info << "-------- ChannelSet info -----------" << std::endl;
    for(int k=0; k<channels.size(); ++k)
    {
      channels[k]->info();
    }
    debug_info << "------------------------------------" << std::endl;
  }

  inline virtual void updateConcentration (const V& conc)
  {
    for(int k=0; k<channels.size(); ++k)
    {
      channels[k]->updateConcentration(conc);
    }
  }

  int getMembraneIndex(int iIndex) const
  {
    return membraneIndices.at(iIndex);
  }

  const std::map<int, int>& getMembraneIndices() const
  {
    return membraneIndices;
  }

  const std::vector<T> serializeChannelStates() const
  {
    int vSize = 0;

    // Determine total number of state values in this channelset
    for(int k=0; k<channels.size(); ++k)
    {
      vSize += channels[k]->numGatingParticles();
    }

    // Multiply this by the number of membrane elements
    vSize *= membraneIndices.size();

    // Initialize state vector
    std::vector<T> states(vSize,0.0);

    //debug_jochen << "Initialized state vector to size " << vSize << std::endl;

    typename std::vector<T>::iterator it = states.begin();


    // Loop over channels
    for(int k=0; k<channels.size(); ++k)
    {
      // Loop over this channel's gating particles
      for(int j=0; j<channels[k]->numGatingParticles(); j++)
      {
        // Get vector of all elements' gating particles at once
        const V& g = channels[k]->getGatingParticle(j);

        it = std::copy(g.begin(),g.end(),it);

        //debug_jochen << "Current state vector:" << std::endl;
        //Output::printVector(states);
      }
    }

    //debug_jochen << "overall number of gating particles: " << states.size() << std::endl;
    return states;
  }

  void deserializeChannelStates(const std::vector<T>& states)
  {
    int vSize = 0;

    // Determine total number of state values
    for(int k=0; k<channels.size(); ++k)
    {
      vSize += channels[k]->numGatingParticles();
    }
    vSize *= membraneIndices.size();

    if(vSize != states.size())
      DUNE_THROW(Dune::Exception, "Error loading channel states, vector size does not match number of states!");

    typename std::vector<T>::const_iterator it = states.begin();

    // Loop over channels
    for(int k=0; k<channels.size(); ++k)
    {
      // Loop over this channel's gating particles
      for(int j=0; j<channels[k]->numGatingParticles(); j++)
      {
        // Get vector of all elements' gating particles at once
        V g(membraneIndices.size());

        std::copy(it,it+membraneIndices.size(),g.begin());
        it += membraneIndices.size();

        //debug_jochen << "Restored gating particle vector: " << g << std::endl;

        channels[k]->setGatingParticle(j,g);
      }
    }

    assert(it == states.end());
    debug_info << "Successfully loaded vector of " << states.size() << " channel states!" << std::endl;
    Output::printVector(states);
  }



private:
  void checkIndex(const int i) const
  {
    if(membraneIndices.find(i) == membraneIndices.end())
    {
      DUNE_THROW(Dune::Exception, "Given element #" << i << " is not a membrane element!");
    }
  }


  std::vector<Dune::shared_ptr<Channel<T,V> > > channels;
  std::map<int,int> membraneIndices;

};

#endif
