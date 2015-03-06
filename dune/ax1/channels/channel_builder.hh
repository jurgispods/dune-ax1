/*
 * channel_builder.hh
 *
 *  Created on: May 11, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_CHANNEL_BUILDER_HH
#define DUNE_AX1_CHANNEL_BUILDER_HH

#include <string>
#include <iostream>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/channels/channelset.hh>
#include <dune/ax1/channels/channel.hh>
#include <dune/ax1/channels/traub_channels.hh>
#include <dune/ax1/channels/ion_pump.hh>

template<class T, typename ChannelSet>
class ChannelBuilder
{
  public:

    typedef Dune::BlockVector<Dune::FieldVector<T,1> > V;

    template<typename AcmeParameters>
    static ChannelSet buildChannelSet(AcmeParameters& params)
    {

      ChannelSet channelSet;

      const std::vector<std::string>& membGroups = params.membrane.getSubKeys();

      int countLeakChannels = 0;

      std::vector<std::string> ignoreKeys = {"start", "width", "stride", "permittivity", "d_memb", "dx",
          "smoothTransition", "geometricRefinement", "dx_transition", "n_transition", "loadFilename"};

      for(int i=0; i<membGroups.size(); ++i)
      {
        debug_jochen << "Group: " << membGroups[i] << std::endl;
        std::vector<std::string> chStrings = params.membrane.sub(membGroups[i]).getValueKeys();
        for(int k=0; k<chStrings.size(); k++)
        {
          std::string ch = chStrings[k];

          // Check if string ch is not one of the ignored ones and not yet present in channelset
          if(not channelSet.hasChannel(ch)
              && std::find(ignoreKeys.begin(), ignoreKeys.end(), ch) == ignoreKeys.end())
          {
            debug_jochen << "Found new channel '" << ch << "'" <<std::endl;

            // Special case for leak channel: Search for ion-specific leak channel ratios in [equilibration]
            // and add one leak channel per ion species
            if(ch == "leak")
            {
              // Do this only once (i.e., as long as count == 0!)
              if(countLeakChannels == 0)
              {

                std::vector<std::string> leakStrings = params.equilibration.getValueKeys();


                for(int l=0; l<leakStrings.size(); l++)
                {
                  int pos = leakStrings[l].find("_leak");
                  if(pos != std::string::npos)
                  {
                    debug_jochen << " - Found leak channel " << leakStrings[l] << std::endl;

                    // We found a leak channel, let's set it up!
                    int ionSpecies = -1;

                    std::string ionName = leakStrings[l].substr(0,pos);
                    for(int j=0; j<ION_NAMES.size(); j++)
                    {
                      if(ION_NAMES[j] == ionName) {
                        ionSpecies = j;
                        break;
                      }
                    }
                    if(ionSpecies == -1) DUNE_THROW(Dune::Exception,
                        "Can not create leak channel for unknown ion species '" << ionName  << "'");

                    typedef LeakChannel<T,V> Channel;
                    Channel* channel = new Channel(ionSpecies);
                    Dune::shared_ptr<Channel> channelPtr(channel);
                    channelSet.addChannel(channelPtr);

                    countLeakChannels++;
                  }
                }
                if(countLeakChannels == 0)
                {
                  DUNE_THROW(Dune::Exception,
                      "Could not find any ion-specific leak channel definitions in [equilibration] part of config file!");
                }
              }

              continue;
            }

            if(ch == "HH-Na")
            {
              DUNE_THROW(Dune::NotImplemented,
                  "Old channel types were deactivated after changing to the concentration-dependent syntax");
//              typedef NaChannel<T,V> Channel;
//              Channel* channel = new Channel();
//              Dune::shared_ptr<Channel> channelPtr(channel);
//              channelSet.addChannel(channelPtr);
//              continue;
            }
            if(ch == "HH-K")
            {
              DUNE_THROW(Dune::NotImplemented,
                  "Old channel types were deactivated after changing to the concentration-dependent syntax");
//              typedef KChannel<T,V> Channel;
//              Channel* channel = new Channel();
//              Dune::shared_ptr<Channel> channelPtr(channel);
//              channelSet.addChannel(channelPtr);
//              continue;
            }
            if(ch == "Nav")
            {
              typedef NavChannel<T,V> Channel;
              Channel* channel = new Channel();
              Dune::shared_ptr<Channel> channelPtr(channel);
              channelSet.addChannel(channelPtr);
              continue;
            }
            if(ch == "Kv")
            {
              typedef KvChannel<T,V> Channel;
              Channel* channel = new Channel();
              Dune::shared_ptr<Channel> channelPtr(channel);
              channelSet.addChannel(channelPtr);
              continue;
            }

            // ************** Traub channels **************************

            if(ch == "Na-Traub")
            {
              typedef NaTraub<T,V> Channel;
              Channel* channel = new Channel();
              Dune::shared_ptr<Channel> channelPtr(channel);
              channelSet.addChannel(channelPtr);
              continue;
            }

            if(ch == "Ca-Traub")
            {
              typedef CaTraub<T,V> Channel;
              Channel* channel = new Channel();
              Dune::shared_ptr<Channel> channelPtr(channel);
              channelSet.addChannel(channelPtr);
              continue;
            }
            if(ch == "Ca-Traub1991")
            {
              typedef CaTraub1991<T,V> Channel;
              Channel* channel = new Channel();
              Dune::shared_ptr<Channel> channelPtr(channel);
              channelSet.addChannel(channelPtr);
              continue;
            }
            if(ch == "K(DR)-Traub")
            {
              typedef K_DRTraub<T,V> Channel;
              Channel* channel = new Channel();
              Dune::shared_ptr<Channel> channelPtr(channel);
              channelSet.addChannel(channelPtr);
              continue;
            }
            if(ch == "K(DR)-Traub1991")
            {
              typedef K_DRTraub1991<T,V> Channel;
              Channel* channel = new Channel();
              Dune::shared_ptr<Channel> channelPtr(channel);
              channelSet.addChannel(channelPtr);
              continue;
            }
            if(ch == "K(AHP)-Traub")
            {
              typedef K_AHPTraub<T,V> Channel;
              Channel* channel = new Channel();
              Dune::shared_ptr<Channel> channelPtr(channel);
              channelSet.addChannel(channelPtr);
              continue;
            }
            if(ch == "K(A)-Traub")
            {
              typedef K_ATraub<T,V> Channel;
              Channel* channel = new Channel();
              Dune::shared_ptr<Channel> channelPtr(channel);
              channelSet.addChannel(channelPtr);
              continue;
            }
            if(ch == "K(C)-Traub")
            {
              typedef K_CTraub<T,V> Channel;
              Channel* channel = new Channel();
              Dune::shared_ptr<Channel> channelPtr(channel);
              channelSet.addChannel(channelPtr);
              continue;
            }
            if(ch == "NaAxon-Traub")
            {
              typedef NaAxonTraub<T,V> Channel;
              Channel* channel = new Channel();
              Dune::shared_ptr<Channel> channelPtr(channel);
              channelSet.addChannel(channelPtr);
              continue;
            }
            if(ch == "K(DR)Axon-Traub")
            {
              typedef K_DRAxonTraub<T,V> Channel;
              Channel* channel = new Channel();
              Dune::shared_ptr<Channel> channelPtr(channel);
              channelSet.addChannel(channelPtr);
              continue;
            }
            if(ch == "NaKPump")
            {
              typedef NaKPump<T,V> Channel;
              Channel* naPump = new Channel(Na);
              Channel* kPump = new Channel(K);
              Dune::shared_ptr<Channel> naPumpPtr(naPump);
              channelSet.addChannel(naPumpPtr);
              Dune::shared_ptr<Channel> kPumpPtr(kPump);
              channelSet.addChannel(kPumpPtr);
              continue;
            }


            DUNE_THROW(Dune::Exception, "No channel type '" << ch << "' could be found!");
          }
        }

      }

      //std::cout << "All channels were added!" << std::endl;
      return channelSet;
    }

    template<class Channel>
    static void createChannel(ChannelSet& channelSet)
    {
      Channel* channel = new Channel();
      Dune::shared_ptr<Channel> channelPtr(channel);
      channelSet.addChannel(channelPtr);
    }

};


#endif /* DUNE_AX1_CHANNEL_BUILDER_HH */
