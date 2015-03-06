#ifndef DUNE_AX1_CHANNELGRIDFUNCTION_HH
#define DUNE_AX1_CHANNELGRIDFUNCTION_HH

#include <dune/ax1/common/constants.hh>


template<typename GV, typename RF, typename PHYSICS>
class ChannelGridFunction
  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::IntersectionGridFunctionTraits<GV,
                                       RF,
                                       5, // This should be enough for every channel (gbar, g and up to 3 gating particles)
                                       Dune::FieldVector<RF,5> >,
      ChannelGridFunction<GV,RF,PHYSICS> >
{
public:
  typedef Dune::PDELab::IntersectionGridFunctionTraits<GV,   // grid view type
                                           RF, // image field type (double)
                                           5,                                        // number of components of image (k)
                                           Dune::FieldVector<RF,5> // image type (Dune::FieldVector<double, k>)
                                           > Traits;

  typedef typename PHYSICS::ChannelSet::ChannelType Channel;

  //! constructor
  ChannelGridFunction (const GV& gv_, PHYSICS& physics_, int channelIndex_)
    : gridView(gv_),
      physics(physics_),
      channelIndex(channelIndex_),
      channel(physics.getMembrane().getChannelSet().getChannel(channelIndex_)),
      membraneIndices(physics.getChannelSet().getMembraneIndices())
  {
  }

  template<typename I>
  inline void evaluate (const I &is,
                        const typename Traits::DomainType &x,
                        typename Traits::RangeType &y) const
  {
    y = 0.0;

    if(physics.isMembraneInterface(is))
    {
      int iIndex = physics.getIntersectionIndex(is);
      int mIndex = membraneIndices.at(iIndex);

      // maximum conductance (HH: 'gbar')
      y[0] = channel.getConductance(mIndex);

      // actual conductance based on mac conductance and current state (gating particles)
      y[1] = channel.getEffConductance(mIndex);

      // This gridfunction has a vector of size 5 as image type, therefore no more than 3 gating particles
      // can be represented.
      assert(channel.numGatingParticles() <= 3);

      // Here come the gating particles
      for(int l=0; l<channel.numGatingParticles(); l++)
      {
        y[2+l] = channel.getGatingParticle(l, mIndex);
      }
    }
  }



  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  const typename Traits::GridViewType& gridView;
  PHYSICS& physics;
  const int channelIndex;

  const Channel& channel;
  const std::map<int,int>& membraneIndices;
};

#endif /* DUNE_AX1_CHANNELGRIDFUNCTION_HH */
