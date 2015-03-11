/*
 * ion_pump.hh
 *
 *  Created on: Nov 28, 2013
 *      Author: jpods
 */

#include <dune/ax1/channels/channel.hh>

#ifndef DUNE_AX1_ION_PUMP_HH
#define DUNE_AX1_ION_PUMP_HH

//! \brief  A sodium-potassium ion pump implemented according to Scriven 1981; use two channels, one for
// potassium with constant ratio = 1, and one for sodium with ratio = -1.5 (and dynamically changing)
template <class T, class V>
class NaKPump : public Channel<T,V> {

public:
  typedef Channel<T,V> BaseT;
  typedef typename BaseT::VCON VCON;

  // constructor
  NaKPump (const int ionSpecies_)
    : Channel<T,V>::Channel(ION_NAMES[ionSpecies_] + "_pump", 0., 0., 6.3),
      a(0.0),
      b1(30.),
      b2(1.),
      c(0.05),
      d(0.0),
      ratio(),
      ratio_new()
  {
    ratio = ionSpecies_ == Na ? 1.5 : 1.0;
    ratio_new = ratio;
    this->ionSpecies = ionSpecies_;
  }

  // constructor
  NaKPump (const T g_m_, const int ionSpecies_, const T temp_)
    : Channel<T,V>::Channel(ION_NAMES[ionSpecies_] + "_pump", 0.0, 0.0, temp_),
      ratio(),
      ratio_new()
  {
    ratio = ionSpecies_ == Na ? 1.5 : 1.0;
    ratio_new = ratio;
    this->ionSpecies = ionSpecies_;
  }

  //! \brief initialize state variable of i-th component to steady-state
  //! with respect to given membrane potential v
  virtual inline void init (const int i, const T& v, const VCON& conCytosol, const VCON& conExtra,
      const std::map<std::string, T>& gOtherChannels)
  {
    T kExt = conExtra[K];
    T naIn = conCytosol[Na];

    // Calculate d such that r = 1.5 at rest
    d = 1.5 - c*naIn;

    // TODO Bring I_K to the correct units
    // Calculate a such that I_K = -I_Kp at rest
    T I_K = (gOtherChannels.at("Kv") + gOtherChannels.at("K_leak")) * v;
    a = I_K * ( (1+(this->b1 / kExt))*(1+(this->b1 / kExt))*(1+(this->b2 / naIn)) );
  }

  //! \brief resize state variable vector
  virtual inline void resize (const int nMembElems)
  {
    ratio.resize(nMembElems);
    ratio_new.resize(nMembElems);
    I_Kp.resize(nMembElems);
    I_Kp_new.resize(nMembElems);
    Channel<T,V>::resize(nMembElems);
  }

  //! \brief get number of state variables/gating particles
  virtual int numGatingParticles () const
  {
    return 2;
  }

  virtual bool isVoltageGated () const
  {
    return false;
  }

  virtual bool isLeakChannel () const
  {
    return false;
  }

  virtual bool isConcentrationDependent () const
  {
    return true;
  }

  virtual bool isSecondary () const
  {
    return true;
  }

  virtual void updateState()
  {
    ratio = ratio_new;
    I_Kp = I_Kp_new;
  }

  //! \brief get conductance for membrane element i
  virtual inline T getConductance (const int i) const
  {
    assert (i<this->size());
    return(this->g[i]);
  }

  //! \brief effective conductance for membrane element i
  virtual inline T getEffConductance (const int i) const
  {
    assert (i<this->size());
    return(this->g[i] * ratio[i] * I_Kp[i]);
  }

  //! \brief new effective conductance for element i (= ratio * max coductance)
  virtual inline T getNewEffConductance (const int i) const
  {
    assert (i<this->size());
    return(this->g[i] * ratio_new[i] * I_Kp_new[i]);
  }

  //! This is the ratio of the maximum conductance which is active
  virtual inline void setRatio(T ratio_, const int i)
  {
    assert (i<this->size());
    ratio[i] = ratio_;
  }

  //! \brief A little awkward: This leak channel's ratio is defined as a gating particle
  virtual T getGatingParticle (int j, int i) const
  {
    assert (j<=1);
    assert (i<this->size());
    if(j==0)
      return ratio[i];
    else
      return I_Kp[i];
  }

  //! \brief A little awkward: This leak channel's ratio is defined as a gating particle
  virtual T getNewGatingParticle (int j, int i) const
  {
    assert (j<=1);
    assert (i<this->size());
    if(j==0)
      return ratio[i];
    else
      return I_Kp_new[i];
  }

  virtual const V& getGatingParticle(int j) const
  {
    assert(j <= 1);
    if(j==0)
      return ratio;
    else
      return I_Kp;
  }

  virtual const V& getNewGatingParticle(int j) const
  {
    assert(j <= 1);
    if(j==0)
      return ratio;
    else
      return I_Kp_new;
  }

  virtual void setGatingParticle(int j, const V& newState)
  {
    assert(j <= 1);
    if(j==0)
      ratio = newState;
    else
      I_Kp = newState;
  }

  // ---------------- Dummy methods to comply with interface -----------------------------
  //! \brief Backward Euler step at position i
  virtual void step_be(const int i, const T dt, const T v, const VCON& conCytosol, const VCON& conExtra)
  {
    T kExt = conExtra[K];
    T naIn = conCytosol[Na];

    // Ratio is dynamic for sodium; for K, it is constant 1
    if(this->ionSpecies == Na)
    {
      ratio_new[i] = this->c * naIn + this->d;
    }

    // TODO Bring I_Kp_new to the correct units
    // Calculate new 'conductance' which actually is a current
    I_Kp_new[i] = this->a / ( (1+(this->b1 / kExt))*(1+(this->b1 / kExt))*(1+(this->b2 / naIn)) );

    // Hack: We actually calculate a _current_, but the interface expects a _conductance_ that is later
    // multiplied by the membrane voltage to yield a current. So we divide by the membrane voltage here
    // such that after multiplication outside the desired current comes out again
    I_Kp_new[i] /= v;
  }

  //! \brief Crank-Nicolson step at position i
  virtual void step_cn(const int i, const T dt, const T v, const VCON& conCytosol, const VCON& conExtra)
  {
    // TODO Calculate new 'conductance'
    DUNE_THROW(Dune::NotImplemented, "Crank-Nicholson not implemented for class IonPump");
  }
  // -------------------------------------------------------------------------------------

private:
  T a, b1, b2, c, d; // Parameters according to Scriven 1981
  V ratio;
  V ratio_new;
  V I_Kp;
  V I_Kp_new;
};






#endif /* DUNE_AX1_ION_PUMP_HH */
