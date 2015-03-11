#ifndef DUNE_AX1_CHANNEL_HH
#define DUNE_AX1_CHANNEL_HH

#include <cassert>
#include <cmath>
#include <string>

#include <dune/ax1/common/constants.hh>


/** \class   Channel
*  \brief   Base class representing an ion channel following the Hpdgkin-Huxley formalism
*           containing up to two gating particles (for activation/inactivation)
*/
template <class T, class V>
class Channel
{
public:
  typedef Dune::FieldVector<T, NUMBER_OF_SPECIES> VCON;

  Channel () :
    name("none"),
    g(0),
    e_r(0.),
    v_rest(0.),
    temp(6.3),
    temp_fac(1.),
    ionSpecies(-1)
  {}

  Channel (const std::string name_, const T e_r_, const T v_rest_, const T temp_) :
    name(name_),
    g(0),
    e_r(e_r_),
    v_rest(v_rest_),
    temp(temp_),
    temp_fac(pow(3,((temp - 6.3)/10))),
    ionSpecies(-1)
  {}

  //! \brief resize internal vectors to fit  given number of elements
  inline virtual void resize (int nMembElems)
  {
    g.resize(nMembElems);
  }

  //! \brief initialize internal vectors
  virtual void init()
  {
    g = 0.;
  }

  /*//! \brief initialize internal vectors depending on a given potential
  virtual void init(const V& v)
  {
    g = 0.;
  }*/

  //! \brief initialize one component i depending on a given potential; version with concentations, but
  // without dependence on other channels; default implementation does nothing
  virtual void init(const int i, const T& v, const VCON& conCytosol, const VCON& conExtra)
  {
  }

  //! \brief initialize one component i depending on a given potential; version with concentations and
  // dependece on other channels; default implementation does what init(i,v,conCytosol,conExtra) does
  virtual void init(const int i, const T& v, const VCON& conCytosol, const VCON& conExtra,
      const std::map<std::string, T>&)
  {
    init(i,v,conCytosol,conExtra);
  }

  //! \brief get channel conductance at position i
  inline T getConductance (int i) const
  {
    assert(i<this->size());
    return(g[i]);
  }

  //! \brief get channel conductance at position i
  inline void setConductance (int i, const T g_)
  {
    assert(i<this->size());
    g[i] = g_;
  }

  //! \brief get channel reversal potential
  inline T getReversalPotential () const
  {
    return(e_r);
  }

  //! \brief get channel name
  inline std::string getName() const
  {
    return name;
  }

  //! \brief channel info
  inline void info() const
  {
    debug_jochen << "--- Channel " << name << " ---" << std::endl;

    debug_jochen << "###" << std::endl;
    for(int i=0; i<g.size(); ++i)
    {
      debug_jochen << "ME[" << i << "] = " << g[i] << std::endl;
    }
  }

  //! \brief Size of the internal conductance vector
  //! (i.e., the number of membrane elements)
  inline unsigned int size() const
  {
    return g.size();
  }

  inline void setV_rest(const T v_rest_)
  {
    v_rest = v_rest_;
  }

  //! \brief update channel state variables
  virtual void timeStep (int i, const T dt, const T v, const VCON& conCytosol, const VCON& conExtra)
  {
    assert(i<this->size());
    this->step_be(i, dt, v, conCytosol, conExtra);
  }

  //! \brief update state
  virtual void updateState() = 0;

  //! \brief get number of gating particles
  virtual int numGatingParticles () const = 0;

  virtual bool isVoltageGated () const = 0;

  virtual bool isLeakChannel () const = 0;

  virtual bool isConcentrationDependent () const = 0;

  virtual bool isSecondary () const = 0;

  //! \brief get i-th component of the channel's j-th channel state variable
  virtual T getGatingParticle (int j, int i) const = 0;

  //! \brief get i-th component of the channel's j-th channel state variable
  virtual T getNewGatingParticle (int j, int i) const = 0;

  //! \brief get vector of channel's j-th channel state variable
  virtual const V& getGatingParticle (int j) const = 0;

  //! \brief get vector of channel's j-th channel state variable
  virtual const V& getNewGatingParticle (int j) const = 0;

  //! \brief set vector of channel's j-th channel state variable
  virtual void setGatingParticle (int j, const V& p) = 0;

  //! \brief factor for ionic currents at position i
  virtual T getEffConductance (int i) const = 0;

  //! \brief factor for ionic currents at position i
  virtual T getNewEffConductance (int i) const = 0;

  //! \brief update ion concentration this channel is depending on
  virtual void updateConcentration (const V& conc) {}

  //! \brief get the ion species this channel is selective for
  inline int getIonSpecies() const
  {
    return ionSpecies;
  }

protected:

  /**
   * Function as used in NEURON, traps for 0 in denominator
   * of rate equations.
   * @see http://www.neuron.yale.edu/phpBB/viewtopic.php?f=15&t=1075
   *
   * @param x
   * @param y
   * @return
   */
  inline T v_trap(T x, T y) const
  {
    T vtrap;
    if (std::fabs(x/y) < 1e-8) {
      vtrap = y*(1 - x/y/2);
    } else {
      vtrap = x/(std::exp(x/y) - 1);
    }
    return vtrap;
  }

  //! \brief Backward Euler step at position i
  virtual void step_be(const int i, const T dt, const T v, const VCON& conCytosol, const VCON& conExtra) = 0;

  //! \brief Crank-Nicolson step at position i
  virtual void step_cn(const int i, const T dt, const T v, const VCON& conCytosol, const VCON& conExtra) = 0;


  const std::string name;  //! name of channel
  V g;                     //! vector of membrane conductances [1/Ohm/cm^2]
  const T e_r;             //! channel reversal potential [mV]
  T v_rest;                //! the cell's resting potential, important for rate functions
  const T temp;            //! temperature of environment, for scaling of rate functions
  T temp_fac;              //! temperature scaling factor for rate equations according to Nernst equation
  int ionSpecies;          //! Ion species gated by this channel
};

//! \brief abstract base class for channels with voltage-dependent gating particles
template <class T, class V>
class VoltageGatedChannel : public Channel<T,V> {

public:
  typedef Channel<T,V> BaseT;
  typedef typename BaseT::VCON VCON;

  VoltageGatedChannel () :
    Channel<T,V>(),
    nGatingParticles(2),
    p(nGatingParticles),
    p_new(nGatingParticles),
    p_exponents(nGatingParticles)
  {}

  VoltageGatedChannel (const std::string c, int nGatingParticles_, const T e_r_,
      const T v_rest_, const T temp_) :
    Channel<T,V>(c,e_r_,v_rest_,temp_),
    nGatingParticles(nGatingParticles_),
    p(nGatingParticles),
    p_new(nGatingParticles),
    p_exponents(nGatingParticles)
  {}

  //! \brief initialize state variable to steady-state with respect to the channel's resting potential
  virtual inline void init ()
  {
    Channel<T,V>::init();
    for(int j=0; j<nGatingParticles; ++j)
    {
      for(int i=0; i<this->size(); ++i)
      {
        p[j][i] = alpha(this->v_rest, j, i) / (alpha(this->v_rest, j, i) + beta(this->v_rest, j, i));
      }
      p_new[j] = p[j];
    }
  }

//  //! \brief initialize all state variables to steady-state with respect to
//  //! given membrane potential v
//  virtual inline void init (const V &v)
//  {
//    Channel<T,V>::init(v);
//    for(int j=0; j<nGatingParticles; ++j)
//    {
//      for(int i=0; i<this->size(); ++i)
//      {
//        p[j][i] = alpha(v[i], j, i) / (alpha(v[i], j, i) + beta(v[i], j, i));
//      }
//      p_new[j] = p[j];
//    }
//  }

  //! \brief initialize state variable of i-th component to steady-state
  //! with respect to given membrane potential v
  virtual inline void init (const int i, const T& v, const VCON& conCytosol, const VCON& conExtra)
  {
    //debug_jochen << "VoltageGatedChannel::init(i,v,conCytosol,conExtra)" << std::endl;
    for(int j=0; j<nGatingParticles; ++j)
    {
      p[j][i] = alpha(v, j, i) / (alpha(v, j, i) + beta(v, j, i));
      p_new[j] = p[j];
      //debug_jochen << "p[" << j << "][" << i << "] = " << p[j][i] << std::endl;
    }
  }


  //! \brief resize state variable vector
  virtual inline void resize (int nMembElems)
  {
    //std::cout << "Channel_2G.resize" << std::endl;
    for(int j=0; j<nGatingParticles; ++j)
    {
      p[j].resize(nMembElems);
      p_new[j].resize(nMembElems);
    }
    Channel<T,V>::resize(nMembElems);
  }

  //! \brief get number of state variables/gating particles
  virtual int numGatingParticles () const
  {
    return nGatingParticles;
  }

  virtual bool isVoltageGated () const
  {
    return true;
  }

  virtual bool isLeakChannel () const
  {
    return false;
  }

  virtual bool isConcentrationDependent () const
  {
    return false;
  }

  virtual bool isSecondary () const
  {
    return false;
  }

  //! \brief update state
  virtual void updateState()
  {
    for(int j=0; j<nGatingParticles; ++j)
    {
      p[j] = p_new[j];
    }
  }

  virtual inline T getGatingParticle (int j, int i) const
  {
    assert (j<nGatingParticles);
    assert (i<this->size());
    return (p[j][i]);
  }

  //! \brief get channel state vector
  virtual inline T getNewGatingParticle (int j, int i) const
  {
   assert (j<nGatingParticles);
   assert (i<this->size());
   return (p_new[j][i]);
}

  //! \brief effective conductance for membrane element i
  virtual inline T getEffConductance (const int i) const
  {
    assert (i<this->size());

    T g_eff = this->getConductance(i);
    for(int j=0; j<nGatingParticles; ++j)
    {
      g_eff *= std::pow(p[j][i], p_exponents[j]);
    }
    return g_eff;
  }

  //! \brief new effective conductance for membrane element i
  inline T getNewEffConductance (const int i) const
  {
    assert (i<this->size());

    T g_eff = this->getConductance(i);
    for(int j=0; j<nGatingParticles; ++j)
    {
      g_eff *= std::pow(p_new[j][i], p_exponents[j]);
    }
    return g_eff;
  }

  virtual const V& getGatingParticle(int j) const
  {
    assert (j<nGatingParticles);
    return p[j];
  }

  virtual const V& getNewGatingParticle(int j) const
  {
    assert (j<nGatingParticles);
    return p_new[j];
  }

  virtual void setGatingParticle(int j, const V& newState)
  {
    assert (j<nGatingParticles);
    p[j] = newState;
    p_new[j] = newState;
  }


  virtual T alpha (const T v, int j, int i=0)   const = 0;
  virtual T beta  (const T v, int j, int i=0)   const = 0;

private:

  //! \brief do euler step for element i
  inline void step_be (const int i, const T dt, const T v, const VCON& conCytosol, const VCON& conExtra)
  {
    assert (i<this->size());
    for(int j=0; j<nGatingParticles; ++j)
    {
      p_new[j][i] = (dt*alpha(v,j,i) + p[j][i])/
                  (1+dt*(alpha(v,j,i)+ beta(v,j,i)));
    }
  }

  //! \brief do crank-nicolson step for element i
  inline void step_cn (const int i, const T dt, const T v, const VCON& conCytosol, const VCON& conExtra)
  {
    assert (i<this->size());
    // note: formally this step is only 0(t^2) convergent, if v has
    //       been computed at staggered time t+0.5*dt
    for(int j=0; j<nGatingParticles; ++j)
    {
      const T rate = 0.5*dt*(alpha(v,j,i)+ beta(v,j,i));
      p_new[j][i] = (dt*alpha(v,j,i) + (1.0-rate)*p[j][i])/
              (1.0+rate);
    }
  }

  const int nGatingParticles;
  std::vector<V> p;
  std::vector<V> p_new;

protected:
  std::vector<int> p_exponents;
};


//! \brief A leak channel as in Hodgkin's and Huxley's model
template <class T, class V>
class LeakChannel : public Channel<T,V> {

public:

  typedef Channel<T,V> BaseT;
  typedef typename BaseT::VCON VCON;

  // constructor
  LeakChannel (const int ionSpecies_)
    : Channel<T,V>::Channel(ION_NAMES[ionSpecies_] + "_leak", 0., 0., 6.3),
      g_m(0.5),
      ratio()
  {
    ratio = 1.0;
    this->ionSpecies = ionSpecies_;
  }

  // constructor
  LeakChannel (const T g_m_, const int ionSpecies_, const T temp_)
    : Channel<T,V>::Channel(ION_NAMES[ionSpecies_] + "_leak", 0.0, 0.0, temp_),
      g_m(g_m_),
      ratio()
  {
    ratio = 1.0;
    this->ionSpecies = ionSpecies_;
  }

  //! \brief initialize state variable to steady-state
  virtual inline void init ()
  {
    Channel<T,V>::init();
    this->g = g_m;
  }

//  //! \brief initialize state variable to steady-state
//  virtual inline void init (const V& v)
//  {
//    Channel<T,V>::init(v);
//    this->g = g_m;
//  }

  //! \brief resize state variable vector
  virtual inline void resize (const int nMembElems)
  {
    ratio.resize(nMembElems);
    Channel<T,V>::resize(nMembElems);
  }

  //! \brief get number of state variables/gating particles
  virtual int numGatingParticles () const
  {
    return 1;
  }

  virtual bool isVoltageGated () const
  {
    return false;
  }

  virtual bool isLeakChannel () const
  {
    return true;
  }

  virtual bool isConcentrationDependent () const
  {
    return false;
  }

  virtual bool isSecondary () const
  {
    return false;
  }

  //! \brief update state; does nothing as leak channels have a fixed conductance
  virtual void updateState()
  {
  }

  //! \brief get total leak conductance for membrane element i
  virtual inline T getConductance (const int i) const
  {
    assert (i<this->size());
    return(this->g[i]);
  }

  //! \brief effective conductance for membrane element i
  virtual inline T getEffConductance (const int i) const
  {
    assert (i<this->size());
    return(this->g[i] * ratio[i]);
  }

  //! \brief new effective conductance for element i (= ratio * max coductance)
  virtual inline T getNewEffConductance (const int i) const
  {
    assert (i<this->size());
    return(this->g[i] * ratio[i]);
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
    assert (j==0);
    assert (i<this->size());
    return ratio[i];
  }

  //! \brief A little awkward: This leak channel's ratio is defined as a gating particle
  virtual T getNewGatingParticle (int j, int i) const
  {
    return getGatingParticle(j, i);
  }

  virtual const V& getGatingParticle(int j) const
  {
    assert(j == 0);
    return ratio;
  }

  virtual const V& getNewGatingParticle(int j) const
  {
    return getGatingParticle(j);
  }

  virtual void setGatingParticle(int j, const V& newState)
  {
    assert(j == 0);
    ratio = newState;
  }

  // ---------------- Dummy methods to comply with interface -----------------------------
  //! \brief Backward Euler step at position i
  virtual void step_be(const int i, const T dt, const T v, const VCON& conCytosol, const VCON& conExtra)
  {
  }

  //! \brief Crank-Nicolson step at position i
  virtual void step_cn(const int i, const T dt, const T v, const VCON& conCytosol, const VCON& conExtra)
  {
  }
  // -------------------------------------------------------------------------------------

private:
  const T g_m;
  V ratio;
};

// A natrium/sodium channel as in Hodgkin's and Huxley's model
template <class T, class V>
class NavChannel : public VoltageGatedChannel<T,V> {

public:

  NavChannel () :
    //VoltageGatedChannel<T,V>("Nav", 2, 115., 0. , 6.3)
    VoltageGatedChannel<T,V>("Nav", 2, 50., -65. , 6.3)
  {
    this->ionSpecies = Na;
    this->p_exponents[0] = 3;
    this->p_exponents[1] = 1;
  }

  //NavChannel (const std::string c, const T e_r, const T v_rest_ = 0., const T temp_ = 6.3) :
  NavChannel (const std::string c, const T e_r, const T v_rest_ = -65., const T temp_ = 6.3) :
    VoltageGatedChannel<T,V>(c, 2, e_r, v_rest_, temp_)
  {
    this->ionSpecies = Na;
  }

  virtual inline T alpha (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<2);

    // activation rate
    if(j==0)
    {
      const T v_ = v-this->v_rest;
      const T a = 0.1*this->v_trap((25-v_),10);
      return(this->temp_fac * a);
    // inactivation rate
    } else {
      const T v_ = v-this->v_rest;
      const T a = 0.07*std::exp(-v_/20);
      return(this->temp_fac * a);
    }
  }

  virtual inline T beta (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<2);

    // activation rate
    if(j==0)
    {
      const T v_ = v-this->v_rest;
      const T b = 4*std::exp(-v_/18);
      return(this->temp_fac * b);
    // inactivation rate
    } else {
      const T v_ = v-this->v_rest;
      const T b = 1./(std::exp((30-v_)/10)+1);
      return(this->temp_fac * b);
    }
  }

};

//! \brief A potassium channel as in Hodgkin's and Huxley's model
template <class T, class V>
class KvChannel : public VoltageGatedChannel<T,V> {

public:

  KvChannel () :
    //VoltageGatedChannel<T,V>("Kv", 1, -12., 0., 6.3)
    VoltageGatedChannel<T,V>("Kv", 1, -77., -65., 6.3)
  {
    this->ionSpecies= K;
    this->p_exponents[0] = 4;
  }

  //KvChannel (std::string c, T e_r, const T v_rest_ = 0., const T temp_ = 6.3) :
  KvChannel (std::string c, T e_r, const T v_rest_ = -65., const T temp_ = 6.3) :
    VoltageGatedChannel<T,V>(c, 1, e_r, v_rest_, temp_)
  {
    this->ionSpecies= K;
  }

  virtual inline T alpha (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j==0);

    const T v_ = v-this->v_rest;
    const T a = 0.01 * this->v_trap((10-v_),10);

    return(this->temp_fac * a);
  }

  virtual inline T beta (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j==0);

    const T v_ = v-this->v_rest;
    const T b = 0.125*std::exp(-v_/80);

    return(this->temp_fac * b);
  }
};

#endif /* DUNE_AX1_CHANNEL_HH */
