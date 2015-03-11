#ifndef DUNE_AX1_TRAUB_CHANNELS_HH
#define DUNE_AX1_TRAUB_CHANNELS_HH

#include <dune/ax1/channels/channel.hh>

/**
 * traub_channels.hh
 *
 * \brief All channels implemented according to Traub (1994)
 *
 *  Created on: Jun 27, 2011
 *      Author: jpods
 */


/**************************** SD compartments **************************/

//! Sodium channel for SD compartments
template <class T, class V>
class NaTraub : public VoltageGatedChannel<T,V> {

public:

  NaTraub () :
    VoltageGatedChannel<T,V>("Na-Traub", 2, 115., 0. , 23)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 3;
    this->p_exponents[1] = 1;
  }

  NaTraub (const std::string c, const T e_r, const T v_rest_ = 0., const T temp_ = 23) :
    VoltageGatedChannel<T,V>(c, 2, e_r, v_rest_, temp_)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 3;
    this->p_exponents[1] = 1;
  }

  virtual inline T alpha (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    if(j==0)
    {
      const T a = 0.32*this->v_trap((13.1-v_),4);
      return(this->temp_fac * a);
    } else
    {
      const T a = 0.128*std::exp((17-v_)/18);
      return(this->temp_fac * a);
    }
  }

  virtual inline T beta (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    if(j==0)
    {
      const T b = 0.28*this->v_trap((v_-40.1),5);
      return(this->temp_fac * b);
    } else {
      const T b = 4./(std::exp((40-v_)/5)+1);
      return(this->temp_fac * b);
    }
  }
};


//! High-threshold calcium channel
template <class T, class V>
class CaTraub : public VoltageGatedChannel<T,V> {

public:

  // constructor
  CaTraub () :
    VoltageGatedChannel<T,V>("Ca-Traub", 1, 140., 0., 23)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 2;
  }

  // constructor
  CaTraub (std::string c, T e_r, const T v_rest_ = 0., const T temp_ = 23) :
    VoltageGatedChannel<T,V>(c, 1, e_r, v_rest_, temp_)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 2;
  }

  // rate functions for activation
  virtual inline T alpha (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    const T a = 1.6 / (1 + std::exp(-0.072*(v_-65)));

    return(this->temp_fac * a);
  }

  virtual inline T beta (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    const T b = 0.02*this->v_trap(v_-51.1, 5);

    return(this->temp_fac * b);
  }
};


//! High-threshold calcium channel (1991 version)
template <class T, class V>
class CaTraub1991 : public VoltageGatedChannel<T,V> {

public:

  CaTraub1991 () :
    VoltageGatedChannel<T,V>("Ca-Traub1991", 2, 140., 0., 23)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 2;
    this->p_exponents[1] = 1;
  }

  CaTraub1991 (const std::string c, const T e_r, const T v_rest_ = 0., const T temp_ = 23)
    : VoltageGatedChannel<T,V>(c, 2, e_r, v_rest_, temp_)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 2;
    this->p_exponents[1] = 1;
  }

  virtual inline T alpha (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    if(j==0)
    {
      const T a = 1.6 / (1 + std::exp(-0.072*(v_-65)) );
      return(this->temp_fac * a);
    } else {
      T a = 0.005;
      if(v_ > 0)
      {
        a = std::exp(-v_/20)/200.;
      }
      return(this->temp_fac * a);
    }
  }

  virtual inline T beta (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    if(j==0)
    {
      const T b = 0.02*this->v_trap(v_-51.1, 5);
      return(this->temp_fac * b);
    } else {
      T b = 0.;
      if(v_ > 0)
      {
        b = 0.005 - this->alpha(v,j,i);
      }
      return(this->temp_fac * b);
    }
  }
};

//! Delayed-rectifier potassium channel
template <class T, class V>
class K_DRTraub : public VoltageGatedChannel<T,V> {

public:

  K_DRTraub () :
    VoltageGatedChannel<T,V>("K(DR)-Traub", 1, -15., 0., 23)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 2;
  }

  K_DRTraub (std::string c, T e_r, const T v_rest_ = 0., const T temp_ = 23) :
    VoltageGatedChannel<T,V>(c, 1, e_r, v_rest_, temp_)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 2;
  }

  virtual inline T alpha (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    const T a = 0.016 * this->v_trap(35.1-v_, 5);

    return(this->temp_fac * a);
  }

  virtual inline T beta (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    const T b = 0.25 * std::exp((20-v_) / 40);

    return(this->temp_fac * b);
  }
};


//! Delayed-rectifier potassium channel (1991 version)
template <class T, class V>
class K_DRTraub1991 : public VoltageGatedChannel<T,V> {

public:

  // constructor
  K_DRTraub1991 ()
    : VoltageGatedChannel<T,V>("K(DR)-Traub1991", 1, -15., 0., 23)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 1;
  }

  // constructor
  K_DRTraub1991 (std::string c, T e_r, const T v_rest_ = 0., const T temp_ = 23)
    : VoltageGatedChannel<T,V>(c, 1, e_r, v_rest_, temp_)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 1;
  }

  virtual inline T alpha (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    const T a = 0.016 * this->v_trap(35.1-v_, 5);

    return(this->temp_fac * a);
  }

  virtual inline T beta (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    const T b = 0.25 * std::exp((20-v_) / 40);

    return(this->temp_fac * b);
  }
};


//! slow AHP potassium channel
template <class T, class V>
class K_AHPTraub : public VoltageGatedChannel<T,V> {

public:

  // constructor
  K_AHPTraub () :
    VoltageGatedChannel<T,V>("K(AHP)-Traub", 1, -15., 0., 23),
    chi(0.)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 1;
  }

  // constructor
  K_AHPTraub (std::string c, T e_r, const T v_rest_ = 0., const T temp_ = 23) :
    VoltageGatedChannel<T,V>::Channel_1G(c, 1, e_r, v_rest_, temp_),
    chi(0.)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 1;
  }

  virtual inline T alpha (const T v, int j, int i) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T a = std::min(0.2 * 1e-4 * (T) chi[i], 0.01);
    return(this->temp_fac * a);
  }

  virtual inline T beta (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T b = 0.001;
    return(this->temp_fac * b);
  }

  inline void resize(int nMembElems)
  {
    chi.resize(nMembElems);
    VoltageGatedChannel<T,V>::resize(nMembElems);
  }

  virtual inline void init () {
    chi = 0.;
    VoltageGatedChannel<T,V>::init();
  }

//  virtual inline void init (const V &v) {
//    chi = 0.;
//    VoltageGatedChannel<T,V>::init(v);
//  }

  inline void updateConcentration(const V& chi_)
  {
    chi = chi_;
  }

private:
  V chi;
};


//! Transient ('A') potassium channel
template <class T, class V>
class K_ATraub : public VoltageGatedChannel<T,V> {

public:

  K_ATraub ()
    : VoltageGatedChannel<T,V>("K(A)-Traub", 2, -15., 0., 23)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 1;
    this->p_exponents[1] = 1;
  }

  K_ATraub (const std::string c, const T e_r, const T v_rest_ = 0., const T temp_ = 23)
    : VoltageGatedChannel<T,V>(c, 2, e_r, v_rest_, temp_)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 1;
    this->p_exponents[1] = 1;
  }

  virtual inline T alpha (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    if(j==0)
    {
      const T a = 0.02*this->v_trap((13.1-v_),10);
      return(this->temp_fac * a);
    } else {
      const T a = 0.0016*std::exp((-13-v_)/18);
      return(this->temp_fac * a);
    }
  }

  virtual inline T beta (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    if(j==0)
    {
      const T b = 0.0175*this->v_trap((v_-40.1),10);
      return(this->temp_fac * b);
    } else {
      const T b = 0.05/(std::exp((10.1-v_)/5)+1);
      return(this->temp_fac * b);
    }
  }
};


//! Rapid voltage- and Ca2+-dependent ('C') potassium channel
template <class T, class V>
class K_CTraub : public VoltageGatedChannel<T,V> {

public:

  K_CTraub () :
    VoltageGatedChannel<T,V>("K(C)-Traub", 1, -15., 0., 23)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 1;
  }

  K_CTraub (std::string c, T e_r, const T v_rest_ = 0., const T temp_ = 23) :
    VoltageGatedChannel<T,V>(c, 1, e_r, v_rest_, temp_)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 1;
  }

  //! \brief effective conductance for membrane element i
  virtual inline T getEffConductance (const int i) const
  {
    assert (i<this->size());

    const T g = this->getConductance(i);
    const T c = this->getGatingParticle(0,i);
    const T g_eff = g*c*std::min(1.0, (T) chi[i]/250.);

    return g_eff;
  }

  //! \brief new effective conductance for membrane element i
  inline T getNewEffConductance (const int i) const
  {
    assert (i<this->size());

    const T g = this->getConductance(i);
    const T c = this->getNewGatingParticle(0,i);
    const T g_eff_new = g*c*std::min(1.0, (T) chi[i]/250.);

    return g_eff_new;
  }

  virtual inline T alpha (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    T alpha_c = 0.;
    if(v_ <= 50)
    {
      alpha_c = std::exp( (((v_-10)/11) - ((v_-6.5)/27)) / 18.975 );
    } else {
      alpha_c = 2*std::exp(-(v_-6.5)/27);
    }
    return(this->temp_fac * alpha_c);
  }

  virtual inline T beta (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    T beta_c = 0.;
    if(v_ <= 50)
    {
      beta_c = 2 * std::exp(-(v_-6.5)/27) - this->alpha(v,j,i);
    }
    return (this->temp_fac * beta_c);
  }

  inline void resize(int nMembElems)
  {
    chi.resize(nMembElems);
    VoltageGatedChannel<T,V>::resize(nMembElems);
  }

  virtual inline void init () {
    chi = 0.;
    VoltageGatedChannel<T,V>::init();
  }

//  virtual inline void init (const V &v) {
//    chi = 0.;
//    VoltageGatedChannel<T,V>::init(v);
//  }

  inline void updateConcentration(const V& chi_)
  {
    chi = chi_;
  }

private:
  V chi;
};


/**************************** Axon compartments **************************/

//! Sodium channel for axonic compartments
template <class T, class V>
class NaAxonTraub : public VoltageGatedChannel<T,V> {

public:

  // constructor
  NaAxonTraub () :
    VoltageGatedChannel<T,V>("NaAxon-Traub", 2, 115., 0. , 23)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 3;
    this->p_exponents[1] = 1;
  }

  // constructor
  NaAxonTraub (const std::string c, const T e_r, const T v_rest_ = 0.,
      const T temp_ = 23) :
    VoltageGatedChannel<T,V>(c, 2, e_r, v_rest_, temp_)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 3;
    this->p_exponents[1] = 1;
  }

  virtual inline T alpha (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    if(j==0)
    {
      const T a = 0.8*this->v_trap((17.2-v_),4);
      return(this->temp_fac * a);
    } else {
      const T a = 0.32*std::exp((42-v_)/18);
      return(this->temp_fac * a);
    }
  }

  virtual inline T beta (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    if(j==0)
    {
      const T b = 0.7*this->v_trap((v_-42.2),5);
      return(this->temp_fac * b);
    } else {
      const T b = 10./(std::exp((42-v_)/5)+1);
      return(this->temp_fac * b);
    }
  }
};


//! Delayed-rectifier potassium channel for axonic compartments
template <class T, class V>
class K_DRAxonTraub : public VoltageGatedChannel<T,V> {

public:

  K_DRAxonTraub () :
    VoltageGatedChannel<T,V>("K(DR)Axon-Traub", 1, -25., 0., 23)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 4;
  }

  K_DRAxonTraub (std::string c, T e_r, const T v_rest_ = 0., const T temp_ = 23) :
    VoltageGatedChannel<T,V>(c, 1, e_r, v_rest_, temp_)
  {
    this->temp_fac = pow(2.3,((this->temp - 23)/10));
    this->p_exponents[0] = 4;
  }

  virtual inline T alpha (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    const T a = 0.03 * this->v_trap(17.2-v_, 5);

    return(this->temp_fac * a);
  }

  virtual inline T beta (const T v, int j, int i=0) const
  {
    assert (i<this->size());
    assert (j<this->numGatingParticles());

    const T v_ = v-this->v_rest;
    const T b = 0.45 * std::exp((12-v_)/40);

    return(this->temp_fac * b);
  }
};

#endif /* DUNE_AX1_TRAUB_CHANNELS_HH */
