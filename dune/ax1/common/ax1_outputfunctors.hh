/*
 * ax1_outputfunctors.hh
 *
 *  Created on: May 6, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_OUTPUTFUNCTORS_HH
#define DUNE_AX1_OUTPUTFUNCTORS_HH

struct Ax1OutputFunctors
{

  struct IdentityFunctor
  {
    template<typename V>
    V operator() (const V& v) const
    {
      return v;
    }
  };

  template<typename PHYSICS>
  class to_millivolt
  {
    public:
     to_millivolt(PHYSICS& physics_)
     : physics(physics_)
     {}

     template<typename V>
     V operator() (const V& v) const
     {
       return physics.convertTo_mV(v);
     }

    private:
     const PHYSICS& physics;
  };

};


#endif /* DUNE_AX1_OUTPUTFUNCTORS_HH */
