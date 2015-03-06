/*
 * ax1_newton_solution_container.hh
 *
 *  Created on: Sep 17, 2014
 *      Author: jpods
 */

#ifndef DUNE_AX1_NEWTON_SOLUTION_CONTAINER_HH
#define DUNE_AX1_NEWTON_SOLUTION_CONTAINER_HH

template<typename UU, typename GFS_P, typename GFS_C>
class Ax1NewtonSolutionContainer
{
  public:
    typedef UU U;
    typedef GFS_P GFS_POT;
    typedef GFS_C GFS_CON;

    Ax1NewtonSolutionContainer(const GFS_POT& gfsPot_, const GFS_CON& gfsCon_)
    : gfsPot(gfsPot_),
      gfsCon(gfsCon_)
    {}

    U* getSolution()
    {
      return u_ptr;
    }

    void setSolution(U* new_u_ptr)
    {
      u_ptr = new_u_ptr;
    }

    U* getSolutionCon()
    {
      return u_ptr;
    }

    U* getSolutionPot()
    {
      return u_ptr;
    }

    U const* getSolutionCon() const
    {
      return u_ptr;
    }

    U const* getSolutionPot() const
    {
      return u_ptr;
    }

    const GFS_POT& getGfsPot() const
    {
      return gfsPot;
    }

    const GFS_CON& getGfsCon() const
    {
      return gfsCon;
    }

  private:
    U* u_ptr;
    const GFS_POT& gfsPot;
    const GFS_CON& gfsCon;

};

#endif /* DUNE_AX1_NEWTON_SOLUTION_CONTAINER_HH */
