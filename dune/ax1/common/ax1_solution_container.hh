/*
 * ax1_solution_container.hh
 *
 *  Created on: Sep 17, 2014
 *      Author: jpods
 */

#ifndef DUNE_AX1_SOLUTION_CONTAINER_HH
#define DUNE_AX1_SOLUTION_CONTAINER_HH

template<typename UU, typename GFS_P, typename GFS_C>
class Ax1SolutionContainer
{
  public:
    typedef UU U;
    typedef GFS_P GFS_POT;
    typedef GFS_C GFS_CON;

    Ax1SolutionContainer(const GFS_POT& gfsPot_, const GFS_CON& gfsCon_)
    : gfsPot(gfsPot_),
      gfsCon(gfsCon_)
    {}

    U* getSolutionConOld()
    {
      return ucon_old_ptr;
    }

    U* getSolutionPotOld()
    {
      return upot_old_ptr;
    }

    U* getSolutionConNew()
    {
      return ucon_new_ptr;
    }

    U* getSolutionPotNew()
    {
      return upot_new_ptr;
    }

    U const* getSolutionConOld() const
    {
      return ucon_old_ptr;
    }

    U const* getSolutionPotOld() const
    {
      return upot_old_ptr;
    }

    U const* getSolutionConNew() const
    {
      return ucon_new_ptr;
    }

    U const* getSolutionPotNew() const
    {
      return upot_new_ptr;
    }

    void setSolutionConOld(U* u_ptr)
    {
      ucon_old_ptr = u_ptr;
    }

    void setSolutionPotOld(U* u_ptr)
    {
      upot_old_ptr = u_ptr;
    }

    void setSolutionConNew(U* u_ptr)
    {
      ucon_new_ptr = u_ptr;
    }

    void setSolutionPotNew(U* u_ptr)
    {
      upot_new_ptr = u_ptr;
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
    U* ucon_old_ptr;
    U* ucon_new_ptr;
    U* upot_old_ptr;
    U* upot_new_ptr;
    const GFS_POT& gfsPot;
    const GFS_CON& gfsCon;

};

#endif /* DUNE_AX1_NEWTON_SOLUTION_CONTAINER_HH */
