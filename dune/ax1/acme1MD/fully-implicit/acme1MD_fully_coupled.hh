/*
 * acme1MD_fully_coupled.hh
 *
 *  Created on: Dec 15, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME1MD_FULLY_COUPLED_HH
#define DUNE_AX1_ACME1MD_FULLY_COUPLED_HH

#include <dune/ax1/common/tools.hh>

class Acme1MDFullyCoupled
{
  public:

    template<typename U, typename U_CON, typename U_POT,
             typename SOLVER, typename GFS,
             typename DIAG_INFO, typename Real>
    static void timeStep(Real time, Real dt, U& uold, U& unew, SOLVER& solver,
        GFS& gfs, U_CON& unewCon, U_POT& uPot, DIAG_INFO& diagInfo)
    {
      //debug_verb << "================= Solving fully implicit system..." << std::endl;
      solver.apply(time,dt,uold,unew);
      // The following two lines are essential and must NOT be removed!
      // Grid functions (e.g. charge density) rely on the current value of child GFS coefficient vectors
      Tools::compositeToChildCoefficientVector(gfs, unew, unewCon, 0);
      Tools::compositeToChildCoefficientVector(gfs, unew, uPot, 1);
      //debug_verb << "================= Solved fully implicit system" << std::endl;
      // Fill diagnosticInfo
      diagInfo.iterations = solver.getPDESolver().result().iterations;
      diagInfo.dt = dt;
    }
};


#endif /* DUNE_AX1_ACME1MD_FULLY_COUPLED_HH */
