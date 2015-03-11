/*
 * acme2_fully_coupled.hh
 *
 *  Created on: Dec 15, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2_FULLY_COUPLED_HH
#define DUNE_AX1_ACME2_FULLY_COUPLED_HH

#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/ax1_boundaryfunction_membranefunction_adapter.hh>

template<typename Traits,typename ACME2_OUTPUT>
class Acme2Simulation
{
  public:

    typedef typename Traits::Real Real;

    template<typename SOLVER>
    static void run(Real& time, Real& dt, Real& dtstart, Real& tend, const Real& tEquilibrium,
        typename Traits::Physics& physics,
        const typename Traits::GridView& gv, const typename Traits::SubGridView& membGV,
        SOLVER& solver,
        typename Traits::PARAMETERS_CON& parametersCon,
        typename Traits::MultiGFS& multigfs,
        typename Traits::U& uold, typename Traits::U& unew,
        //typename Traits::U_CON& uoldCon, typename Traits::U_CON& unewCon,
        //typename Traits::U_POT& uoldPot, typename Traits::U_POT& unewPot,
        typename Traits::GF_MEMB_FLUX& gfMembFlux,
        typename Traits::DGF_CON& dgfCon,
        typename Traits::DGF_POT& dgfPot, typename Traits::DGF_POT_GRAD& dgfGradPot,
        ACME2_OUTPUT& acme2Output,
        Acme2SolutionVectors<Traits>& solutionVectors)
    {
      const Acme2Parameters& params = physics.getParams();

      typename ACME2_OUTPUT::DiagnosticInfo& diagInfo = acme2Output.getDiagInfo();
      typename ACME2_OUTPUT::DiagnosticInfo lastDiagInfo(diagInfo); // copy of diagInfo


      Real dpot(0.0);
      Real old_dt = dt;
      Real last_dt = dt;
      Real last_dpot(0.0);

      diagInfo.tEquilibrium = tEquilibrium;
      diagInfo.registerDebugData(std::string("abs max pot change"), 0.0);
      diagInfo.registerDebugData(std::string("abs max con change"), 0.0);
      diagInfo.registerDebugData(std::string("rel max pot change"), 0.0);
      diagInfo.registerDebugData(std::string("rel max con change"), 0.0);

      diagInfo.registerDebugData(std::string("old_dt"), old_dt);
      diagInfo.registerDebugData(std::string("last_dt"), last_dt);
      diagInfo.registerDebugData(std::string("dpot"), dpot);
      diagInfo.registerDebugData(std::string("last_dpot"), last_dpot);
      diagInfo.registerDebugData(std::string("dpot_dt"), 0.0);
      diagInfo.registerDebugData(std::string("last_dpot_dt"), 0.0);
      diagInfo.registerDebugData(std::string("rel_dpot_dt"), 0.0);

      diagInfo.registerDebugData(std::string("flux_factor"), 1.0);

      // This creates the file with all the header information from the
      // debug data dynamically added via 'registerDebugData()'
      acme2Output.initDiagInfoFile();

      typedef Ax1GoldmanEquationGridFunction<typename Traits::SubGridView,typename Traits::DGF_CON,
          typename Traits::Physics> GoldmanGF;
      GoldmanGF goldmanGF(membGV,dgfCon,physics);


      // Don't do initial output when a previous simulation is continued
      if(not (physics.getParams().doLoadState() && physics.getParams().doContinueSimulation()))
      {
        acme2Output.writeStep(time);
        Tools::pecletNumber(gv, parametersCon, dgfGradPot);

        typename Traits::DGF_POT::Traits::RangeType minMembPot = 1e100;
        typename Traits::DGF_POT::Traits::RangeType maxMembPot = -1e100;

        for(typename Traits::Physics::SubDomainElementIterator sdeit = membGV.template begin<0>();
            sdeit != membGV.template end<0>(); ++sdeit)
        {
          typename Traits::DGF_POT::Traits::RangeType membPot;
          physics.getMembranePotential(membGV.grid().multiDomainEntity(*sdeit), dgfPot, membPot);

          if(membPot < minMembPot) minMembPot = membPot;
          if(membPot > maxMembPot) maxMembPot = membPot;
        }
        debug_info << "MAX membrane potential = " << physics.convertTo_mV(maxMembPot) << " mV" << std::endl;
        debug_info << "MIN membrane potential = " << physics.convertTo_mV(minMembPot) << " mV" << std::endl;

        debug_info  << std::endl << "########## initial output done" << std::endl << std::endl;
      }
      // ==============================================================================================


      // ========= time loop ==========================================================================
      Real outputTimeInterval = physics.getParams().getOutputTimeInterval();
      int outputCounter = 1;
      bool printEveryTimeStep = (dt > outputTimeInterval && !physics.getParams().useAdaptiveTimeStep());

      double tInj_start = physics.getParams().tInj_start();
      double tInj_end   = physics.getParams().tInj_end();

      debug_jochen << "tInj_start = " << tInj_start << ", tInj_end = " << tInj_end << std::endl;

      if(physics.getParams().doStimulation())
      {
        assert(tInj_start > tEquilibrium);
      }

      typename Traits::U uChange = unew;
      //typename Traits::U_CON uConChange(unewCon);
      //typename Traits::U_POT uPotChange(unewPot);
      //typename Traits::U_CON uConChange_Rel(unewCon);
      //typename Traits::U_POT uPotChange_Rel(unewPot);

      // Total number of time steps
      int totalTimeSteps = 0;

      // Number of iterations/time steps since last output
      int iterations = 0;
      int timeSteps = 0;

      while (time<tend-1e-8)
      {
        debug_verb << "TIME = " << time << std::endl;
        debug_verb << "Last dt = " << dt << std::endl;

        last_dt = old_dt;
        old_dt = dt;

        // Use default time step when equilibration phase is over
        if(std::abs(time-tEquilibrium) < 1e-6)
        {
          dt = dtstart;
          bool printEveryTimeStep = (dt > outputTimeInterval && !physics.getParams().useAdaptiveTimeStep());
        }

        // Use default time step when injection starts
        if(params.doStimulation() && (time+dt > tInj_start && time < tInj_start))
        {
          debug_jochen << "== Injection starts, resetting dt to start value " << dtstart << std::endl;
          dt = dtstart;
        }

        bool acceptTimeStep = not physics.getParams().useAdaptiveTimeStep()
            || (time < tEquilibrium   // Fixed time step during equilibration
                  || std::abs(time - tEquilibrium) < 1e-6);
                  // Do not modify resetted time step when switching on active channels

        typename Traits::GF_MEMB_FLUX::Traits::RangeType oldMaxFlux = gfMembFlux.getMaxFlux();
        // Update channels
        gfMembFlux.updateChannels(membGV, time, dt);
        gfMembFlux.updateFlux(membGV, time, dt);
        while(not acceptTimeStep)
        {
          //Real dpot_dt = uPotChange.base().infinity_norm() / dt;
          //Real last_dpot_dt = last_dpot / last_dt;

          typename Traits::GF_MEMB_FLUX::Traits::RangeType maxFlux = gfMembFlux.getMaxFlux();

          debug_jochen << "==========================================" << std::endl;
          debug_jochen << "dt = " << dt << std::endl;
          typename Traits::GF_MEMB_FLUX::Traits::RangeFieldType oldMaxTotalFlux(0.0);
          typename Traits::GF_MEMB_FLUX::Traits::RangeFieldType maxTotalFlux(0.0);
          for(int j=0; j<NUMBER_OF_SPECIES; j++)
          {
            //debug_jochen << "Last time step / this time step (flux*dt) = "
            //    << (oldMaxFlux[j] * dt) << " / " << (maxFlux[j] * dt) << std::endl;
            oldMaxTotalFlux += oldMaxFlux[j];
            maxTotalFlux += maxFlux[j];
          }
          // Factor representing the ratio of the current flux with respect to the one from the last time step
          typename Traits::GF_MEMB_FLUX::Traits::RangeFieldType fluxFactor = maxFlux.one_norm() / oldMaxFlux.one_norm();

          //debug_jochen << std::endl << "Old / this (max flux * dt) = "
          //    << (oldMaxTotalFlux * dt) << " / " << (maxTotalFlux * dt)
          //    << " [factor " << fluxFactor << "]" << std::endl;

          //diagInfo.addDebugData(std::string("flux_factor"), fluxFactor);

          // (1) 'Soft' criterion: Try to adjust time step according to number Newton iterations
          if(diagInfo.iterations < 3 && diagInfo.iterations <= lastDiagInfo.iterations)
          {
            dt *= 1.1; // Carefully increase time step
          }
          if(diagInfo.iterations >= 5)
          {
            dt /= 1.2;  // Use smaller time step when we need too many Newton iterations
          }

          // (2) 'Harder' criterion: Bound time step according to change in potential
          //Real rel_dpot_dt = std::abs((dpot_dt-last_dpot_dt)/(std::max(dpot_dt,last_dpot_dt)));
          //diagInfo.addDebugData(std::string("rel_dpot_dt"), rel_dpot_dt);

//          if(time+dt > tInj_start && time+dt < tInj_end)
//          {
//            if(rel_dpot_dt > 0.05)
//            {
//              dt /= 2;
//            }
//            if(rel_dpot_dt < 0.001)
//            {
//              dt *= 1.2;
//            }
//          }


          // (3) 'Hardest' criterion: Bound maximum/minimum time step
          dt = std::min(dt, 0.05e-3 / physics.getTimeScale()); // maximum 0.05 ms
          dt = std::max(dt, 0.05e-6 / physics.getTimeScale()); // minimum 0.05 µs

          // This should be the _maximum_ membrane potential in order to work!
          // Update: This is true, previously the time step is only limited if membrane potential
          // was above -50mV _everywhere_!!
          typename Traits::DGF_POT::Traits::RangeType maxMembPot = -1e100;
          for(typename Traits::Physics::SubDomainElementIterator sdeit = membGV.template begin<0>();
              sdeit != membGV.template end<0>(); ++sdeit)
          {
            typename Traits::DGF_POT::Traits::RangeType membPot;
            physics.getMembranePotential(membGV.grid().multiDomainEntity(*sdeit), dgfPot, membPot);
            maxMembPot = std::max(maxMembPot, membPot);
          }
          // Get maximum membrane potential on all processors
          maxMembPot = gv.comm().max(maxMembPot);

          // Limit time step during AP
          if((params.doStimulation() && time+dt > tInj_start && time+dt < tInj_end)
              //|| (time+dt > 15e3+tInj_start && time+dt < 15e3+tInj_end)
              || physics.convertTo_mV(maxMembPot) > -50.) // Check if membrane potential is above -50mV
          {
            dt = std::min(dt, 0.01e-3 / physics.getTimeScale()); // maximum 0.01 ms during AP
          }

          debug_jochen << "Last time step #iterations = " << diagInfo.iterations
              << " => new time step dt = " << dt << std::endl;

          debug_jochen << "==========================================" << std::endl;

          acceptTimeStep = true;
        }
        debug_info << "Calculated timestep: " << dt << std::endl;

        // Update channels
        debug_jochen << "# gfMembFlux.updateChannels" << std::endl;
        gfMembFlux.updateChannels(membGV, time, dt);
        debug_jochen << "# gfMembFlux.updateFlux" << std::endl;
        gfMembFlux.updateFlux(membGV, time, dt);
        debug_jochen << "# gfMembFlux.acceptTimeStep" << std::endl;
        gfMembFlux.acceptTimeStep();

        // Update values representing last iteration
        //last_dpot = dpot;
        lastDiagInfo = diagInfo; // save diagnostic information for usage in next time step

        diagInfo.clear();

        bool converged = false;
        const int numberRestarts = 3;
        int countNumberRestarts = 0;

        // Restart mechanism
        while(not converged)
        {
          // Hack to let PDE solver be able to get the current time
          diagInfo.time = (time+dt);

          // !! Time-dependent boundary condition TYPES are currently not implemented in PDELab !!
          // evaluate constraints for current time step
          // bctypePot.setTime(time+dt);
          // bctypeCon.setTime(time+dt);
          // ccCon.clear();
          // Dune::PDELab::constraints( bctypeCon, gfsCon, ccCon );
          // ccPot.clear();
          // Dune::PDELab::constraints( bctypePot, gfsPot, ccPot );
          // TODO As far as I know, also time-dependent Dirichlet values are not supported by dune-multidomain!

          physics.setTimeStep(dt);
          diagInfo.dt = dt;

          try{

            // ========= DO ONE TIME STEP ==============
            //debug_verb << "================= Solving fully implicit system..." << std::endl;
            solver.apply(time,dt,uold,unew);
            //debug_verb << "================= Solved fully implicit system" << std::endl;

          } catch(Dune::PDELab::NewtonError& e)
          {
            countNumberRestarts++;

            if(countNumberRestarts <= params.maxNumberNewtonRestarts())
            {
              debug_warn << "====================================================" << std::endl;
              debug_warn << e << std::endl;
              debug_warn << "====================================================" << std::endl;

              debug_warn << "Newton did not converge for dt = " << dt;
              dt *= 0.5;
              debug_warn  << ", trying new time step dt = " << dt << "!" << std::endl;

              // Recalculate membrane flux
              gfMembFlux.updateChannels(membGV, time, dt);
              gfMembFlux.updateFlux(membGV, time, dt);
              gfMembFlux.acceptTimeStep();

              continue;
            }
            else
            {
              // Write out non-converged solution
              time += dt;
              acme2Output.writeStep(time);
              Tools::pecletNumber(gv, parametersCon, dgfGradPot);
              //Output::printMultipleComponentCoefficientVector(uCon, NUMBER_OF_SPECIES);
              throw e;
            }
          }

          // The following two lines are essential and must NOT be removed!
          // Grid functions (e.g. charge density) rely on the current value of child GFS coefficient vectors
          //Tools::compositeToChildCoefficientVector(multigfs, unew, unewCon, 0);
          //Tools::compositeToChildCoefficientVector(multigfs, unew, unewPot, 1);

          // Fill diagnosticInfo
          diagInfo.iterations = solver.getPDESolver().result().iterations;
          diagInfo.dt = dt;
          converged = true;
        }

        uChange = unew;
        uChange -= uold;

        // TODO Re-implement calculation of equation-specific changes!
        //uConChange = unewCon;
        //uConChange -= uoldCon;
        //uPotChange = unewPot;
        //uPotChange -= uoldPot;

//        dpot = uPotChange.base().infinity_norm();
//        uPotChange_Rel = unewPot;
//        int n = unewPot.flatsize();
//        for(int i=0; i<n; ++i)
//        {
//          Traits::VBE::access(uPotChange_Rel, i) /= Traits::VBE::access(uoldPot, i);
//          Traits::VBE::access(uPotChange_Rel, i) -= 1;
//
//          //debug_jochen << uoldPot[i] << " --> " << unewPot[i] << " (" << uPotChange_Rel[i] << std::endl;
//        }
//        uConChange_Rel = unewCon;
//        for(int i=0; i<uConChange_Rel.flatsize(); i++)
//        {
//          Traits::VBE::access(uConChange_Rel,i) /= Traits::VBE::access(uoldCon, i);
//          Traits::VBE::access(uConChange_Rel, i) -= 1;
//        }
        // Open new scope so that Dune::ios_base_all_saver's destructor can restore stream setting after output
        {
          Dune::ios_base_all_saver hackbraten(std::cout);
          debug_info << std::scientific << std::setprecision(16);
          debug_info << "L2 solution change: " << uChange.base().two_norm() << std::endl;
          debug_info << "MAX solution change: " << uChange.base().infinity_norm() << std::endl;
//          debug_info << "L2  con change  = " << uConChange.base().two_norm() << std::endl;
//          debug_info << "MAX con change  = " << uConChange.base().infinity_norm() << std::endl;
//          debug_info << "L2  pot change  = " << uPotChange.base().two_norm() << std::endl;
//          debug_info << "MAX pot change  = " << uPotChange.base().infinity_norm() << std::endl;
//
//          diagInfo.addDebugData(std::string("abs max pot change"), uPotChange.base().infinity_norm());
//          diagInfo.addDebugData(std::string("rel max pot change"), (uPotChange_Rel.base().infinity_norm()));
//          diagInfo.addDebugData(std::string("abs max con change"), uConChange.base().infinity_norm());
          diagInfo.setDebugData(std::string("old_dt"), old_dt);
          diagInfo.setDebugData(std::string("last_dt"), last_dt);
//          diagInfo.addDebugData(std::string("dpot"), dpot);
//          diagInfo.addDebugData(std::string("last_dpot"), last_dpot);
//          diagInfo.addDebugData(std::string("rel max con change"), (uConChange_Rel.base().infinity_norm()));
//          diagInfo.addDebugData(std::string("dpot_dt"), uPotChange.base().infinity_norm()/dt);
//          diagInfo.addDebugData(std::string("last_dpot_dt"), last_dpot/last_dt);
        }

        uold = unew;
        //uoldCon = unewCon;
        //uoldPot = unewPot;

        time += dt;
        iterations += diagInfo.iterations;
        timeSteps++;
        totalTimeSteps++;

        double diffToNextTimeStep = time - outputCounter * outputTimeInterval;
        if (printEveryTimeStep || std::abs(diffToNextTimeStep) < 1e-8 || diffToNextTimeStep > 0)
        {
          ++outputCounter;

          //Output::printSingleCoefficientVectorDG(uPot, "pot");
          //Output::printMultipleComponentCoefficientVectorDG(unewCon, NUMBER_OF_SPECIES);

          acme2Output.writeStep(time);
          Tools::pecletNumber(gv, parametersCon, dgfGradPot);

          typename Traits::DGF_POT::Traits::RangeType minMembPot = 1e100;
          typename Traits::DGF_POT::Traits::RangeType maxMembPot = -1e100;
          for(typename Traits::Physics::SubDomainElementIterator sdeit = membGV.template begin<0>();
              sdeit != membGV.template end<0>(); ++sdeit)
          {
            typename Traits::DGF_POT::Traits::RangeType membPot;
            physics.getMembranePotential(membGV.grid().multiDomainEntity(*sdeit), dgfPot, membPot);

            if(membPot < minMembPot) minMembPot = membPot;
            if(membPot > maxMembPot) maxMembPot = membPot;

            if(std::abs(time-physics.getParams().tEquilibrium()) < 1e-6)
            {
              typename GoldmanGF::Traits::RangeType goldmanPotential;
              // Evaluate at cell center
              typename GoldmanGF::Traits::DomainType x(0.5);

              goldmanGF.evaluate(*sdeit, x, goldmanPotential);

              debug_verb << sdeit->geometry().global(x) << std::endl;
              debug_verb << "  Membrane potential: " << physics.convertTo_mV(membPot) << std::endl;
              debug_verb << "  Goldman potential:  " << goldmanPotential << std::endl;
            }
          }
          debug_info << "MAX membrane potential = " << physics.convertTo_mV(maxMembPot) << " mV" << std::endl;
          debug_info << "MIN membrane potential = " << physics.convertTo_mV(minMembPot) << " mV" << std::endl;

          debug_info << std::endl << "########## output done, time: " << time
              << " [" << timeSteps << " steps, average #iterations: " << (iterations/timeSteps) << "] ###########"
              << std::endl << std::endl;

          iterations = 0;
          timeSteps = 0;
        }

        // Check if a checkpoint output will be created. A checkpoint is created for the very first
        // time steps and then every k time steps [with k == params.checkpointInterval()]
        if(params.doCheckpointing() &&
            (totalTimeSteps % params.checkpointInterval() == 0  || totalTimeSteps == 1))
        {
          acme2Output.saveState(time,dt);
        }


      }
      // ==============================================================================================

      debug_info << "Total number of time steps: " << totalTimeSteps << std::endl;
    }
};


#endif /* DUNE_AX1_ACME2_FULLY_COUPLED_HH */
