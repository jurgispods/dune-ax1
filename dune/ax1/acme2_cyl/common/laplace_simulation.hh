/*
 * acme2_cyl_fully_coupled.hh
 *
 *  Created on: Dec 15, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_LAPLACE_SIMULATION_HH
#define DUNE_AX1_LAPLACE_SIMULATION_HH

#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/ax1_boundaryfunction_membranefunction_adapter.hh>

template<typename Traits,typename ACME2CYL_OUTPUT>
class LaplaceSimulation
{
  public:

    typedef typename Traits::Real Real;

    template<typename SOLVER, typename ELEC_SUBPROBLEM>
    static void run(Real& time, Real& dt, Real& dtstart, Real& tend, const Real& tEquilibrium,
        typename Traits::Physics& physics,
        const typename Traits::GridView& gv, const typename Traits::SubGridView& membGV,
        SOLVER& solver,
        typename Traits::MultiGFS& multigfs,
        typename Traits::U& uold, typename Traits::U& unew,
        typename Traits::MultiGFS::template ConstraintsContainer<Real>::Type& cc,
        typename Traits::MultiGFS::template ConstraintsContainer<Real>::Type& ccWithoutOverlap,
        typename Traits::DGF_POT& dgfPot, typename Traits::DGF_POT_GRAD& dgfGradPot,
        ACME2CYL_OUTPUT& acme2_cylOutput,
        Acme2CylSolutionVectors<Traits>& solutionVectors,
        typename Traits::INITIAL_ELEC& initialElec, ELEC_SUBPROBLEM& elecSubProblem)
    {
      const Acme2CylParameters& params = physics.getParams();

      typename ACME2CYL_OUTPUT::DiagnosticInfo& diagInfo = acme2_cylOutput.getDiagInfo();
      typename ACME2CYL_OUTPUT::DiagnosticInfo lastDiagInfo(diagInfo); // copy of diagInfo

      Real time0 = time;
      if(dt <= 0.0)
        dt = 1.0;

      Real dpot(0.0);
      Real old_dt = dt;
      Real last_dt = dt;
      Real last_dpot(0.0);

      diagInfo.tEquilibrium = tEquilibrium;

      diagInfo.registerDebugData("first_defect",0.0);
      diagInfo.registerDebugData("defect",0.0);
      diagInfo.registerDebugData("reduction",0.0);
      diagInfo.registerDebugData("conv_rate",0.0);
      diagInfo.registerDebugData("t_assembler",0.0);
      diagInfo.registerDebugData("t_lin_solv",0.0);
      diagInfo.registerDebugData("t_elapsed",0.0);
      diagInfo.registerDebugData("avg_lin_it",0.0);

      diagInfo.registerDebugData(std::string("newton_restarts"), 0.0);
      diagInfo.registerDebugData(std::string("avg_lin_it"), 0.0,
          ACME2CYL_OUTPUT::DiagnosticInfo::ReductionFunction::Function_Mean);
      diagInfo.registerDebugData(std::string("old_dt"), old_dt);
      diagInfo.registerDebugData(std::string("last_dt"), last_dt);

      // Only initialize (i.e., delete) diagnostics file if simulation is not continued
      if(not (physics.getParams().doLoadState() && physics.getParams().doContinueSimulation()))
      {
        // This creates the file with all the header information from the
        // debug data dynamically added via 'registerDebugData()'
        acme2_cylOutput.initDiagInfoFile();
      }

      // Don't do initial output when a previous simulation is continued
      if(not (physics.getParams().doLoadState() && physics.getParams().doContinueSimulation()))
      {
        acme2_cylOutput.writeStep(time);

        Real minMembPot = std::numeric_limits<Real>::max();
        Real maxMembPot = std::numeric_limits<Real>::lowest();
        for(typename Traits::Physics::MIterator mit = physics.mInteriorBegin(); mit != physics.mInteriorEnd(); ++mit)
        {
          typename Traits::DGF_POT::Traits::RangeType membPot;
          physics.getMembranePotential(*mit, dgfPot, membPot);
          if(membPot < minMembPot) minMembPot = membPot;
          if(membPot > maxMembPot) maxMembPot = membPot;
        }
        minMembPot = gv.comm().min(minMembPot);
        maxMembPot = gv.comm().max(maxMembPot);
        debug_info << "MAX membrane potential = " << physics.convertTo_mV(maxMembPot) << " mV" << std::endl;
        debug_info << "MIN membrane potential = " << physics.convertTo_mV(minMembPot) << " mV" << std::endl;
        debug_info  << std::endl << "########## initial output done" << std::endl << std::endl;
      }
      // ==============================================================================================


      Real outputTimeInterval = physics.getParams().getOutputTimeInterval();
      int outputCounter = 1;
      bool printEveryTimeStep = (dt > outputTimeInterval && !physics.getParams().useAdaptiveTimeStep());

      double tInj_start = physics.getParams().tInj_start();
      double tInj_end   = physics.getParams().tInj_end();

      if(params.doStimulation())
      {
        debug_jochen << "tInj_start = " << tInj_start << ", tInj_end = " << tInj_end << std::endl;
      }

      // These variables are needed when reading in timesteps from an external file
      std::ifstream timestepFile(physics.getParams().getTimeStepFile());
      int timestepEvery = params.general.get("timeStepEvery",1);
      double next_time_last = 0.0;
      double next_time = 0.0;
      if(physics.getParams().getTimeStepFile() != "")
      {
        debug_info << "Using fixed time points from file '" << physics.getParams().getTimeStepFile()
            << "'!" << std::endl;

        assert(timestepFile.good());

        double next_time = -1;
        if(physics.getParams().doContinueSimulation() || physics.getParams().general.get("startTime", 0.0) > 0.0)
        {
          // Try to find time value in order to resume simulation
          std::string prefix("");

          if(not physics.getParams().doContinueSimulation())
          {
            time = physics.getParams().general.get("startTime", 0.0);
          }

          next_time = time;

          // Find value in 2nd column of file
          int nLine = Tools::findValue(timestepFile, next_time, 2, prefix, true);

          debug_info << "Found start time=" << next_time << " in line " << nLine << " in '"
              << physics.getParams().getTimeStepFile() << "', resuming from there..." << std::endl;

          if(nLine < 0)
            DUNE_THROW(Dune::Exception, "Could not find the requested time " << time
                << " to continue simulation from in timestep file '"
                << physics.getParams().getTimeStepFile() << "!");
        } else {

          // Simply read first timestep value in file; make sure it is equal to our start time!
          std::string line;
          std::getline(timestepFile, line);
          std::stringstream line_str(line);
          double dummy = -1;
          next_time = 0.0;
          // Extract first (unused) timestep
          line_str >> dummy >> next_time;
        }

        assert(std::abs(time-next_time) < 1e-4);
      }

      if(physics.getParams().doStimulation() && physics.getParams().doEquilibration())
      {
        assert(tInj_start > tEquilibrium);
      }

      // Check if time-dependent Dirichlet constraints shall be interpolated in each time step
      std::string loadBoundaryLocation = params.boundary.get("loadBoundary","bottom");
      // Loading time-dependent Dirichlet values only implemented for potential for now!
      bool isLoadedBoundaryDirichlet = params.isBoundaryDirichlet_Potential(loadBoundaryLocation);

      bool useTimeDependentBoundaryValues = params.boundary.get("useTimeDependentBoundaryValues",false);

      typename Traits::U uChange = unew;

      // Total number of time steps
      int totalTimeSteps = 0;

      // Number of iterations/time steps since last output
      int iterations = 0;
      int timeSteps = 0;

      // This flag only becomes true once when the equilibration is over, so that the time step can be resetted
      bool equilibrationOver = false;

      // ========= time loop ==========================================================================
      while (time<tend-1e-8)
      {
        debug_verb << "TIME = " << time << std::endl;
        debug_verb << "Last dt = " << dt << std::endl;

        // Check if this is the first time point after completing the equilibration
        if(physics.getParams().doEquilibration()
            && (time > tEquilibrium || std::abs(time-tEquilibrium) < 1e-6) //time >= tEquilibrium
            && not equilibrationOver)
        {
          equilibrationOver = true;
        }

        // ======================= Choose an appropriate time step dt ====================================
        last_dt = old_dt;
        old_dt = dt;

        bool acceptTimeStep = not physics.getParams().useAdaptiveTimeStep() // Fixed time step
            || (physics.getParams().doEquilibration() && time < tEquilibrium);   // Fixed time step during equilibration

        // Reset to default time step when equilibration was completed
        if(equilibrationOver)
        {
          debug_info << "== Equilibration over, resetting dt to start value " << dtstart << std::endl;
          dt = dtstart;
          acceptTimeStep = true;
          bool printEveryTimeStep = (dt > outputTimeInterval && !physics.getParams().useAdaptiveTimeStep());
        }

        // Use default time step when injection starts
        if(params.doStimulation() && (time+dt > tInj_start && time < tInj_start))
        {
          debug_info << "== Injection starts, resetting dt to start value " << dtstart << std::endl;
          acceptTimeStep = true;
          dt = dtstart;
        }

        // When loading and continuing an old simulation, take the last used time step as the first dt!
        if(params.doLoadState() && params.doContinueSimulation()
            && std::abs(time-time0) < 1e-6)
        {
          debug_info << "== Continuing old simulation with from loaded time " << time << " with dt=" << dt << std::endl;
          acceptTimeStep = true;
        }

        // Use fixed time step to match next time value from timestep file
        if(not physics.getParams().useAdaptiveTimeStep() && physics.getParams().getTimeStepFile() != "")
        {
          // Read next line in timestep file
          double dummy = -1;
          double next_time = 0.0;
          std::string line;

          // Read the timestepEvery-th next timestep (default: 1)
          for(int i=0; i<timestepEvery; i++)
          {
            bool hasLine = std::getline(timestepFile, line);
            // std::getline returned false => end of timestep file reached, we are done
            if(! hasLine) break;
          }

          std::stringstream line_str(line);
          line_str >> dummy >> next_time;

          // Set dt such that the next prescribed timestep is hit
          dt = next_time - time;
          debug_jochen << "Next time: " << next_time << " => dt = " << dt << std::endl;

          acceptTimeStep = true;
        }

        while(not acceptTimeStep)
        {
          //Real dpot_dt = uPotChange.base().infinity_norm() / dt;
          //Real last_dpot_dt = last_dpot / last_dt;

          debug_jochen << "==========================================" << std::endl;
          debug_jochen << "dt = " << dt << std::endl;

          // (1) 'Soft' criterion: Try to adjust time step according to number Newton iterations
          const int lowerLimit = params.getTimeStepLowerLimitNewtonIt(); // default 3
          const int upperLimit = params.getTimeStepUpperLimitNewtonIt(); // default 5
          if(diagInfo.iterations < lowerLimit && diagInfo.iterations <= lastDiagInfo.iterations)
          {
            dt *= 1.1; // Carefully increase time step
          }
          if(diagInfo.iterations >= upperLimit)
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
          dt = std::min(dt, params.getMaxTimeStep() / physics.getTimeScale()); // maximum timestep (default 0.05 ms)
          dt = std::max(dt, params.getMinTimeStep() / physics.getTimeScale()); // minimum timestep (default 0.05 µs)

          Real minMembPot = std::numeric_limits<Real>::max();
          Real maxMembPot = std::numeric_limits<Real>::lowest();
          for(typename Traits::Physics::MIterator mit = physics.mInteriorBegin(); mit != physics.mInteriorEnd(); ++mit)
          {
            typename Traits::DGF_POT::Traits::RangeType membPot;
            physics.getMembranePotential(*mit, dgfPot, membPot);
            if(membPot < minMembPot) minMembPot = membPot;
            if(membPot > maxMembPot) maxMembPot = membPot;
          }
          minMembPot = gv.comm().min(minMembPot);
          maxMembPot = gv.comm().max(maxMembPot);

          // (4) Additionally limit time step during AP
          if((params.doStimulation() && time+dt > tInj_start && time+dt < tInj_end)
              //|| (time+dt > 15e3+tInj_start && time+dt < 15e3+tInj_end)
              || physics.convertTo_mV(maxMembPot) > -50.) // Check if membrane potential is _anywhere_ above -50mV
          {
            dt = std::min(dt, params.getMaxTimeStepAP() / physics.getTimeScale()); // maximum timestep (default: 0.01ms) during AP
          }

          debug_jochen << "Last time step #iterations = " << diagInfo.iterations
              << " => new time step dt = " << dt << std::endl;

          debug_jochen << "==========================================" << std::endl;

          // Find best matching timestep from timestep file and adjust dt such that time value
          // is matched exactly!
          if(physics.getParams().getTimeStepFile() != "")
          {
            // Read next line in timestep file
            double dummy = -1;
            double desired_time = time+dt;
            debug_info << "Desired next time value: " << desired_time << std::endl;

            assert(timestepFile.good());
            while(next_time < desired_time && timestepFile.good())
            {
              // Save last found value
              next_time_last = next_time;

              std::string line;
              std::getline(timestepFile, line);
              std::stringstream line_str(line);
              line_str >> dummy >> next_time;

              debug_jochen << "Found time value " << next_time << std::endl;
            }
            // Found value next_time should now be the first occurring value >= desired_time,
            // so next_time_last should be the last value < desired_time
            if(next_time < desired_time || next_time_last >= desired_time)
            {
              // Handle case when desired_time is > tend (last time step)
              if(desired_time >= tend)
              {
                debug_jochen << "Last time step detected!" << std::endl;
              } else {
                DUNE_THROW(Dune::Exception, "Something went wrong when trying to find a matching time step from '"
                  << physics.getParams().getTimeStepFile() << "', desired time: " << desired_time
                  << ", best match: " << next_time_last);
              }
            }

            // Choose new time value to be the one with the smallest difference to the desired time value;
            // set dt such that the next prescribed timestep is hit
            if(std::abs(next_time_last-desired_time) < std::abs(next_time-desired_time))
            {
              dt = next_time_last - time;
            } else {
              dt = next_time - time;
            }
            debug_info << "Next time: " << (time+dt) << " => dt = " << dt << std::endl;
          }

          acceptTimeStep = true;
        }
        debug_info << "Calculated timestep: " << dt << std::endl;

        Real dt_max = gv.comm().max(dt);
        Real dt_min = gv.comm().min(dt);
        if(std::abs(dt - dt_max)>1e-6 || std::abs(dt - dt_min)>1e-6)
        {
          DUNE_THROW(Dune::Exception, "Calculated time steps on processors don't match! min_dt = "
              << dt_min << ", dt_max = " << dt_max << "!");
        }

        // ==================================================================================================

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

          physics.setTimeStep(dt);
          diagInfo.dt = dt;

          // Set time also in initialGF (and implicitly also in parameter class!)
          initialElec.setTime(time+dt);

          try{
            // ========= DO ONE TIME STEP ==============
            //debug_verb << "================= Solving fully implicit system..." << std::endl;

            // Time-dependent Dirichlet values:
            // Only do interpolation if the boundary values to be loaded are Dirichlet values (i.e. constraints!),
            // otherwise, the parameter class will provide the time-dependent Neumann values directly to the operator!
            if(useTimeDependentBoundaryValues && isLoadedBoundaryDirichlet)
            {
              // Save unew: It might have been modified in a failed Newton iteration, so unew != uold might
              // hold at this point. Use the modified unew rather than uold, as it might already contain
              // better starting values!
              typename Traits::U unew_save(unew);
              Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,unew,initialElec,elecSubProblem);

              //Dune::PDELab::copy_nonconstrained_dofs(cc, unew_save, unew);
              // IMPORTANT: Use the modified contraints container ccWithoutOverlap to make this work in parallel!
              // This will also copy overlap DOFs, even though they are considered as constrained in the original cc
              Dune::PDELab::copy_nonconstrained_dofs(ccWithoutOverlap, unew_save, unew);
            }
            // Reset unew in this case of a stationary problem; previous solution might not be a good
            // starting value, especially when the time-dependent boundary condition source is taken
            // every k timesteps!
            if(! params.general.get("reuseSolution",true))
            {
              unew = 0.0;
              // Make sure to copy Dirichlet constraints nevertheless!!
              Dune::PDELab::copy_constrained_dofs(ccWithoutOverlap, uold, unew);
            }
            solver.apply(unew);
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

              continue;
            }
            else
            {
              debug_warn << "====================================================" << std::endl;
              debug_warn << e << std::endl;
              debug_warn << "====================================================" << std::endl;

              debug_warn << "Newton did not converge for dt = " << dt
                  <<  " maximum number of Newton restarts (" << params.maxNumberNewtonRestarts()
                  << ") exceeded!" << std::endl;

              debug_warn << "  first defect: " << solver.result().first_defect << std::endl;
              debug_warn << "        defect: " << solver.result().defect << std::endl;

              // Write out non-converged solution
              time += dt;
              acme2_cylOutput.writeStep(time);
              throw e;
            }
          }

          double avgLinearIterations = solver.result().linear_solver_iterations;
          double avgLinearSolverTime = solver.result().linear_solver_time;
          if(solver.result().iterations > 0)
          {
            avgLinearIterations /= solver.result().iterations;
            avgLinearSolverTime /= solver.result().iterations;
          }
          if(! solver.result().converged)
          {
            debug_warn << "===================================================================" << std::endl;
            debug_warn << "Solver did not converge! Achieved reduction: " << solver.result().reduction << std::endl;
            debug_warn << "===================================================================" << std::endl;
            if(params.general.get("linearProblemEnforceConvergence",true))
              DUNE_THROW(Dune::Exception, "Stationary linear problem solver did not converge!");
          } else {
            debug_info << "  Achieved reduction: " << solver.result().reduction << std::endl;
          }

          {
            Dune::ios_base_all_saver wurst(std::cout);
            debug_info << "  Total number of linear iterations: " << solver.result().linear_solver_iterations
                      << std::fixed << " (average " << avgLinearIterations << ")"
                      << std::endl;
            debug_info << "  Total linear solver time: " << solver.result().linear_solver_time
                      << std::fixed << " (average " << avgLinearSolverTime << ")"
                      << std::endl;
          }
          //Output::printRawCoefficientVector(unew, "unew");


          // Fill diagnosticInfo
          diagInfo.setDebugData("first_defect",solver.result().first_defect);
          diagInfo.setDebugData("defect",solver.result().defect);
          diagInfo.setDebugData("reduction",solver.result().reduction);
          diagInfo.setDebugData("conv_rate",solver.result().conv_rate);
          diagInfo.setDebugData("t_assembler",solver.result().assembler_time);
          diagInfo.setDebugData("t_lin_solv",solver.result().linear_solver_time);
          diagInfo.setDebugData("t_elapsed",solver.result().elapsed);
          diagInfo.setDebugData("avg_lin_it",avgLinearIterations);

          diagInfo.iterations = solver.result().iterations;
          diagInfo.dt = dt;
          diagInfo.setDebugData(std::string("newton_restarts"), countNumberRestarts);
          converged = true;
        }

        uChange = unew;
        uChange -= uold;

        // Open new scope so that Dune::ios_base_all_saver's destructor can restore stream setting after output
        {
          Dune::ios_base_all_saver hackbraten(std::cout);
          debug_info << std::scientific << std::setprecision(16);
          debug_info << "L2 solution change: " << uChange.base().two_norm() << std::endl;
          debug_info << "MAX solution change: " << uChange.base().infinity_norm() << std::endl;
          diagInfo.setDebugData(std::string("old_dt"), old_dt);
          diagInfo.setDebugData(std::string("last_dt"), last_dt);
        }

        uold = unew;

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

          acme2_cylOutput.writeStep(time);

          Real minMembPot = std::numeric_limits<Real>::max();
          Real maxMembPot = std::numeric_limits<Real>::lowest();
          for(typename Traits::Physics::MIterator mit = physics.mInteriorBegin(); mit != physics.mInteriorEnd(); ++mit)
          {
            typename Traits::DGF_POT::Traits::RangeType membPot;
            physics.getMembranePotential(*mit, dgfPot, membPot);

            if(membPot < minMembPot) minMembPot = membPot;
            if(membPot > maxMembPot) maxMembPot = membPot;

            if(std::abs(time-physics.getParams().tEquilibrium()) < 1e-6)
            {
              Dune::FieldVector<double,1> x(0.5);
              debug_verb << mit->geometry().global(x) << std::endl;
              debug_verb << "  Membrane potential: " << physics.convertTo_mV(membPot) << std::endl;
            }
            typename Traits::Physics::ElementPointer outside = mit->outside();
            typename Traits::DGF_POT::Traits::DomainType x(mit->geometryInOutside().center());
            typename Traits::DGF_POT::Traits::RangeType extPot;
            dgfPot.evaluate(*outside, x, extPot);
            debug_verb << "  Extracellular membrane interface potential @" << outside->geometry().global(x)
              << ": " << physics.convertTo_mV(extPot) << std::endl;
          }
          minMembPot = gv.comm().min(minMembPot);
          maxMembPot = gv.comm().max(maxMembPot);

          debug_info << "MAX membrane potential = " << physics.convertTo_mV(maxMembPot) << " mV" << std::endl;
          debug_info << "MIN membrane potential = " << physics.convertTo_mV(minMembPot) << " mV" << std::endl;

          debug_info << std::endl << "########## timestep #" << acme2_cylOutput.getNTimeSteps()
              << ": output done, time: " << time
              << " [" << timeSteps << " steps, average #iterations: " << (iterations/timeSteps) << "] ###########"
              << std::endl << std::endl;

//          debug_jochen << "-------------------------------------------------------------" << std::endl;
//          debug_jochen << "CHANNEL STATES: " << std::endl;
//
//          Dune::ios_base_all_saver hack(std::cout);
//          std::cout << std::setprecision(32);
//          Output::printVector(physics.getMembrane().getChannelSet().serializeChannelStates());
//          debug_jochen << "-------------------------------------------------------------" << std::endl;

          iterations = 0;
          timeSteps = 0;
        }

        // Check if a checkpoint output will be created. A checkpoint is created for the very first
        // time steps and then every k time steps [with k == params.checkpointInterval()]
        if(params.doCheckpointing() &&
            (totalTimeSteps % params.checkpointInterval() == 0  || totalTimeSteps == 1))
        {
          std::string checkpointFilename("");
          acme2_cylOutput.saveState(time,dt,checkpointFilename);
        }

      }
      // ==============================================================================================

      debug_info << "Total number of time steps: " << totalTimeSteps << std::endl;
    }


};


#endif /* DUNE_AX1_ACME2CYL_FULLY_COUPLED_HH */
