/*
 * acme2_cyl_fully_coupled.hh
 *
 *  Created on: Dec 15, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_MORI_SIMULATION_HH
#define DUNE_AX1_ACME2CYL_MORI_SIMULATION_HH

#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/ax1_boundaryfunction_membranefunction_adapter.hh>

template<typename Traits,typename ACME2CYL_OUTPUT>
class Acme2CylMoriSimulation
{
  public:

    typedef typename Traits::Real Real;

    template<typename SOLVER_CON, typename SOLVER_POT, typename SolutionContainer,
      typename ELEC_SUBPROBLEM_CON, typename ELEC_SUBPROBLEM_POT>
    static void run(Real& time, Real& dt, Real& dtstart, Real& tend, const Real& tEquilibrium,
        typename Traits::Physics& physics,
        const typename Traits::GridView& gv, const typename Traits::SubGridView& membGV,
        SOLVER_CON& solverCon, SOLVER_POT& solverPot,
        typename Traits::PARAMETERS_POT& parametersPot,
        typename Traits::PARAMETERS_CON& parametersCon,
        typename Traits::MultiGFS& multigfs,
        typename Traits::U& uold, typename Traits::U& unew,
        typename Traits::MultiGFS::template ConstraintsContainer<Real>::Type& cc,
        typename Traits::MultiGFS::template ConstraintsContainer<Real>::Type& ccWithoutOverlap,
        typename Traits::GF_MEMB_FLUX& gfMembFlux,
        typename Traits::GF_MORI_FLUX& gfMoriFlux,
        typename Traits::DGF_CON& dgfCon,
        typename Traits::DGF_POT& dgfPot, typename Traits::DGF_POT_GRAD& dgfGradPot,
        ACME2CYL_OUTPUT& acme2_cylOutput,
        Acme2CylSolutionVectors<Traits>& solutionVectors,
        SolutionContainer& solutionContainer,
        typename Traits::INITIAL_ELEC& initialElec,
        ELEC_SUBPROBLEM_CON& elecSubProblemCon, ELEC_SUBPROBLEM_POT& elecSubProblemPot)
    {
      const Acme2CylParameters& params = physics.getParams();
      const Real timeScale = physics.getTimeScale();
      const Real lengthScale = physics.getLengthScale();

      // Calculate total electrolyte domain volume (used for termination criterion in iteration)
      CylinderGridVector<Real> yCyl(params.xMax()-params.xMin());
      Real totalVolume = yCyl.volume_r0_rend(params.yMin(), params.yMax());
      Real membVolume = 0.0;
      const std::vector<Real>& yMemb = params.yMemb();
      Real dMemb = params.dMemb();
      for(int i=0; i<yMemb.size(); ++i)
      {
        membVolume += yCyl.volume_r0_h(yMemb[i], dMemb);
      }
      Real elecVolume = totalVolume - membVolume;

      // tolerance = eps_tol * e * c_0 (with eps_tol = 1e-5, c_0 = 100 mmol/l = mol/m^3 = mM)
      Real tolerance = 1e-5 * con_e * 100;

      const bool fullyImplicit = params.boundary.get("fullyImplicitMembraneFlux",false);
      const bool useElectroneutralityCriterion = params.general.get("useElectroneutralityTerminationCriterion",true);
      const bool doConcentrationPostprocessing = params.general.get("doConcentrationPostprocessing",false);

      const int intorderPot = acme2_cylOutput.getIntOrderPot();
      const int intorderCon = acme2_cylOutput.getIntOrderCon();
      const int intorderadd = 2;

      const int maxIT = params.general.get("newtonMaxIt",50);

      // Create vector of DGFs, one for each concentration
      std::vector<typename ACME2CYL_OUTPUT::Traits::DGF_SINGLE_CON_MD> dgfSingleConVec;
      typedef MoriChargeGridFunction<typename Traits::GF_MORI_FLUX, typename Traits::DGF_POT, typename Traits::Physics>
        GF_MORI_CHARGE;
      GF_MORI_CHARGE gfMoriCharge(gfMoriFlux, dgfPot, physics);

      std::vector<typename ACME2CYL_OUTPUT::Traits::DGF_SINGLE_CON_MD::Traits::RangeType> initConcIntegralBulk;
      std::vector<typename ACME2CYL_OUTPUT::Traits::DGF_SINGLE_CON_MD::Traits::RangeType> initConcIntegralChargeLayer;
      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        typename ACME2CYL_OUTPUT::Traits::DGF_SINGLE_CON_MD dgfSingleCon(dgfCon, j);
        dgfSingleConVec.push_back(dgfSingleCon);
      }
      std::vector<typename ACME2CYL_OUTPUT::Traits::DGF_SINGLE_CON_MD::Traits::RangeType>
        concIntegralBulk(NUMBER_OF_SPECIES, -1);
      std::vector<typename ACME2CYL_OUTPUT::Traits::DGF_SINGLE_CON_MD::Traits::RangeType>
        concIntegralChargeLayer(NUMBER_OF_SPECIES, -1);

      const std::vector<int> subDomains = {CYTOSOL, ES};

      // scalingFactorBulk misses the valence 'z_i', which will be added in-place below
      const Real scalingFactorBulk = con_e * con_mol * lengthScale * lengthScale * lengthScale;
      const Real scalingFactorChargeLayer = lengthScale * lengthScale;

      typename ACME2CYL_OUTPUT::DiagnosticInfo& diagInfo = acme2_cylOutput.getDiagInfo();
      typename ACME2CYL_OUTPUT::DiagnosticInfo lastDiagInfo(diagInfo); // copy of diagInfo

      Real time0 = time;

      Real dpot(0.0);
      Real old_dt = dt;
      Real last_dt = dt;
      Real last_dpot(0.0);

      diagInfo.tEquilibrium = tEquilibrium;

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

      typedef Ax1GoldmanEquationGridFunction<typename Traits::DGF_CON,
          typename Traits::Physics> GoldmanGF;
      GoldmanGF goldmanGF(dgfCon,physics);

      // Initialize membrane flux classes once and for all
      if(! gfMembFlux.isCompletelyInitialized())
      {
        debug_jochen << "[init] # gfMembFlux.updateState" << std::endl;
        gfMembFlux.updateState(time, dt);
        debug_jochen << "# gfMembFlux.acceptTimeStep" << std::endl;
        gfMembFlux.acceptTimeStep();
      }
      if(! gfMoriFlux.isCompletelyInitialized())
      {
        debug_jochen << "[init] # gfMoriFlux.updateState" << std::endl;
        gfMoriFlux.updateState(time, dt);
        debug_jochen << "# gfMoriFlux.acceptTimeStep" << std::endl;
        gfMoriFlux.acceptTimeStep();
      }

      // Don't do initial output when a previous simulation is continued
      if(not (physics.getParams().doLoadState() && physics.getParams().doContinueSimulation()))
      {
        acme2_cylOutput.writeStep(time);
        Tools::pecletNumber(gv, parametersCon, dgfGradPot);

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
      bool isLoadedBoundaryDirichletPot = params.isBoundaryDirichlet_Potential(loadBoundaryLocation);

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
          std::getline(timestepFile, line);
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
        debug_info << "Calculated timestep: " << dt << " -> new time: " << (time+dt) << std::endl;

        Real dt_max = gv.comm().max(dt);
        Real dt_min = gv.comm().min(dt);
        if(std::abs(dt - dt_max)>1e-6 || std::abs(dt - dt_min)>1e-6)
        {
          DUNE_THROW(Dune::Exception, "Calculated time steps on processors don't match! min_dt = "
              << dt_min << ", dt_max = " << dt_max << "!");
        }

        // ==================================================================================================

        if(time > time0)
        {
          // (0) GATING STEP (outside of Mori operator-split iteration)
          // We use the uold=unew values for the update!
          debug_jochen << "# gfMembFlux.updateState" << std::endl;
          gfMembFlux.updateState(time, dt);
        }

        // Update boundary values
        debug_info << "# Calling prepareNextTimeStep(" << (time+dt) << ") on parameter classes" << std::endl;
        parametersPot.prepareNextTimeStep(time+dt);
        for(int j=0; j<parametersCon.size(); ++j)
        {
          parametersCon[j]->prepareNextTimeStep(time+dt);
        }

        // Calculate initial concentration integrals; we have to do this here in order to make sure that
        // both membrane and Mori flux GFs have been initialized to steady-state with respect to initial conditions
        // TODO: Add possibility to load initial concentration integrals from config file; this would be
        // mandatory when continuing a previous simulation!
        if(initConcIntegralBulk.size() == 0)
        {
          // Calculate initial concentration integral (charge layer)
          typename GF_MORI_CHARGE::Traits::RangeType concIntegralChargeLayer_temp;
          Acme2CylGeometryTools::integrateIntersectionGridFunctionOverCylinderSubdomain(gfMoriCharge, physics,
              subDomains, concIntegralChargeLayer_temp, intorderCon + intorderadd, false);

          debug_jochen << "^^^ Initial concentration integrals (bulk / charge layer)" << std::endl;
          Real sumBulk = 0;
          Real sumChargeLayer = 0;
          for(int j=0; j<NUMBER_OF_SPECIES; ++j)
          {
            // Calculate initial concentration integral (bulk)
            typename ACME2CYL_OUTPUT::Traits::DGF_SINGLE_CON_MD::Traits::RangeType concIntegralBulk_temp;
            Acme2CylGeometryTools::integrateGridFunctionOverCylinderSubdomain(dgfSingleConVec[j], physics,
                subDomains, concIntegralBulk_temp, intorderCon + intorderadd);

            // Bring to units of [C]
            concIntegralBulk_temp *= scalingFactorBulk * physics.getValence(j);
            initConcIntegralBulk.push_back(concIntegralBulk_temp);
            sumBulk += concIntegralBulk_temp;

            // Bring to units of [C]
            concIntegralChargeLayer_temp[j] *= scalingFactorChargeLayer;
            initConcIntegralChargeLayer.push_back(concIntegralChargeLayer_temp[j]);
            sumChargeLayer += concIntegralChargeLayer_temp[j];

            debug_jochen << "^^^ [" << ION_NAMES[j] << "] " << initConcIntegralBulk[j] << " / "
                << initConcIntegralChargeLayer[j] << std::endl;
          }
          sumBulk = gv.comm().sum(sumBulk);
          sumChargeLayer = gv.comm().sum(sumChargeLayer);
          debug_jochen << "^^^ [TOTAL] " << sumBulk << " / " << sumChargeLayer << std::endl;
        }


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

          int m = 0;
          try
          {
            typename Traits::U ulast = uold;

            // FIXME Implemement correct termination criterion
            bool iterationConverged = false;
            for(m = 0; m<maxIT && !iterationConverged; ++m)
            {
              debug_info << "==================== Mori operator-split, iteration #" << m << " ===================="
                  << std::endl;

              // Time-dependent Dirichlet values (supported for potential only);
              if(params.boundary.get("useTimeDependentBoundaryValuesPot",false) && isLoadedBoundaryDirichletPot)
              {
                debug_info << "== Applying time-dependent Dirichlet values..." << std::endl;

                // Save unew: It might have been modified in a failed Newton iteration, so unew != uold might
                // hold at this point. Use the modified unew rather than uold, as it might already contain
                // better starting values!
                typename Traits::U unew_save(unew);

                //std::vector<typename Traits::U> uu;
                //std::vector<std::string> uu_names(3);
                //uu.push_back(unew_save);
                //uu_names[0] = "unew_save";

                Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,unew,initialElec.template child<0>(),
                    elecSubProblemCon);
                Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,unew,initialElec.template child<1>(),
                    elecSubProblemPot);
                //uu.push_back(unew);
                //uu_names[1] = "unew after interpolation";

                //Dune::PDELab::copy_nonconstrained_dofs(cc, unew_save, unew);
                // IMPORTANT: Use the modified contraints container ccWithoutOverlap to make this work in parallel!
                // This will also copy overlap DOFs, even though they are considered as constrained in the original cc
                Dune::PDELab::copy_nonconstrained_dofs(ccWithoutOverlap, unew_save, unew);
                //uu.push_back(unew);
                //uu_names[2] = "unew after copy_nonconstrained_dofs";
                //Output::printRawCoefficientVector(uu, uu_names);

                //Output::printRawCoefficientVector(unew, "unew");
              }
              ulast = unew;

              // Populate solution container with solution of previous iteration; changes in unew will _not_
              // be reflected in ulast in order to implement Mori's operator-split approach

              // (1) POTENTIAL UPDATE
              // Remember: Since the potential Neumann boundary condition does not depend on the gamma
              // states of the Mori charge layer (their sum is always 1), it doesn't matter that those
              // states have not yet been updated at this point! All that matter is the old and new potential.
              // The old potential is known via 'uold' in the charge layer class, the new is included implicitly!
              solutionContainer.setSolutionConNew(&ulast);
              solutionContainer.setSolutionPotNew(&unew);
              // Update potential
              debug_info << "===== Solving potential subsystem... =========" << std::endl;
              solverPot.apply(unew);
              debug_info << "===== DONE. ==================================" << std::endl;

              // (2) MORI CHARGE LAYER GAMMA UPDATE
              // Gamma variables are updated using _new_ (== last, doesn't matter here) concentration
              // and, of course, _new_ potential values
              solutionContainer.setSolutionConNew(&ulast);
              solutionContainer.setSolutionPotNew(&unew);

              // Update Mori (gamma) states; this is now safe even at the initial time step
              debug_jochen << "# gfMoriFlux.updateState" << std::endl;
              gfMoriFlux.updateState(time, dt);

              // (3) CONCENTRATION UPDATE
              // Concentration operator needs _last_ concentration and potential values to evaluate the drift flux!
              solutionContainer.setSolutionConNew(&ulast);
              solutionContainer.setSolutionPotNew(&ulast);
              // Solve instationary linear concentration system
              debug_info << "===== Solving concentration subsystem... =====" << std::endl;
              solverCon.apply(time,dt,uold,unew);
              debug_info << "===== DONE. ==================================" << std::endl;

              // Make sure we use _new_ values once the iteration is terminated
              solutionContainer.setSolutionConNew(&unew);
              solutionContainer.setSolutionPotNew(&unew);

              debug_info << "^^^ Concentration integrals (bulk)" << std::endl;
              Real sumBulk = 0;
              for(int j=0; j<NUMBER_OF_SPECIES; ++j)
              {
                // Calculate concentration integral (bulk)
                Acme2CylGeometryTools::integrateGridFunctionOverCylinderSubdomain(dgfSingleConVec[j], physics,
                    subDomains, concIntegralBulk[j], intorderCon + intorderadd);
                // Bring to units of [C]
                concIntegralBulk[j] *= scalingFactorBulk * physics.getValence(j);
                sumBulk += concIntegralBulk[j];

                debug_info << "^^^ [" << ION_NAMES[j] << "] initial / current " << initConcIntegralBulk[j] << " / "
                    << concIntegralBulk[j] << std::endl;
              }
              sumBulk = gv.comm().sum(sumBulk);
              debug_jochen << "^^^ [TOTAL] " << sumBulk << std::endl;

              if(useElectroneutralityCriterion)
              {
                iterationConverged = std::abs(sumBulk / elecVolume) < tolerance;
                debug_info << "^^^ Mori termination criterion: " << std::abs(sumBulk / elecVolume) << " < "
                  << tolerance << " ? --> " << iterationConverged << std::endl;
              }
            }
            // Update solution containers with new values after iteration converged
            solutionContainer.setSolutionConNew(&unew);
            solutionContainer.setSolutionPotNew(&unew);

            debug_info << "==================== Mori operator-split converged after " << m
                << " iterations ====================" << std::endl;

          } catch(Dune::PDELab::NewtonError& e)
          {
            countNumberRestarts++;

            if(countNumberRestarts <= params.maxNumberNewtonRestarts())
            {
              debug_warn << "====================================================" << std::endl;
              debug_warn << e << std::endl;
              debug_warn << "====================================================" << std::endl;

              debug_warn << "Iteration did not converge for dt = " << dt;
              dt *= 0.5;
              debug_warn  << ", trying new time step dt = " << dt << "!" << std::endl;

              // Recalculate membrane flux
              debug_jochen << "# gfMembFlux.discardTimeStep" << std::endl;
              gfMembFlux.discardTimeStep(); // this doesn't actually do anything for now
              debug_jochen << "# gfMembFlux.updateState" << std::endl;
              gfMembFlux.updateState(time, dt);

              // Mori flux will be recalculated within the iteration loop
              debug_jochen << "# gfMoriFlux.discardTimeStep" << std::endl;
              gfMoriFlux.discardTimeStep(); // this doesn't actually do anything for now

              // Update boundary values
              debug_info << "# Calling prepareNextTimeStep(" << (time+dt) << ") on parameter classes" << std::endl;
              parametersPot.prepareNextTimeStep(time+dt);
              for(int j=0; j<parametersCon.size(); ++j)
              {
                parametersCon[j]->prepareNextTimeStep(time+dt);
              }

              continue;
            }
            else
            {
              debug_warn << "====================================================" << std::endl;
              debug_warn << e << std::endl;
              debug_warn << "====================================================" << std::endl;

              debug_warn << "Iteration did not converge for dt = " << dt
                  <<  " maximum number of Newton restarts (" << params.maxNumberNewtonRestarts()
                  << ") exceeded!" << std::endl;

              debug_warn << "  [con] first defect: " << solverCon.getPDESolver().result().first_defect << std::endl;
              debug_warn << "  [con]       defect: " << solverCon.getPDESolver().result().defect << std::endl;
              debug_warn << "  [pot] first defect: " << solverPot.result().first_defect << std::endl;
              debug_warn << "  [pot]       defect: " << solverPot.result().defect << std::endl;


              // Write out non-converged solution
              time += dt;
              acme2_cylOutput.writeStep(time);
              Tools::pecletNumber(gv, parametersCon, dgfGradPot);
              throw e;
            }
          } catch(Dune::Exception& e)
          {
            debug_warn << "====================================================" << std::endl;
            debug_warn << "DUNE Exception caught: " << e.what() << std::endl;
            debug_warn << "====================================================" << std::endl;

            // Write out non-converged solution
            time += dt;
            acme2_cylOutput.writeStep(time);
            Tools::pecletNumber(gv, parametersCon, dgfGradPot);
            throw e;
          } catch(std::exception& e)
          {
            debug_warn << "====================================================" << std::endl;
            debug_warn << "Exception caught: " << e.what() << std::endl;
            debug_warn << "====================================================" << std::endl;

            // Write out non-converged solution
            time += dt;
            acme2_cylOutput.writeStep(time);
            Tools::pecletNumber(gv, parametersCon, dgfGradPot);
            throw e;
          } catch(...)
          {
            debug_warn << "====================================================" << std::endl;
            debug_warn << "Unknown exception!" << std::endl;
            debug_warn << "====================================================" << std::endl;

            // Write out non-converged solution
            time += dt;
            acme2_cylOutput.writeStep(time);
            Tools::pecletNumber(gv, parametersCon, dgfGradPot);
            throw;
          }

          // Fill diagnosticInfo
          diagInfo.iterations = m;
          diagInfo.dt = dt;
          diagInfo.setDebugData(std::string("newton_restarts"), countNumberRestarts);
          converged = true;
        }

        uChange = 0.0;
        uChange -= uold;
        if(doConcentrationPostprocessing)
        {
          // Save non-postprocessed solution vector
          typename Traits::U unew_orig = unew;

          // Initial state will be wrong when resuming simulation
          if(params.doLoadState() && params.doContinueSimulation())
            DUNE_THROW(Dune::NotImplemented, "Resuming simulation with concentration postprocessing not implemented!");

          debug_info << "=================================================================================" << std::endl;
          debug_info << "   Postprocessing ion concentrations to avoid global charge accumulation..." << std::endl;
          debug_info << "=================================================================================" << std::endl;

          std::vector<Real> corrections(3, 1.0);
          debug_info << "^^^ Concentration integrals (bulk / charge layer)" << std::endl;

          // Calculate concentration integral (charge layer)
          typename GF_MORI_CHARGE::Traits::RangeType concIntegralChargeLayer_temp;
          Acme2CylGeometryTools::integrateIntersectionGridFunctionOverCylinderSubdomain(gfMoriCharge, physics,
              subDomains, concIntegralChargeLayer_temp, intorderCon + intorderadd, false);

          Real sumBulk = 0;
          Real sumChargeLayer = 0;
          for(int j=0; j<NUMBER_OF_SPECIES; ++j)
          {
            // Bring to units of [C]
            concIntegralChargeLayer_temp[j] *= scalingFactorChargeLayer;
            concIntegralChargeLayer[j] = concIntegralChargeLayer_temp[j];
            sumChargeLayer += concIntegralChargeLayer[j];

            // Ion concentration correction to avoid global charge accumulation when using Neumann-0 boundary condition
            // (MoriPeskin 2008)
            corrections[j] = (initConcIntegralBulk[j] + initConcIntegralChargeLayer[j] - concIntegralChargeLayer[j])
                / concIntegralBulk[j];

            debug_info << "^^^ [" << ION_NAMES[j] << "] initial: " << initConcIntegralBulk[j] << " / "
                << initConcIntegralChargeLayer[j] << std::endl;
            debug_info << "^^^ [" << ION_NAMES[j] << "] current: " << concIntegralBulk[j] << " / "
                << concIntegralChargeLayer[j] << std::endl;
            debug_info << "^^^  => corrections[" << ION_NAMES[j] << "] = " << corrections[j] << std::endl;
          }
          sumChargeLayer = gv.comm().sum(sumChargeLayer);
          debug_jochen << "^^^ [TOTAL] " << sumBulk << " / " << sumChargeLayer << std::endl;

          // Here comes the ugly hack: As we don't have a global mapping to extract a single GFS component from
          // a composite GFS coefficient vector, we need to assume we know the ordering of DOFs in the container
          // and do it manually!
          assert(params.nMembraneElements() <= 1);

          const int blockSize = AX1_BLOCKSIZE;
          assert(unew.flatsize() % blockSize == 0);

          //std::size_t size = unew.flatsize() / blockSize;
          std::size_t i = 0;
          for (typename Traits::U::iterator it = unew.begin(); it != unew.end(); ++it)
          {
            // Component is out of the interval [0 NUMBER_OF_SPECIES], the last value represents the potential
            // (which we don't need here)
            std::size_t component = i % blockSize;

            // Apply correction for the respective ion species
            if(component < NUMBER_OF_SPECIES)
            {
              *it *= corrections[component];
            }
            ++i;
          }
          debug_info << "=================================================================================" << std::endl;

          // Debug output (remove this later)
          std::vector<typename Traits::U> uu;
          uu.push_back(unew_orig);
          uu.push_back(unew);
          std::vector<std::string> uu_str;
          uu_str.push_back("unew before postprocessing");
          uu_str.push_back("unew after postprocessing ");
          Output::printRawCoefficientVector(uu, uu_str);
        }
        uChange += unew;

        // Open new scope so that Dune::ios_base_all_saver's destructor can restore stream setting after output
        {
          Dune::ios_base_all_saver hackbraten(std::cout);
          debug_info << std::scientific << std::setprecision(16);
          debug_info << "[pot] L2 solution change: " << uChange.base()[1].two_norm() << std::endl;
          debug_info << "[con] L2 solution change: " << uChange.base()[0].two_norm() << std::endl;
          debug_info << "[pot] MAX solution change: " << uChange.base()[1].infinity_norm() << std::endl;
          debug_info << "[con] MAX solution change: " << uChange.base()[0].infinity_norm() << std::endl;
          diagInfo.setDebugData(std::string("old_dt"), old_dt);
          diagInfo.setDebugData(std::string("last_dt"), last_dt);
        }

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
          Tools::pecletNumber(gv, parametersCon, dgfGradPot);

          Real minMembPot = std::numeric_limits<Real>::max();
          Real maxMembPot = std::numeric_limits<Real>::lowest();
          typename Traits::DGF_CON::Traits::RangeType minConCytosol(std::numeric_limits<Real>::max());
          typename Traits::DGF_CON::Traits::RangeType maxConCytosol(std::numeric_limits<Real>::lowest());
          typename Traits::DGF_CON::Traits::RangeType minConExtra(std::numeric_limits<Real>::max());
          typename Traits::DGF_CON::Traits::RangeType maxConExtra(std::numeric_limits<Real>::lowest());
          for(typename Traits::Physics::MIterator mit = physics.mInteriorBegin(); mit != physics.mInteriorEnd(); ++mit)
          {
            typename Traits::DGF_POT::Traits::RangeType membPot;
            physics.getMembranePotential(*mit, dgfPot, membPot);

            typename Traits::DGF_CON::Traits::RangeType conCytosol, conExtra;
            physics.getMembraneConcentrationJump(*mit, dgfCon, conCytosol, conExtra);

            if(membPot < minMembPot) minMembPot = membPot;
            if(membPot > maxMembPot) maxMembPot = membPot;

            for(int j=0; j<conCytosol.size(); j++)
            {
              if(conCytosol[j] < minConCytosol[j]) minConCytosol[j] = conCytosol[j];
              if(conCytosol[j] > maxConCytosol[j]) maxConCytosol[j] = conCytosol[j];
              if(conExtra[j] < minConExtra[j]) minConExtra[j] = conExtra[j];
              if(conExtra[j] > maxConExtra[j]) maxConExtra[j] = conExtra[j];
            }

            if(std::abs(time-physics.getParams().tEquilibrium()) < 1e-6)
            {
              typename GoldmanGF::Traits::RangeType goldmanPotential;
              // Evaluate at cell center
              typename GoldmanGF::Traits::DomainType x(0.5);

              //Dune::PDELab::IntersectionGeometry<typename Traits::Physics::ElementIntersection> ig(*mit,-1);
              goldmanGF.evaluate(*mit, x, goldmanPotential);

              debug_verb << mit->geometry().global(x) << std::endl;
              debug_verb << "  Membrane potential: " << physics.convertTo_mV(membPot) << std::endl;
              debug_verb << "  Goldman potential:  " << goldmanPotential << std::endl;
            }
          }
          minMembPot = gv.comm().min(minMembPot);
          maxMembPot = gv.comm().max(maxMembPot);

          debug_info << "MAX membrane potential = " << physics.convertTo_mV(maxMembPot) << " mV" << std::endl;
          debug_info << "MIN membrane potential = " << physics.convertTo_mV(minMembPot) << " mV" << std::endl;
          debug_info << "MIN/MAX membrane concentrations (CY) = " << minConCytosol << " / "
              << maxConCytosol << std::endl;
          debug_info << "MIN/MAX membrane concentrations (ES) = " << minConExtra << " / "
              << maxConExtra << std::endl;

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

        // Do not overwrite uold before output, otherwise Mori flux will be calculated wrongly.
        // So we update the vector as the very last instruction in the time loop
        uold = unew;

        // Same goes for the membrane flux states: The old values are needed (at least for Mori flux)
        // during output in order to correctly calculated fluxes. Do not overwrite old values before
        // output has completed!
        debug_jochen << "# gfMembFlux.acceptTimeStep" << std::endl;
        gfMembFlux.acceptTimeStep();
        debug_jochen << "# gfMoriFlux.acceptTimeStep" << std::endl;
        gfMoriFlux.acceptTimeStep();

      }
      // ==============================================================================================

      debug_info << "Total number of time steps: " << totalTimeSteps << std::endl;
    }


};


#endif /* DUNE_AX1_ACME2CYL_MORI_SIMULATION_HH */
