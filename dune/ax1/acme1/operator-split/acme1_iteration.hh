/*
 * acme1_operator_split.hh
 *
 *  Created on: Dec 15, 2011
 *      Author: jpods
 */

#ifndef ACME1_ITERATION_HH
#define ACME1_ITERATION_HH

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>

#include <dune/ax1/common/constants.hh>


template<typename Acme1Traits,
         typename DIRICHLET_GF_CON, typename DIRICHLET_GF_POT,
         typename GO_CON, typename GO_POT,
         typename SOLVER_CON, typename SOLVER_POT
         >
class Acme1Iteration
{
  public:

    typedef Acme1Traits Traits;

    typedef Ax1VectorDiscreteGridFunction
        <typename Traits::Physics,typename Traits::GFS_CON,typename Traits::U_CON> DGF_CON;

    typedef Dune::PDELab::DiscreteGridFunction
        <typename Traits::GFS_POT,typename Traits::U_POT> DGF_POT;

    static const bool testPotentialOnly = true;

    Acme1Iteration(typename Traits::U_CON& uCon_, typename Traits::U_POT& uPot_,
          typename Traits::USUB_CON& uoldCon_Inside_, typename Traits::USUB_CON& unewCon_Inside_,
          typename Traits::USUB_CON& uoldCon_Outside_, typename Traits::USUB_CON& unewCon_Outside_,
          typename Traits::USUB_POT& uPot_Inside_, typename Traits::USUB_POT& uPot_Outside_,
          DIRICHLET_GF_CON& dirichletValuesCon_Inside_, DIRICHLET_GF_CON& dirichletValuesCon_Outside_,
          DIRICHLET_GF_POT& dirichletValuesPot_,
          GO_CON& goCon_Inside_, GO_CON& goCon_Outside_, GO_POT& goPot_,
          SOLVER_CON& solverCon_Inside_, SOLVER_CON& solverCon_Outside_, SOLVER_POT& solverPot_,
          typename Traits::GFS_CON& gfsCon_, typename Traits::GFS_POT& gfsPot_,
          typename Traits::SUB_GFS_CON& gfsCon_Inside_, typename Traits::SUB_GFS_CON& gfsCon_Outside_,
          typename Traits::SUB_GFS_POT& gfsPot_Inside_, typename Traits::SUB_GFS_POT& gfsPot_Outside_,
          typename Traits::Physics& physics_, const typename Traits::GridView& gv_,
          const typename Traits::SubGridView& subGridView_Inside_,
          const typename Traits::SubGridView& subGridView_Outside_,
          typename Traits::ACME1_OUTPUT& acme1Output_)
      : uCon(uCon_),
        uPot(uPot_),
        uoldCon_Inside(uoldCon_Inside_),
        unewCon_Inside(unewCon_Inside_),
        uoldCon_Outside(uoldCon_Outside_),
        unewCon_Outside(unewCon_Outside_),
        uPot_Inside(uPot_Inside_),
        uPot_Outside(uPot_Outside_),
        uConPrevious_Inside(uoldCon_Inside_),
        uConPrevious_Outside(uoldCon_Outside_),
        uPotPrevious(uPot_),
        uConPrevious(uCon_),
        uoldPot(uPot_),
        dirichletValuesCon_Inside(dirichletValuesCon_Inside_),
        dirichletValuesCon_Outside(dirichletValuesCon_Outside_),
        dirichletValuesPot(dirichletValuesPot_),
        gridoperatorCon_Inside(goCon_Inside_),
        gridoperatorCon_Outside(goCon_Outside_),
        gridoperatorPot(goPot_),
        solverCon_Inside(solverCon_Inside_),
        solverCon_Outside(solverCon_Outside_),
        solverPot(solverPot_),
        gfsCon(gfsCon_),
        gfsPot(gfsPot_),
        gfsCon_Inside(gfsCon_Inside_),
        gfsCon_Outside(gfsCon_Outside_),
        gfsPot_Inside(gfsPot_Inside_),
        gfsPot_Outside(gfsPot_Outside_),
        physics(physics_),
        gv(gv_),
        subGV_Inside(subGridView_Inside_),
        subGV_Outside(subGridView_Outside_),
        acme1Output(acme1Output_),
        diagInfo(acme1Output_.getDiagInfo()),
        dgfCon(physics_,gfsCon_,uCon_),
        dgfConPrevious(physics_,gfsCon_,uConPrevious),
        dgfPot(gfsPot_,uPot_),
        dgfPotPrevious(gfsPot_,uPotPrevious),
        dgfOldPot(gfsPot_,uoldPot),
        reduction(1e-8),
        absLimit(1e-12),
        maxIterations(40),
        // we need at least one iteration when we're using changes
        // in the coefficients as the termination criterion
        useDefect(true),
        forceIteration(!useDefect),
        initialDefectCon_Inside(1e100),
        initialDefectCon_Outside(1e100),
        initialDefectPot(1e100),
        initialMaxDiffCon(1e100),
        initialMaxDiffPot(1e100),
        iterations(0),
        debugOutput(false)
    {
    }


    //! \brief Übelst geile Iteration
    typename Traits::Real timeStep(typename Traits::Real time, typename Traits::Real dt)
    {
      iterations = 0;
      if(Traits::SETUP::writeIfNotConverged &&
          (time > diagInfo.tEquilibrium || std::abs(time-diagInfo.tEquilibrium) < 1e-8))
      {
        debugOutput = true;
      }

      uConPrevious_Inside = unewCon_Inside;
      uConPrevious_Outside = unewCon_Outside;
      uConPrevious = uCon;
      uPotPrevious = uPot;
      uoldPot = uPot;

      // Handle potential constraints; concentrations constraints are handled by instationary gridoperator
      dirichletValuesPot.setTime(time+dt);
      // Apply computed constraints to coefficient vector
      Dune::PDELab::interpolate(dirichletValuesPot,gfsPot,uPot);
      // Copy non-constrained DOFs (although this is actually unnecessary, they are overwritten anyways)
      Dune::PDELab::copy_nonconstrained_dofs(solverPot.getGridOperator().localAssembler().trialConstraints(),
          uPotPrevious, uPot);


      if(useDefect)
      {
        // TODO Find a better way how to get the initial residual
        // Possible solution: Call instationary GO's preStep() and preStage() methods once before
        // calculating the initial residual and clean up with postStage(), postStep() afterwards!
        //debug_verb << "=== Do one dummy step for NernstPlanck solver in order to get initial defect ===" << std::endl;

        // Workaround to call necessary assembler calls in instationary gridoperator
        solverCon_Inside.apply(time,dt,uoldCon_Inside,dirichletValuesCon_Inside,unewCon_Inside);
        solverCon_Outside.apply(time,dt,uoldCon_Outside,dirichletValuesCon_Outside,unewCon_Outside);

        unewCon_Inside = uoldCon_Inside; // reset unewCon_Inside
        unewCon_Outside = uoldCon_Outside; // reset unewCon_Inside
        //debug_verb << "================================================================================" << std::endl;
      }

      // Damping parameters
      typename Traits::Real gammaCon = physics.getParams().getGammaCon();
      typename Traits::Real gammaPot = physics.getParams().getGammaPot();

      // Begin solving decoupled system
      while(not converged())
      {
        ++iterations;

        //debug_verb << "Iteration #" << iterations << std::endl;

        uConPrevious_Inside = unewCon_Inside;
        uConPrevious_Outside = unewCon_Outside;
        uConPrevious = uCon;
        uPotPrevious = uPot;



        //debug_verb << "================= Solving Poisson equation..." << std::endl;
        typename Traits::U_POT tempPot(uPot);
        tempPot *= (1-gammaPot);

        // Apply computed constraints to coefficient vector; forget other DOFs, this is a stationary problem
        Dune::PDELab::interpolate(dirichletValuesPot,gfsPot,uPot);
        solverPot.apply(uPot);

        uPot *= gammaPot;
        uPot += tempPot;

        // Check change in potential with respect to last time step and adapt time step
        /*
        if(iterations==1 && time>diagInfo.tEquilibrium)
        {
          dt = calculateTimeStep(dt);
          debug_info << "================================" << std::endl;
          debug_info << "Using time step dt = " << dt << std::endl;
          debug_info << "================================" << std::endl;
        }*/

        //debug_verb << "================= Solved Poisson equation" << std::endl;

        // Transfer potential solution to subgrids
        Traits::UPOT_Restrictor::restrict(gv, subGV_Inside, gfsPot, gfsPot_Inside, uPot, uPot_Inside);
        Traits::UPOT_Restrictor::restrict(gv, subGV_Outside, gfsPot, gfsPot_Outside, uPot, uPot_Outside);


        //debug_verb << "================= Solving Nernst-Planck equation..." << std::endl;
        typename Traits::USUB_CON tempCon_Inside(uConPrevious_Inside);
        tempCon_Inside *= (1-gammaCon);
        typename Traits::USUB_CON tempCon_Outside(uConPrevious_Outside);
        tempCon_Outside *= (1-gammaCon);

        solverCon_Inside.apply(time,dt,uoldCon_Inside,dirichletValuesCon_Inside,unewCon_Inside);
        solverCon_Outside.apply(time,dt,uoldCon_Outside,dirichletValuesCon_Outside,unewCon_Outside);

        // Damped update: unewCon = (gamma-1)*uOldCon + gamma*unewCon;
        unewCon_Inside *= gammaCon;
        unewCon_Inside += tempCon_Inside;
        unewCon_Outside *= gammaCon;
        unewCon_Outside += tempCon_Outside;
        //debug_verb << "================= Solved Nernst-Planck equation" << std::endl;

        // "Real" defect for potential using the old concentrations
        //typename Traits::U_POT rPot(gridoperatorPot.testGridFunctionSpace());
        //rPot = 0.0;
        //gridoperatorPot.residual(uPot, rPot);
        //typename Traits::Real defectPot = solverPot.getLinearSolver().norm(rPot);
        //debug_verb << "BEFORE:                   "    << defectPot << std::endl;

        // Transfer concentration solutions to host grid
        Traits::UCON_Interpolator::interpolate(subGV_Inside, gv, gfsCon_Inside, gfsCon, unewCon_Inside, uCon);
        Traits::UCON_Interpolator::interpolate(subGV_Outside, gv, gfsCon_Outside, gfsCon, unewCon_Outside, uCon);


//        Output::printSingleCoefficientVector(uPotPrevious,"uPotPrevious");
//        Output::printSingleCoefficientVector(uPot,"uPot");
//        Output::printSingleCoefficientVector(uPot_Inside,"uPot_Inside");
//        Output::printSingleCoefficientVector(uPot_Outside,"uPot_Outside");
//
//        Output::printMultipleComponentCoefficientVector(uConPrevious,NUMBER_OF_SPECIES);
//        Output::printMultipleComponentCoefficientVector(uConPrevious_Inside,NUMBER_OF_SPECIES);
//        Output::printMultipleComponentCoefficientVector(uConPrevious_Outside,NUMBER_OF_SPECIES);
//        Output::printMultipleComponentCoefficientVector(uCon,NUMBER_OF_SPECIES);
//        Output::printMultipleComponentCoefficientVector(unewCon_Inside,NUMBER_OF_SPECIES);
//        Output::printMultipleComponentCoefficientVector(unewCon_Outside,NUMBER_OF_SPECIES);

        //acme1Output.writeStep(time+dt);
      }
      debug_verb << "-----------------------------------------------------------------" << std::endl;

      // Fill diagInfo
      diagInfo.iterations = iterations;
      diagInfo.dt = dt;

      return dt;
    }

    bool converged()
    {
      if(iterations == 0 && forceIteration)
      {
        debug_verb << "                  CHANGE IN SOLUTION                           DEFECT            " << std::endl;
        debug_verb << "   abs pot  (rel pot)        abs con  (rel con)              pot  (con[in],      con[out])" << std::endl;

        //typename Traits::Real defectCon_Inside = 1e100, defectCon_Outside = 1e100, defectPot = 1e100;
        //calcDefect(defectCon_Inside, defectCon_Outside, defectPot);

        //debug_verb << "FIRST:                                                "
        //  << defectPot << "  (" << defectCon_Inside << ", " << defectCon_Outside << ")" << std::endl;

        return false;
      }

      // This should actually never be called; iterations==maxIterations is tested for below
      if(iterations > maxIterations)
      {
        DUNE_THROW(Dune::Exception, "Did not converge after " << iterations << " iterations");
      }

      if(useDefect)
      {
        typename Traits::Real defectCon_Inside = 1e100, defectCon_Outside = 1e100, defectPot = 1e100;

        calcDefect(defectCon_Inside, defectCon_Outside, defectPot);
        if(iterations == 0)
        {
          initialDefectCon_Inside = defectCon_Inside;
          initialDefectCon_Outside = defectCon_Outside;
          initialDefectPot = defectPot;
        }
        //debug_verb << "  defectCon = " << defectCon << std::endl;
        //debug_verb << "  defectPot = " << defectPot << std::endl;

        // Both concentrations and potential need to satisfy either the
        // absolute or the relative error tolerance criterion
        if((defectCon_Inside < absLimit || defectCon_Inside < initialDefectCon_Inside * reduction)
            && (defectCon_Outside < absLimit || defectCon_Outside < initialDefectCon_Outside * reduction)
            && (defectPot < absLimit || defectPot < initialDefectPot * reduction))
        {
          return true;
        } else {
          if(iterations == maxIterations)
          {
            DUNE_THROW(Dune::Exception, "Did not converge after " << iterations << " iterations");
          }
          return false;
        }
      }
      else // Use max change in any unknown over the entire domain
      {
        Dune::ios_base_all_saver ding(std::cout);
        debug_verb << std::setprecision(4) << std::scientific;

        typename Traits::Real absDiffCon(0.0), absDiffPot(0.0), relDiffCon(0.0), relDiffPot(0.0);
        calcMaxDiffConPot(absDiffCon, absDiffPot, relDiffCon, relDiffPot);

        // Potential changes with respect to last time step
        typename Traits::Real absPotChange(0.0), relPotChange(0.0);
        typename DGF_POT::Traits::DomainType xAbs(0.0), xRel(0.0);
        calcMaxDiffPot(absPotChange, relPotChange, xAbs, xRel);

        // For testing purposes only
        typename Traits::Real defectCon_Inside = 1e100, defectCon_Outside = 1e100, defectPot = 1e100;
        //calcDefect(defectCon_Inside, defectCon_Outside, defectPot);

        // Changes with respect to last iteration
        debug_verb << absDiffPot << "  (" << relDiffPot << ")  ";
        debug_verb << absDiffCon << "  (" << relDiffCon << ")";

        debug_verb << " || " << defectPot << "  (" << defectCon_Inside << ", " << defectCon_Outside << ")"
            << std::endl;

        // Changes with respect to last time step
        debug_verb << absPotChange << "  (" << relPotChange << ")"
            << " [" << xAbs << " (" << xRel << ")]" << std::endl;

        // maximum change in the solution is given after iteration 0, in the beginning of the 1st
        //if(iterations == 1)
        //{
        //  initialMaxDiffCon = absDiffCon;
        //  initialMaxDiffPot = absDiffPot;
        //}

        //debug_verb << "maxDiffCon = " << maxDiffCon << std::endl;
        //debug_verb << "maxDiffPot = " << maxDiffPot << std::endl;

        // It should be sufficient to test for changes in the potential only!
        //if(testPotentialOnly &&
        //    (relDiffPot < reduction || absDiffPot < absLimit))
        //{
        //  return true;
        //}

        // Use 'reduction' as error tolerance for relative changes, 'absLimit' for absolute changes
        if((relDiffCon < reduction || absDiffCon < absLimit)
          && (relDiffPot < reduction || absDiffPot < absLimit))
        {
          return true;
        } else {
          if(iterations == maxIterations)
          {
            DUNE_THROW(Dune::Exception, "Did not converge after " << iterations << " iterations");
          }
          return false;
        }
      }
      return false;
    }

    void calcDefect(typename Traits::Real& defectCon_Inside, typename Traits::Real& defectCon_Outside,
        typename Traits::Real& defectPot) const
    {
      typename Traits::USUB_CON rCon_Inside(gridoperatorCon_Inside.testGridFunctionSpace());
      rCon_Inside = 0.0;
      typename Traits::USUB_CON rCon_Outside(gridoperatorCon_Outside.testGridFunctionSpace());
      rCon_Outside = 0.0;

      // Internally already uses the updated potential
      gridoperatorCon_Inside.residual(unewCon_Inside, rCon_Inside);
      gridoperatorCon_Outside.residual(unewCon_Outside, rCon_Outside);

      // It doesn't matter which linear solver's norm we take, as the
      // current implementation of seqistlsolverbackend always take the 2-norm
      //Real defectCon = solverCon.getPDESolver().getLinearSolver().norm(rCon);
      defectCon_Inside = solverPot.getLinearSolver().norm(rCon_Inside);
      defectCon_Outside = solverPot.getLinearSolver().norm(rCon_Outside);

      typename Traits::U_POT rPot(gridoperatorPot.testGridFunctionSpace());
      rPot = 0.0;
      gridoperatorPot.residual(uPot, rPot);

      defectPot = solverPot.getLinearSolver().norm(rPot);

      if (!std::isfinite(defectCon_Inside) || !std::isfinite(defectCon_Outside)
            ||!std::isfinite(defectPot))
          DUNE_THROW(Dune::Exception,
                     "Acme1OperatorSplit::calcDefect(): Non-linear defect is NaN or Inf");

    }


    //! This method calculated the maximum difference between old and new solution of potential/concentrations
    //! It uses the whole domain for calculations, i.e. it doesn't distinguish between concentration subdomains!
    void calcMaxDiffConPot(typename Traits::Real& absMaxDiffCon, typename Traits::Real& absMaxDiffPot,
        typename Traits::Real& relMaxDiffCon, typename Traits::Real& relMaxDiffPot) const
    {
      // ================== Calculate max conc/pot change (works in 1D only!) ===========================
      typedef typename Traits::GridView::template Codim<0>::Entity Entity;
      typedef typename Traits::GridView::template Codim<0>::EntityPointer EntityPointer;
      typedef typename Traits::GridView::template Codim<0>::Iterator ElementLeafIterator;

      typename DGF_CON::Traits::DomainType x;

      typename DGF_CON::Traits::RangeType yoldCon;
      typename DGF_CON::Traits::RangeType ynewCon;

      typename DGF_POT::Traits::RangeType yoldPot;
      typename DGF_POT::Traits::RangeType ynewPot;

      absMaxDiffCon = 0.0;
      absMaxDiffPot = 0.0;
      relMaxDiffCon = 0.0;
      relMaxDiffPot = 0.0;

      for (ElementLeafIterator eit=gv.template begin<0>(); eit!=gv.template end<0>(); ++eit)
      {
        Entity& e = *eit;

        x = e.geometry().local(e.geometry().center());

        // =============== Concentration part ========================================
        // Calculate concentrations differences only on non-membrane regions!
        if(not physics.isMembrane(e))
        {
          dgfConPrevious.evaluate(e, x, yoldCon);
          dgfCon.evaluate(e, x, ynewCon);
          typename Traits::Real absDiffCon = 0.0;
          typename Traits::Real relDiffCon = 0.0;


          for(int j=0; j<NUMBER_OF_SPECIES; ++j)
          {
            // Relative change in concentrations
            // Check for division by 0
            if(std::abs(yoldCon[j]) > 1e-12)
            {
              relDiffCon = std::abs(ynewCon[j] / yoldCon[j] - 1.0);
            } else {
              // Exclude this relative change as it can not be calculated
              //diffCon = std::abs(ynewCon[j] - yoldCon[j]);
            }

            // Absolute change
            absDiffCon = std::abs(ynewCon[j] - yoldCon[j]);

            absMaxDiffCon = std::max(absMaxDiffCon, absDiffCon);
            relMaxDiffCon = std::max(relMaxDiffCon, relDiffCon);
          }
        }

        // =============== Potential part ============================================
        dgfPotPrevious.evaluate(e, x, yoldPot);
        dgfPot.evaluate(e, x, ynewPot);

        typename Traits::Real absDiffPot = 0.0;
        typename Traits::Real relDiffPot = 0.0;

        // Relative change in potential
        // Check for division by 0
        if(std::abs(yoldPot[0]) > 1e-12)
        {
          relDiffPot = std::abs(ynewPot[0] / yoldPot[0] - 1.0);
        } else {
          // Exlude this relative change as it can not be calculated
          //diffPot = std::abs(ynewPot[0] - yoldPot[0]);
        }

        // Absolute errors
        absDiffPot = std::abs(ynewPot[0] - yoldPot[0]);

        absMaxDiffPot = std::max(absMaxDiffPot, absDiffPot);
        relMaxDiffPot = std::max(relMaxDiffPot, relDiffPot);
      }

      if(debugOutput)
      {
        diagInfo.setMaxDiffConPot(absMaxDiffCon, absMaxDiffPot);
      }
    }

    typename Traits::Real calculateTimeStep(typename Traits::Real dt)
    {
      typename Traits::Real absPotChange(0.0), relPotChange(0.0);
      typename DGF_POT::Traits::DomainType xAbs(0.0), xRel(0.0);

      calcMaxDiffPot(absPotChange, relPotChange, xAbs, xRel);
      debug_verb << "== max pot change: " << absPotChange << " [" << xAbs << "]"
          << " (" << relPotChange << " [" << xRel << "])" << std::endl;

      // Branch "increase time step", maximum 100ns
      if(relPotChange < reduction /* || absPotChange < absLimit*/)
      {
        dt *= 2;
        dt = std::min(dt, 100.0);
      // Branch "decrease time step", minimum 0.1ns
      } else {
        dt /= 2;
        dt = std::max(dt, 0.1);
      }

      return dt;
    }

    //! This method calculated the maximum difference between old time step and new solution of potential
    void calcMaxDiffPot(typename Traits::Real& absMaxDiffPot, typename Traits::Real& relMaxDiffPot,
        typename DGF_POT::Traits::DomainType& xAbs, typename DGF_POT::Traits::DomainType& xRel) const
    {
      // ================== Calculate max conc/pot change (works in 1D only!) ===========================
      typedef typename Traits::GridView::template Codim<0>::Entity Entity;
      typedef typename Traits::GridView::template Codim<0>::EntityPointer EntityPointer;
      typedef typename Traits::GridView::template Codim<0>::Iterator ElementLeafIterator;

      typename DGF_POT::Traits::DomainType x;

      typename DGF_POT::Traits::RangeType yoldPot;
      typename DGF_POT::Traits::RangeType ynewPot;

      absMaxDiffPot = 0.0;
      relMaxDiffPot = 0.0;

      for (ElementLeafIterator eit=gv.template begin<0>(); eit!=gv.template end<0>(); ++eit)
      {
        Entity& e = *eit;
        x = e.geometry().local(e.geometry().center());

        // =============== Potential part ============================================
        dgfOldPot.evaluate(e, x, yoldPot);
        dgfPot.evaluate(e, x, ynewPot);

        typename Traits::Real absDiffPot = 0.0;
        typename Traits::Real relDiffPot = 0.0;

        // Relative change in potential
        // Check for division by 0
        if(std::abs(yoldPot[0]) > 1e-12)
        {
          relDiffPot = std::abs(ynewPot[0] / yoldPot[0] - 1.0);
        }

        // Absolute errors
        absDiffPot = std::abs(ynewPot[0] - yoldPot[0]);

        if(absDiffPot > absMaxDiffPot)
        {
          absMaxDiffPot = absDiffPot;
          xAbs = e.geometry().center();
        }
        if(relDiffPot > relMaxDiffPot)
        {
          relMaxDiffPot = relDiffPot;
          xRel = e.geometry().center();
        }
      }
    }



    void setReduction(typename Traits::Real reduction_)
    {
      reduction = reduction_;
    }

    void setAbsLimit(typename Traits::Real absLimit_)
    {
      absLimit = absLimit_;
    }

    void setMaxIterations(int maxIterations_)
    {
      maxIterations = maxIterations_;
    }

    void setUseDefect(bool useDefect_)
    {
      useDefect = useDefect_;
      if(!useDefect_) forceIteration = true;
    }

    void setForceIteration(bool forceIteration_)
    {
      if(useDefect) forceIteration = forceIteration_;
    }


  private:

    // References to solution vectors (whole domain)
    typename Traits::U_CON& uCon;
    typename Traits::U_POT& uPot;

    // Reference to solution vectors (subdomains)
    typename Traits::USUB_CON& uoldCon_Inside;
    typename Traits::USUB_CON& unewCon_Inside;
    typename Traits::USUB_CON& uoldCon_Outside;
    typename Traits::USUB_CON& unewCon_Outside;
    typename Traits::USUB_POT& uPot_Inside;
    typename Traits::USUB_POT& uPot_Outside;

    // Solution vectors for saving the state of last the iteration / time step =>  copies, no references!
    typename Traits::USUB_CON uConPrevious_Inside;   // copy of uoldCon, updated each iteration
    typename Traits::USUB_CON uConPrevious_Outside;   // copy of uoldCon, updated each iteration
    typename Traits::U_POT uPotPrevious;   // copy of uPot, updated each iteration
    typename Traits::U_CON uConPrevious;   // copy of uCon, updated each iteration
    typename Traits::U_POT uoldPot;        // copy of uPot representing value after last time step

    DIRICHLET_GF_CON& dirichletValuesCon_Inside;
    DIRICHLET_GF_CON& dirichletValuesCon_Outside;
    DIRICHLET_GF_POT& dirichletValuesPot;

    GO_CON& gridoperatorCon_Inside;
    GO_CON& gridoperatorCon_Outside;
    GO_POT& gridoperatorPot;

    SOLVER_CON& solverCon_Inside;
    SOLVER_CON& solverCon_Outside;
    SOLVER_POT& solverPot;

    typename Traits::GFS_CON& gfsCon;
    typename Traits::GFS_POT& gfsPot;
    typename Traits::SUB_GFS_CON& gfsCon_Inside;
    typename Traits::SUB_GFS_CON& gfsCon_Outside;
    typename Traits::SUB_GFS_POT& gfsPot_Inside;
    typename Traits::SUB_GFS_POT& gfsPot_Outside;

    typename Traits::Physics& physics;
    const typename Traits::GridView& gv;
    const typename Traits::SubGridView& subGV_Inside;
    const typename Traits::SubGridView& subGV_Outside;

    typename Traits::ACME1_OUTPUT& acme1Output;
    typename Traits::ACME1_OUTPUT::DiagnosticInfo& diagInfo;

    DGF_CON dgfCon;
    DGF_POT dgfPot;
    DGF_CON dgfConPrevious;
    DGF_POT dgfPotPrevious;
    DGF_POT dgfOldPot;

    typename Traits::Real reduction;
    typename Traits::Real absLimit;

    int maxIterations;
    bool useDefect;
    bool forceIteration;

    typename Traits::Real initialDefectCon_Inside;
    typename Traits::Real initialDefectCon_Outside;
    typename Traits::Real initialDefectPot;

    typename Traits::Real initialMaxDiffCon;
    typename Traits::Real initialMaxDiffPot;

    int iterations;

    bool debugOutput;
};


#endif /* ACME1_ITERATION_HH */
