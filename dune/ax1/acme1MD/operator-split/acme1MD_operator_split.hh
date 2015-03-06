/*
 * acme1MD_operator_split.hh
 *
 *  Created on: Dec 15, 2011
 *      Author: jpods
 */

#ifndef ACME1MD_OPERATOR_SPLIT_HH
#define ACME1MD_OPERATOR_SPLIT_HH

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>

#include <dune/ax1/common/constants.hh>


template<typename Acme1MDTraits,
         typename DIRICHLET_GF_CON, typename DIRICHLET_GF_POT,
         typename GO_CON, typename GO_POT,
         typename SOLVER_CON, typename SOLVER_POT
         >
class Acme1MDOperatorSplit
{
  public:

    typedef Acme1MDTraits Traits;

    typedef Ax1VectorDiscreteGridFunction
        <typename Traits::Physics,typename Traits::GFS_CON,typename Traits::U_CON> DGF_CON;

    typedef Dune::PDELab::DiscreteGridFunction
        <typename Traits::GFS_POT,typename Traits::U_POT> DGF_POT;

    static const bool testPotentialOnly = true;
    static const bool useAdaptiveTimeStep = false;

    Acme1MDOperatorSplit(typename Traits::U_CON& uCon_, typename Traits::U_CON& uoldCon_,
          typename Traits::U_POT& uPot_,
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
          typename Traits::ACME1MD_OUTPUT& acme1MDOutput_)
      : uCon(uCon_),
        uoldCon(uoldCon_),
        uPot(uPot_),
        uoldCon_Inside(uoldCon_Inside_),
        unewCon_Inside(unewCon_Inside_),
        uoldCon_Outside(uoldCon_Outside_),
        unewCon_Outside(unewCon_Outside_),
        uPot_Inside(uPot_Inside_),
        uPot_Outside(uPot_Outside_),
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
        acme1MDOutput(acme1MDOutput_),
        diagInfo(acme1MDOutput_.getDiagInfo()),
        dgfCon(physics_,gfsCon_,uCon_),
        dgfPot(gfsPot_,uPot_),
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

      uoldPot = uPot;

      // Handle potential constraints; concentrations constraints are handled by instationary gridoperator
      dirichletValuesPot.setTime(time+dt);
      // Apply computed constraints to coefficient vector
      Dune::PDELab::interpolate(dirichletValuesPot,gfsPot,uPot);
      // Copy non-constrained DOFs (although this is actually unnecessary, they are overwritten anyways)
      Dune::PDELab::copy_nonconstrained_dofs(solverPot.getGridOperator().localAssembler().trialConstraints(),
          uoldPot, uPot);

      //Output::printMultipleComponentCoefficientVector(uoldCon, NUMBER_OF_SPECIES);
      //Output::printMultipleComponentCoefficientVector(uCon, NUMBER_OF_SPECIES);


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
      ++iterations;

      //debug_verb << "Iteration #" << iterations << std::endl;

      //debug_verb << "================= Solving Poisson equation..." << std::endl;
      typename Traits::U_POT tempPot(uPot);
      tempPot *= (1-gammaPot);

      // Apply computed constraints to coefficient vector; forget other DOFs, this is a stationary problem
      Dune::PDELab::interpolate(dirichletValuesPot,gfsPot,uPot);
      solverPot.apply(uPot);

      uPot *= gammaPot;
      uPot += tempPot;

      // Check change in potential with respect to last time step and adapt time step
      if(useAdaptiveTimeStep)
      {
        dt = calculateTimeStep(dt);
      }
      //debug_verb << "================= Solved Poisson equation" << std::endl;

      // Transfer potential solution to subgrids
      Traits::UPOT_Restrictor::restrict(gv, subGV_Inside, gfsPot, gfsPot_Inside, uPot, uPot_Inside);
      Traits::UPOT_Restrictor::restrict(gv, subGV_Outside, gfsPot, gfsPot_Outside, uPot, uPot_Outside);

      //debug_verb << "================= Solving Nernst-Planck equation..." << std::endl;
      typename Traits::USUB_CON tempCon_Inside(uoldCon_Inside);
      tempCon_Inside *= (1-gammaCon);
      typename Traits::USUB_CON tempCon_Outside(uoldCon_Outside);
      tempCon_Outside *= (1-gammaCon);

      solverCon_Inside.apply(time,dt,uoldCon_Inside,dirichletValuesCon_Inside,unewCon_Inside);
      solverCon_Outside.apply(time,dt,uoldCon_Outside,dirichletValuesCon_Outside,unewCon_Outside);

      // Damped update: unewCon = (gamma-1)*uOldCon + gamma*unewCon;
      unewCon_Inside *= gammaCon;
      unewCon_Inside += tempCon_Inside;
      unewCon_Outside *= gammaCon;
      unewCon_Outside += tempCon_Outside;
      //debug_verb << "================= Solved Nernst-Planck equation" << std::endl;

      uoldCon = uCon;

      // Transfer concentration solutions to host grid
      Traits::UCON_Interpolator::interpolate(subGV_Inside, gv, gfsCon_Inside, gfsCon, unewCon_Inside, uCon);
      Traits::UCON_Interpolator::interpolate(subGV_Outside, gv, gfsCon_Outside, gfsCon, unewCon_Outside, uCon);

      //Output::printMultipleComponentCoefficientVector(uCon, NUMBER_OF_SPECIES);

      debug_verb << "-----------------------------------------------------------------" << std::endl;

      // Fill diagInfo
      diagInfo.iterations = iterations;
      diagInfo.dt = dt;

      return dt;
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
                     "Acme1MDOperatorSplit::calcDefect(): Non-linear defect is NaN or Inf");

    }


    typename Traits::Real calculateTimeStep(typename Traits::Real dt)
    {
      typename Traits::Real absPotChange(0.0), relPotChange(0.0);
      typename DGF_POT::Traits::DomainType xAbs(0.0), xRel(0.0);

      calcMaxDiffPot(absPotChange, relPotChange, xAbs, xRel);

      debug_info << "================================" << std::endl;
      debug_verb << "== max pot change: " << absPotChange << " [" << xAbs << "]"
          << " (" << relPotChange << " [" << xRel << "])" << std::endl;

      debug_info << "== dt: " << dt;
      // Branch "increase time step", maximum 100ns
      if(relPotChange < reduction /* || absPotChange < absLimit*/)
      {
        dt *= 1.7;
        dt = std::min(dt, 100.0);
      // Branch "decrease time step", minimum 0.1ns
      } else {
        dt /= 2;
        dt = std::max(dt, 0.1);
      }

      debug_verb  << " ==> " << dt << std::endl;
      debug_info << "================================" << std::endl;

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
    typename Traits::U_CON& uoldCon;
    typename Traits::U_POT& uPot;

    // Reference to solution vectors (subdomains)
    typename Traits::USUB_CON& uoldCon_Inside;
    typename Traits::USUB_CON& unewCon_Inside;
    typename Traits::USUB_CON& uoldCon_Outside;
    typename Traits::USUB_CON& unewCon_Outside;
    typename Traits::USUB_POT& uPot_Inside;
    typename Traits::USUB_POT& uPot_Outside;

    // Solution vectors for saving the state of last the iteration / time step =>  copies, no references!
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

    typename Traits::ACME1MD_OUTPUT& acme1MDOutput;
    typename Traits::ACME1MD_OUTPUT::DiagnosticInfo& diagInfo;

    DGF_CON dgfCon;
    DGF_POT dgfPot;
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


#endif /* ACME1MD_OPERATOR_SPLIT_HH */
