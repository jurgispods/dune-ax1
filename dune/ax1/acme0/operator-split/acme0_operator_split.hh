/*
 * acme0_operator_split.hh
 *
 *  Created on: Dec 15, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME0_OPERATOR_SPLIT_HH
#define DUNE_AX1_ACME0_OPERATOR_SPLIT_HH

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>

#include <dune/ax1/common/constants.hh>

template<typename U_CON, typename U_POT,
         typename DIRICHLET_GF_CON, typename DIRICHLET_GF_POT,
         typename SOLVER_CON, typename SOLVER_POT,
         typename GFS_CON, typename GFS_POT,
         typename PHYSICS, typename GV,
         typename GO_CON, typename GO_POT, typename Real>
class Acme0OperatorSplit
{
  public:

    typedef Dune::PDELab::VectorDiscreteGridFunction<GFS_CON,U_CON> DGF_CON;
    typedef Dune::PDELab::DiscreteGridFunction<GFS_POT,U_POT>       DGF_POT;

    static const bool calcRelativeErrors = false;

    Acme0OperatorSplit(U_CON& uoldCon_, U_CON& unewCon_,
         U_POT& uPot_, DIRICHLET_GF_CON& dirichletValuesCon_, DIRICHLET_GF_POT& dirichletValuesPot_,
         SOLVER_CON& solverCon_, SOLVER_POT& solverPot_,
         GFS_CON& gfsCon_, GFS_POT& gfsPot_,
         PHYSICS& physics_, GV& gv_, GO_CON& gridoperatorCon_,
         GO_POT& gridoperatorPot_)
      : uoldCon(uoldCon_),
        unewCon(unewCon_),
        uConPrevious(uoldCon),
        uPot(uPot_),
        uPotPrevious(uPot),
        dirichletValuesCon(dirichletValuesCon_),
        dirichletValuesPot(dirichletValuesPot_),
        solverCon(solverCon_),
        solverPot(solverPot_),
        gfsCon(gfsCon_),
        gfsPot(gfsPot_),
        dgfCon(gfsCon,unewCon),
        dgfConPrevious(gfsCon,uConPrevious),
        dgfPot(gfsPot,uPot),
        dgfPotPrevious(gfsPot,uPotPrevious),
        physics(physics_),
        gv(gv_),
        gridoperatorCon(gridoperatorCon_),
        gridoperatorPot(gridoperatorPot_),
        reduction(1e-8),
        absLimit(1e-12),
        maxIterations(40),
        // we need at least one iteration when we're using changes
        // in the coefficients as the termination criterion
        useDefect(true),
        forceIteration(!useDefect),
        initialDefectCon(1e100),
        initialDefectPot(1e100),
        iterations(0)
    {
      // Relative errors are not implemented in a satisfactory way
      assert(calcRelativeErrors == false);
    }


    //! \brief Übelst geile Iteration (nicht)
    template<typename DIAG_INFO>
    void timeStep(Real time, Real dt, DIAG_INFO& diagInfo)
    {
      iterations = 0;

      if(useDefect)
      {
        // TODO Find a better way how to get the initial residual
        // Workaround to call necessary assembler calls in instationary gridoperator
        solverCon.apply(time,dt,uoldCon,dirichletValuesCon,unewCon);
        unewCon = uoldCon; // reset unewCon
        //debug_verb << "Initial defect: " << initialDefect << std::endl;
      }

      // Damping parameters
      Real gammaCon = physics.getParams().getGammaCon();
      Real gammaPot = physics.getParams().getGammaPot();

      // Begin solving decoupled system
      while(not converged())
      {
        ++iterations;

        uConPrevious = unewCon;
        uPotPrevious = uPot;

        //debug_info << "iteration # " << iterations;

        //debug_verb << "================= Solving Poisson equation..." << std::endl;
        U_POT tempPot(uPot);
        tempPot *= (1-gammaPot);

        // Apply computed constraints to coefficient vector
        Dune::PDELab::interpolate(dirichletValuesPot,gfsPot,uPot);
        solverPot.apply(uPot);

        uPot *= gammaPot;
        uPot += tempPot;
        //debug_verb << "================= Solved Poisson equation" << std::endl;

        //debug_verb << "================= Solving Nernst-Planck equation..." << std::endl;
        U_CON tempCon(uConPrevious);
        tempCon *= (1-gammaCon);

        solverCon.apply(time,dt,uoldCon,dirichletValuesCon,unewCon);

        // Damped update: unewCon = (gamma-1)*uOldCon + gamma*unewCon;
        unewCon *= gammaCon;
        unewCon += tempCon;
        //debug_verb << "================= Solved Nernst-Planck equation" << std::endl;

        //Output::printCoefficientVectorDG(uPot,1);
        //Output::printCoefficientVectorDG(unewCon,NUMBER_OF_SPECIES);
      }
      //debug_verb << "-----------------------------------------------------------------" << std::endl;

      // Fill diagInfo
      diagInfo.iterations = iterations;
      diagInfo.dt = dt;
    }

    bool converged()
    {
      if(iterations == 0 && forceIteration) return false;

      if(iterations > maxIterations)
      {
        // This should actually never be called; iterations==maxIterations is tested for below
        DUNE_THROW(Dune::Exception, "Did not converge after " << iterations << " iterations");
      }

      if(useDefect)
      {
        Real defectCon = 1e100, defectPot = 1e100;
        calcDefect(defectCon, defectPot);
        if(iterations == 0)
        {
          initialDefectCon = defectCon;
          initialDefectPot = defectPot;
        }
        //debug_verb << "  defectCon = " << defectCon << std::endl;
        //debug_verb << "  defectPot = " << defectPot << std::endl;

        // Both concentrations and potential need to satisfy either the
        // absolute or the relative error tolerance criterion
        if((defectCon < absLimit || defectCon < initialDefectCon * reduction)
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
        Real maxDiffCon(0.0), maxDiffPot(0.0);
        calcMaxDiffConPot(maxDiffCon, maxDiffPot);

        //debug_verb << "maxDiffCon = " << maxDiffCon << std::endl;
        //debug_verb << "maxDiffPot = " << maxDiffPot << std::endl;

        // Use 'reduction' as error tolerance for absolute changes
        if(maxDiffCon < reduction && maxDiffPot < reduction)
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

    void calcDefect(Real& defectCon, Real& defectPot) const
    {
      U_CON rCon(gridoperatorCon.testGridFunctionSpace());
      rCon = 0.0;
      gridoperatorCon.residual(unewCon, rCon);
      // It doesn't matter which linear solver's norm we take, as the
      // current implementation of seqistlsolverbackend always take the 2-norm
      //Real defectCon = solverCon.getPDESolver().getLinearSolver().norm(rCon);
      defectCon = solverPot.getLinearSolver().norm(rCon);
      //debug_verb << "defectCon = " << defectCon << std::endl;

      U_POT rPot(gridoperatorPot.testGridFunctionSpace());
      rPot = 0.0;
      gridoperatorPot.residual(uPot, rPot);
      defectPot = solverPot.getLinearSolver().norm(rPot);
      //debug_verb << "defectPot = " << defectPot << std::endl;

      if (!std::isfinite(defectCon) || !std::isfinite(defectPot))
          DUNE_THROW(Dune::Exception,
                     "Acme0OperatorSplit::calcDefect(): Non-linear defect is NaN or Inf");

    }

    void calcMaxDiffConPot(Real& maxDiffCon, Real& maxDiffPot) const
    {
      // ================== Calculate max conc/pot change (works in 1D only!) ===========================
      typedef typename GV::template Codim<0>::Entity Entity;
      typedef typename GV::template Codim<0>::EntityPointer EntityPointer;
      typedef typename GV::template Codim<0>::Iterator ElementLeafIterator;

      typename DGF_CON::Traits::DomainType x_l;
      typename DGF_CON::Traits::DomainType x_r;

      typename DGF_CON::Traits::RangeType yoldCon_l;
      typename DGF_CON::Traits::RangeType yoldCon_r;
      typename DGF_CON::Traits::RangeType ynewCon_l;
      typename DGF_CON::Traits::RangeType ynewCon_r;

      typename DGF_POT::Traits::RangeType yoldPot_l;
      typename DGF_POT::Traits::RangeType yoldPot_r;
      typename DGF_POT::Traits::RangeType ynewPot_l;
      typename DGF_POT::Traits::RangeType ynewPot_r;

      maxDiffCon = 0.0;
      maxDiffPot = 0.0;

      for (ElementLeafIterator eit=gv.template begin<0>(); eit!=gv.template end<0>(); ++eit)
      {
        Entity& e = *eit;

        x_l = e.geometry().local(e.geometry().corner(0));
        x_r = e.geometry().local(e.geometry().corner(1));

        dgfConPrevious.evaluate(e, x_l, yoldCon_l);
        dgfConPrevious.evaluate(e, x_r, yoldCon_r);
        dgfCon.evaluate(e, x_l, ynewCon_l);
        dgfCon.evaluate(e, x_r, ynewCon_r);

        dgfPotPrevious.evaluate(e, x_l, yoldPot_l);
        dgfPotPrevious.evaluate(e, x_r, yoldPot_r);
        dgfPot.evaluate(e, x_l, ynewPot_l);
        dgfPot.evaluate(e, x_r, ynewPot_r);

        Real diffCon_l = 0.0, diffCon_r = 0.0;
        Real diffPot_l = 0.0, diffPot_r = 0.0;

        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          if(calcRelativeErrors)
          {
            // Check for division by 0
            if(std::abs(yoldCon_l[j]) > 1e-12)
            {
              diffCon_l = std::abs(ynewCon_l[j] / yoldCon_l[j] - 1.0);
            } else {
              // TODO Find a better way to handle this
              diffCon_l = std::abs(ynewCon_l[j] - yoldCon_l[j]);
            }
            if(std::abs(yoldCon_r[j]) > 1e-12)
            {
              diffCon_r = std::abs(ynewCon_r[j] / yoldCon_r[j] - 1.0);
            } else {
              // TODO Find a better way to handle this
              diffCon_r = std::abs(ynewCon_r[j] - yoldCon_r[j]);
            }
            if(std::abs(yoldPot_l[0]) > 1e-12)
            {
              diffPot_l = std::abs(ynewPot_l[0] / yoldPot_l[0] - 1.0);
            } else {
              // TODO Find a better way to handle this
              diffPot_l = std::abs(ynewPot_l[0] - yoldPot_l[0]);
            }
            if(std::abs(yoldPot_r[0]) > 1e-12)
            {
              diffPot_r = std::abs(ynewPot_r[0] / yoldPot_r[0] - 1.0);
            } else {
              // TODO Find a better way to handle this
              diffPot_r = std::abs(ynewPot_r[0] - yoldPot_r[0]);
            }
          } else { // Absolute errors
            diffCon_l = std::abs(ynewCon_l[j] - yoldCon_l[j]);
            diffCon_r = std::abs(ynewCon_r[j] - yoldCon_r[j]);
            diffPot_l = std::abs(ynewPot_l[0] - yoldPot_l[0]);
            diffPot_r = std::abs(ynewPot_r[0] - yoldPot_r[0]);
          }
          //debug_verb << e.geometry().global(xCon_l) << ": "
          //    << yoldCon_l[j] << " -> " << ynewCon_l[j] << " diffCon_l = " << diffCon_l << std::endl;
          //debug_verb << e.geometry().global(xCon_r) << ": "
          //    << ynewCon_r[j] << " -> " << ynewCon_r[j] << " diffCon_r = " << diffCon_r << std::endl;
          maxDiffCon = std::max(maxDiffCon, std::max(diffCon_l, diffCon_r));
          maxDiffPot = std::max(maxDiffPot, std::max(diffPot_l, diffPot_r));
        }
      }
      //debug_verb << "----------------------------------------" << std::endl;
      //debug_info << std::setprecision(4) << std::scientific << "  DiffPot = " << maxDiffPot;
      //debug_info << std::setprecision(4) << std::scientific << "  DiffCon = " << maxDiffCon << std::endl;
    }

    void setReduction(Real reduction_)
    {
      reduction = reduction_;
    }

    void setAbsLimit(Real absLimit_)
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

    U_CON& uoldCon;
    U_CON& unewCon;
    U_CON uConPrevious;   // copy of uoldCon, updated each iteration
    U_POT& uPot;
    U_POT uPotPrevious;   // copy of uPot, updated each iteration

    DIRICHLET_GF_CON& dirichletValuesCon;
    DIRICHLET_GF_POT& dirichletValuesPot;

    SOLVER_CON& solverCon;
    SOLVER_POT& solverPot;

    GFS_CON& gfsCon;
    GFS_POT& gfsPot;

    DGF_CON dgfCon;
    DGF_POT dgfPot;
    DGF_CON dgfConPrevious;
    DGF_POT dgfPotPrevious;

    PHYSICS& physics;
    GV& gv;
    GO_CON& gridoperatorCon;
    GO_POT& gridoperatorPot;

    Real reduction;
    Real absLimit;
    int maxIterations;
    bool useDefect;
    bool forceIteration;

    Real initialDefectCon;
    Real initialDefectPot;

    int iterations;

};


#endif /* DUNE_AX1_ACME0_OPERATOR_SPLIT_HH */
