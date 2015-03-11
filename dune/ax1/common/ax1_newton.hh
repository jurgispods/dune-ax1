/*
 * newton.hh
 *
 *  Created on: Aug 2, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_NEWTON_HH
#define DUNE_AX1_NEWTON_HH

#include <iostream>
#include <fstream>

#include <dune/common/ios_state.hh>

#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/instationary/onestep.hh>

#include <dune/ax1/common/ax1_solution_container.hh>
#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/rownorm_preconditioner.hh>

template<class AcmeOutput, class LOP, class GOS, class S, class TrlV, class TstV = TrlV>
class Ax1Newton : public Dune::PDELab::NewtonSolver<GOS,S,TrlV,TstV>,
                  public Dune::PDELab::NewtonTerminate<GOS,TrlV,TstV>,
                  public Dune::PDELab::NewtonLineSearch<GOS,TrlV,TstV>,
                  public Dune::PDELab::NewtonPrepareStep<GOS,TrlV,TstV>
{
  typedef GOS GridOperator;
  typedef S Solver;
  typedef TrlV TrialVector;
  typedef TstV TestVector;
  typedef typename GridOperator::Traits::TrialGridFunctionSpace::Traits::GridViewType GV;
  typedef typename TstV::ElementType RFType;
  //typedef typename GOS::template MatrixContainer<RFType>::Type Matrix;
  typedef typename GOS::Traits::Jacobian Matrix;

  typedef Dune::PDELab::NewtonLineSearch<GOS,TrlV,TstV> LineSearchBase;

  typedef Ax1SolutionContainer<TrialVector,typename AcmeOutput::Traits::GFS_POT_SUB,
      typename AcmeOutput::Traits::GFS_CON_SUB> SolutionContainer;

  public:
    enum LineSearchStrategy { noLineSearch,
                            hackbuschReusken,
                            hackbuschReuskenAcceptBest,
                            hackbuschReuskenAcceptBestNoThrow };

    Ax1Newton(AcmeOutput& acmeOutput_, LOP& lop_, GridOperator& go, TrialVector& u_, Solver& solver_, const char* basename_,
        SolutionContainer& solutionContainer_)
      : Dune::PDELab::NewtonBase<GOS,TrlV,TstV>(go,u_)
      //, Dune::PDELab::Newton<GOS,S,TrlV,TstV>(go, u_, solver_)
      , Dune::PDELab::NewtonSolver<GOS,S,TrlV,TstV>(go,u_,solver_)
      , Dune::PDELab::NewtonTerminate<GOS,TrlV,TstV>(go,u_)
      , Dune::PDELab::NewtonLineSearch<GOS,TrlV,TstV>(go,u_)
      , Dune::PDELab::NewtonPrepareStep<GOS,TrlV,TstV>(go,u_)
      , acmeOutput(acmeOutput_)
      , lop(lop_)
      , solutionContainer(solutionContainer_)
      , fn_matrix(std::string(basename_) + "_mat")
      , fn_vector(std::string(basename_) + "_rhs")
      , fn_vector2(std::string(basename_) + "_sol")
      , row_preconditioner(false)
      , printMatrix(false)
      , printRhs(false)
      , reorderMatrix(false)
      , printResidual(false)
      , nDOFsCytosol(0)
      , nDOFsMembrane(0)
      , lastReduction(0.0)
      , residual(this->gridoperator.testGridFunctionSpace())
      , subGfsPot(this->gridoperator.testGridFunctionSpace())
      , subGfsCon(this->gridoperator.testGridFunctionSpace())
      , dgfResidualPot(subGfsPot, residual)
      , dgfResidualCon(subGfsCon, residual)
      , fullyImplicit(false)
    {
      init();
    }

    Ax1Newton(AcmeOutput& acmeOutput_, LOP& lop_, GridOperator& go, Solver& solver_, const char* basename_,
        SolutionContainer& solutionContainer_)
      : Dune::PDELab::NewtonBase<GOS,TrlV,TstV>(go)
      //, Dune::PDELab::Newton<GOS,S,TrlV,TstV>(go, solver_)
      , Dune::PDELab::NewtonSolver<GOS,S,TrlV,TstV>(go,solver_)
      , Dune::PDELab::NewtonTerminate<GOS,TrlV,TstV>(go)
      , Dune::PDELab::NewtonLineSearch<GOS,TrlV,TstV>(go)
      , Dune::PDELab::NewtonPrepareStep<GOS,TrlV,TstV>(go)
      , acmeOutput(acmeOutput_)
      , lop(lop_)
      , solutionContainer(solutionContainer_)
      , fn_matrix(std::string(basename_) + "_mat")
      , fn_vector(std::string(basename_) + "_rhs")
      , fn_vector2(std::string(basename_) + "_sol")
      , row_preconditioner(false)
      , printMatrix(false)
      , printRhs(false)
      , reorderMatrix(false)
      , printResidual(false)
      , nDOFsCytosol(0)
      , nDOFsMembrane(0)
      , lastReduction(0.0)
      , residual(this->gridoperator.testGridFunctionSpace())
      , subGfsPot(this->gridoperator.testGridFunctionSpace())
      , subGfsCon(this->gridoperator.testGridFunctionSpace())
      , dgfResidualPot(subGfsPot, residual)
      , dgfResidualCon(subGfsCon, residual)
      , fullyImplicit(false)
    {
      init();
    }

    //! Add some debug data to the output class' diagInfo structure
    void init()
    {
      // Register debug data from Newton result
      // Most of these values have already been communicated; others like timings are subjective
      // anway; so don't to any communcation of these values, each processor keeps its own values!
      acmeOutput.getDiagInfo().registerDebugData("first_defect",0.0);
      acmeOutput.getDiagInfo().registerDebugData("defect",0.0);
      acmeOutput.getDiagInfo().registerDebugData("reduction",0.0);
      acmeOutput.getDiagInfo().registerDebugData("conv_rate",0.0);
      acmeOutput.getDiagInfo().registerDebugData("t_assembler",0.0);
      acmeOutput.getDiagInfo().registerDebugData("t_lin_solv",0.0);
      acmeOutput.getDiagInfo().registerDebugData("t_elapsed",0.0);
      acmeOutput.getDiagInfo().registerDebugData("avg_lin_it",0.0);

      std::stringstream header;
      header << "# time: -1" << std::endl
             << "0.000000000000e+00 0.000000000000e+00 -0.000000000000e+00"
             << std::endl << std::endl;
      acmeOutput.initGnuplotFile("residual_pot_init", header.str());
      acmeOutput.initGnuplotFile("residual_con_init", header.str());
      acmeOutput.initGnuplotFile("residual_pot", header.str());
      acmeOutput.initGnuplotFile("residual_con", header.str());

      //      debug_jochen << "MultiGFS: " << std::endl;
      //      debug_jochen << Tools::getTypeName(this->gridoperator.testGridFunctionSpace()) << std::endl;
      //      debug_jochen << Tools::getTypeName(this->gridoperator.testGridFunctionSpace().gridView()) << std::endl;
      //      debug_jochen << "Residual GF: " << std::endl;
      //      debug_jochen << Tools::getTypeName(dgfResidualPot.getGridView()) << std::endl;
    }

    void setRowPreconditioner(bool row_preconditioner_)
    {
     row_preconditioner = row_preconditioner_;
    }

    void setPrintMatrix(bool printMatrix_)
    {
        printMatrix = printMatrix_;
    }

    void setPrintRhs(bool printRhs_)
    {
        printRhs = printRhs_;
    }

    void setReorderMatrix(bool reorder)
    {
      reorderMatrix = reorder;
    }

    void setPrintResidual(bool printResidual_)
    {
      printResidual = printResidual_;
    }

    void setFullyImplicit(bool fullyImplicit_)
    {
      fullyImplicit = fullyImplicit_;
    }

    std::vector<int>& getPermutation()
    {
      return permute;
    }

    void setPermutation(std::vector<int>& permute_)
    {
      permute = permute_;
    }

    std::vector<int>& getInversePermutation()
    {
      return inv_permute;
    }

    void setInversePermutation(std::vector<int>& inv_permute_)
    {
      inv_permute = inv_permute_;
    }

    virtual void defect(TestVector& r)
    {
      //debug_jochen << "Ax1Newton::defect()" << std::endl;
      Dune::PDELab::NewtonSolver<GOS,S,TrlV,TstV>::defect(r);

      debug_verb << "  [Ax1Newton::defect] calculated defect is " << this->res.defect << std::endl;

      //Output::printRawCoefficientVector(r, "r after defect()");

      // Store residual for later output
      residual = r;
    }


    /*!
     * \brief This method overwrites the corresponding one in Dune::PDELab::Newton
     * by adding a row preconditioner to the matrix and optionally writing the
     * resulting matrix to a file for later inspection/evaluation.
     *
     * \param A The system matrix
     * \param r The residual
     */
    virtual void prepare_step(Matrix& A, TestVector& r)
    {
      // Write out initial residual
      if(printResidual && this->res.iterations == 0)
      {
        //debug_jochen << "writing out initial residual" << std::endl;
        //Output::printSingleCoefficientVector(residual, "residual");
        assert(this->res.defect == this->res.first_defect);

        std::stringstream info;
        info << "   total defect=" << this->res.defect;
        acmeOutput.writeGFToGnuplot(dgfResidualPot, "residual_pot_init", info.str());
        acmeOutput.writeGFToGnuplot(dgfResidualCon, "residual_con_init", info.str());
      }

      // Restore original sparsity pattern which was possibly destroyed by reordering in the last iteration
      if(reorderMatrix)
      {
        Matrix A_init(this->gridoperator);
        A = A_init;
      }

      if(fullyImplicit)
      {
        // Inform solutionContainer about the current Newton solution
        // before doing any assembly in local operator!
        solutionContainer.setSolutionConNew(this->u);
        solutionContainer.setSolutionPotNew(this->u);

        debug_info << "Calling lop.preAssembly() for Newton iteration #" << this->res.iterations << "..." << std::endl;
        // Use the version without handing over a solution vector; prefer the helper GFs to extract the relevant solution
        // vectors from the solutionContainer
        lop.preAssembly();
      }
      Dune::PDELab::NewtonPrepareStep<GOS,TrlV,TstV>::prepare_step(A, r);

      if(this->reassembled)
      {
        // Add row-norm preconditioner
        if (row_preconditioner)
        {
          RowNormPreconditioner<GV,AX1_BLOCKSIZE> preconditioner(Dune::PDELab::istl::raw(A),Dune::PDELab::istl::raw(r));
          //DUNE_THROW(Dune::NotImplemented, "Mach et, Junge!");
        }

        // Now: Reorder matrix!
        if(reorderMatrix)
        {
          debug_info << "  Custom matrix reordering in progress..." << std::endl;
          Dune::ios_base_all_saver nuss(std::cout);
          Dune::Timer timer_reordering;
          this->reorderMatrixAndRHS(A,r);
          debug_info << "  Finished custom matrix reordering, time:"
              << std::setw(12) << std::setprecision(4) << std::scientific
              << timer_reordering.elapsed() << std::endl;
        }
      }

      // Write matrix to file ((in first iteration))
      if (/*this->verbosity_level >= 4*/ printMatrix /*&& this->res.iterations == 0*/)
      {
        std::string filename = getMatFileName();
        debug_verb << "Writing matrix to file [" << filename << "]" << std::endl;
        Dune::writeMatrixToMatlab(A.base(), filename);
        if(printRhs)
        {
          std::string filename = getRhsFileName();
          std::ofstream filestream(filename);
          Dune::printvector(filestream, r.base(), std::string("%rhs"), std::string(""), 1, 20, 12);
        }
      }

//      debug_verb << "==== ITERATION #" << this->res.iterations << "===================" << std::endl;
//      debug_verb << "u_" << this->res.iterations << std::endl;
//      Output::printMultipleComponentCoefficientVectorDG(*this->u, NUMBER_OF_SPECIES+1);
//      debug_verb << "r_" << this->res.iterations << std::endl;
//      Output::printMultipleComponentCoefficientVectorDG(r, NUMBER_OF_SPECIES+1);
    }

    void setLineSearchStrategy(LineSearchStrategy strategy_)
    {
      strategy = strategy_;
    }


    /**!
     * This is the only method of any father class where the member 'u' (pointer to
     * full coefficient vector) is modified.
     *
     * @param z
     * @param r
     */
    virtual void line_search(TrialVector& z, TestVector& r)
    {
      //Output::printRawCoefficientVector(z, "z after linear_solve()");

//      if(printMatrix && printRhs)
//      {
//        std::string filename = getLinearSolutionFileName();
//        filename += ".unordered";
//        std::ofstream filestream(filename);
//        Dune::printvector(filestream, z.base(), std::string("%sol_unordered"), std::string(""), 1, 20, 12);
//      }

      if(reorderMatrix)
      {
        debug_jochen << "[Ax1Newton::line_search()] Reordering solution..." << std::endl;
        TestVector z_new(z);
        for(int i=0; i<z.flatsize(); i++)
        {
          // i is the permuted index, get back the original one
          int i_orig = permute[i];

          Dune::PDELab::istl::raw(z_new)[i_orig] = Dune::PDELab::istl::raw(z)[i];
        }
        z = z_new;
        debug_jochen << "[Ax1Newton::line_search()] Done." << std::endl;
      }

      //Output::printSingleCoefficientVector(z, "z_new");

      // TODO For reordered matrix, better print out the reordered solution (see above)!
      if(printMatrix && printRhs)
      {
        std::string filename = getLinearSolutionFileName();
        std::ofstream filestream(filename);
        Dune::printvector(filestream, z.base(), std::string("%sol"), std::string(""), 1, 20, 12);
      }

      // ============ Following: Modified copy of PDELab Newton's line_search ===================================
      /*
       * Note: My modifications with respect to the original Newton line search contain updating the solution
       * container each time (*this->u) is changed in the line seach process. Additionally, the preAssembly()
       * routine is called which updates membrane states in order to yield a consistent, fully-implicit numerical
       * scheme.
       *
       * One further addition is the new line search strategy 'hackbuschReuskenAcceptBestNoThrow', which basically
       * works in the same way as 'hackbuschReuskenAcceptBest', but instead of throwing an exception when it does
       * not converge, it simply accepts the full Newton step anyway. This is because in some cases, a Newton
       * iteration resutls in an _increase_ of the defect rather than a reduction. In subsequent iterations however,
       * the defect is reduced as desired. In such a case, the line search algorithm can not converge, so we simply
       * accept that fact and hope that it will do so again in the following iteration.
       */
      if (strategy == LineSearchStrategy::noLineSearch)
      {
        // u = u - z;
        this->u->axpy(-1.0, z);                     // TODO: vector interface
        solutionContainer.setSolutionConNew(this->u);
        solutionContainer.setSolutionPotNew(this->u);
        lop.preAssembly();
        this->defect(r);
        return;
      }

      if (this->verbosity_level >= 4)
        std::cout << "      Performing line search..." << std::endl;
      RFType lambda = 1.0;
      RFType best_lambda = 0.0;
      RFType best_defect = this->res.defect;

      if (this->verbosity_level >= 4)
        std::cout << "      Initial defect: " << best_defect << std::endl;

      TrialVector prev_u(*this->u);  // TODO: vector interface
      unsigned int i = 0;
      Dune::ios_base_all_saver restorer(std::cout); // store old ios flags

      while (1)
      {
        if (this->verbosity_level >= 4)
        {
          std::cout << "          trying line search damping factor:   "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << lambda
                    << std::endl;
        }

        this->u->axpy(-lambda, z);                  // TODO: vector interface
        // Inform solutionContainer about the current Newton solution
        solutionContainer.setSolutionConNew(this->u);
        solutionContainer.setSolutionPotNew(this->u);
        lop.preAssembly();
        try {
          this->defect(r);
        }
        catch (Dune::PDELab::NewtonDefectError)
        {
          if (this->verbosity_level >= 4)
            std::cout << "          Nans detected" << std::endl;
        }       // ignore NaNs and try again with lower lambda

        if (this->res.defect <= (1.0 - lambda/4) * this->prev_defect)
        {
          if (this->verbosity_level >= 4)
            std::cout << "          line search converged" << std::endl;
          break;
        }

        if (this->res.defect < best_defect)
        {
          best_defect = this->res.defect;
          best_lambda = lambda;
        }

        if (++i >= LineSearchBase::maxit)
        {
          if (this->verbosity_level >= 4)
            std::cout << "          max line search iterations exceeded" << std::endl;

          switch (strategy)
          {
            case LineSearchStrategy::hackbuschReusken:
              *this->u = prev_u;
              solutionContainer.setSolutionConNew(this->u);
              solutionContainer.setSolutionPotNew(this->u);
              lop.preAssembly();
              this->defect(r);
              DUNE_THROW(Dune::PDELab::NewtonLineSearchError,
                         "NewtonLineSearch::line_search(): line search failed, "
                         "max iteration count reached, "
                         "defect did not improve enough");
            case LineSearchStrategy::hackbuschReuskenAcceptBest:
              if (best_lambda == 0.0)
              {
                *this->u = prev_u;
                solutionContainer.setSolutionConNew(this->u);
                solutionContainer.setSolutionPotNew(this->u);
                lop.preAssembly();
                this->defect(r);
                DUNE_THROW(Dune::PDELab::NewtonLineSearchError,
                           "NewtonLineSearch::line_search(): line search failed, "
                           "max iteration count reached, "
                           "defect did not improve in any of the iterations");
              }
              if (best_lambda != lambda)
              {
                *this->u = prev_u;
                this->u->axpy(-best_lambda, z);
                // Inform solutionContainer about the current Newton solution
                solutionContainer.setSolutionConNew(this->u);
                solutionContainer.setSolutionPotNew(this->u);
                lop.preAssembly();
                this->defect(r);
              }
              break;
            case LineSearchStrategy::hackbuschReuskenAcceptBestNoThrow:
              if (best_lambda == 0.0)
              {
                // Experimental handling of non-converging line search: Accept full step anyway,
                // even if the defect went up instead of down; hope that a subsequent Newton iteration
                // will again result in a defect reduction
                best_lambda = 1.0;
                *this->u = prev_u;
                this->u->axpy(-best_lambda, z);
                solutionContainer.setSolutionConNew(this->u);
                solutionContainer.setSolutionPotNew(this->u);
                lop.preAssembly();
                this->defect(r);
                debug_warn << "NewtonLineSearch::line_search(): line search failed, "
                  << "max iteration count reached, defect did not improve in any of the iterations. "
                  << "Accepting full Newton step anyway!" << std::endl;
                // No throw here!
              }
              if (best_lambda != lambda)
              {
                *this->u = prev_u;
                this->u->axpy(-best_lambda, z);
                // Inform solutionContainer about the current Newton solution
                solutionContainer.setSolutionConNew(this->u);
                solutionContainer.setSolutionPotNew(this->u);
                lop.preAssembly();
                this->defect(r);
              }
              break;
            case LineSearchStrategy::noLineSearch:
              break;
          }
          break;
        }

        lambda *= LineSearchBase::damping_factor;
        *this->u = prev_u;                          // TODO: vector interface
      }
      if (this->verbosity_level >= 4)
        std::cout << "          line search damping factor:   "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << lambda << std::endl;
    }

//    virtual void defect(TestVector& r)
//    {
//      Dune::PDELab::NewtonSolver<GOS,S,TrlV,TstV>::defect(r);
//      debug_verb << "Defect vor : "<< this->res.defect << std::endl;
//
//      Tools::compositeToChildCoefficientVector(this->gridoperator.trialGridFunctionSpace(), *this->u, uCon, 0);
//      Tools::compositeToChildCoefficientVector(this->gridoperator.trialGridFunctionSpace(), *this->u, uPot, 1);
//      Dune::PDELab::NewtonSolver<GOS,S,TrlV,TstV>::defect(r);
//      debug_verb << "Defect nach: "<< this->res.defect << std::endl;
//    }


//    virtual bool terminate()
//    {
//      bool terminated = Dune::PDELab::NewtonTerminate<GOS,TrlV,TstV>::terminate();
//
//      double thisReduction = this->res.defect/this->prev_defect;
//
//      debug_verb << "lastReduction = " << lastReduction << std::endl;
//      debug_verb << "thisReduction = " << thisReduction << std::endl;
//
//      if(!terminated && std::abs(lastReduction-1.0)<1e-6 && std::abs(thisReduction-1.0)<1e-6)
//      {
//        debug_warn << "======================" << std::endl;
//        debug_warn << "Newton reduction was 1.0 is two succeeding iterations; assuming minimum possible defect was reached! Terminating."
//            << std::endl;
//        debug_warn << "======================" << std::endl;
//
//        terminated = true;
//      }
//      lastReduction = thisReduction;
//
//      return terminated;
//    }

    //! Overwrite apply(TrialVector) method from Dune::PDELab::NewtonSolver
    void apply(TrialVector& u_)
    {
      this->u = &u_;
      apply();
    }

    //! Overwrite apply() method from Dune::PDELab::NewtonSolver
    void apply()
    {
      try
      {
        lastReduction = 0.0;

        // Call base class apply()
        Dune::PDELab::NewtonSolver<GOS,S,TrlV,TstV>::apply();

      } catch (Dune::PDELab::NewtonError& e)
      {
        // Fill diagnostic info for debug purposes, even if Newton did not converge
        acmeOutput.getDiagInfo().setDebugData("first_defect",this->res.first_defect);
        acmeOutput.getDiagInfo().setDebugData("defect",this->res.defect);
        acmeOutput.getDiagInfo().setDebugData("reduction",this->res.reduction);
        acmeOutput.getDiagInfo().setDebugData("conv_rate",this->res.conv_rate);
        acmeOutput.getDiagInfo().setDebugData("t_assembler",this->res.assembler_time);
        acmeOutput.getDiagInfo().setDebugData("t_lin_solv",this->res.linear_solver_time);
        acmeOutput.getDiagInfo().setDebugData("t_elapsed",this->res.elapsed);
        acmeOutput.getDiagInfo().setDebugData("avg_lin_it", (this->res.iterations > 0 ?
            (this->res.linear_solver_iterations / this->res.iterations) : this->res.linear_solver_iterations));

        // Write out final residual
        if(printResidual)
        {
          std::stringstream info;
          info << "   total defect=" << this->res.defect;
          acmeOutput.writeGFToGnuplot(dgfResidualPot, "residual_pot", info.str());
          acmeOutput.writeGFToGnuplot(dgfResidualCon, "residual_con", info.str());
        }
        throw e;
      }

      double avgLinearIterations = this->res.linear_solver_iterations;
      double avgLinearSolverTime = this->res.linear_solver_time;
      if(this->res.iterations > 0)
      {
        avgLinearIterations /= this->res.iterations;
        avgLinearSolverTime /= this->res.iterations;
      }

      bool converged = this->res.converged;
      if(converged && this->verbosity_level > 2)
      {
        debug_info << "  Total number of linear iterations: " << this->res.linear_solver_iterations
                  << std::fixed << " (average " << avgLinearIterations << ")"
                  << std::endl;

        debug_info << "  Total linear solver time: " << this->res.linear_solver_time
                  << std::fixed << " (average " << avgLinearSolverTime << ")"
                  << std::endl;
      }
      acmeOutput.getDiagInfo().setDebugData("first_defect",this->res.first_defect);
      acmeOutput.getDiagInfo().setDebugData("defect",this->res.defect);
      acmeOutput.getDiagInfo().setDebugData("reduction",this->res.reduction);
      acmeOutput.getDiagInfo().setDebugData("conv_rate",this->res.conv_rate);
      acmeOutput.getDiagInfo().setDebugData("t_assembler",this->res.assembler_time);
      acmeOutput.getDiagInfo().setDebugData("t_lin_solv",this->res.linear_solver_time);
      acmeOutput.getDiagInfo().setDebugData("t_elapsed",this->res.elapsed);
      acmeOutput.getDiagInfo().setDebugData("avg_lin_it",avgLinearIterations);

      // Write out final residual
      if(printResidual)
      {
        std::stringstream info;
        info << "   total defect=" << this->res.defect;
        acmeOutput.writeGFToGnuplot(dgfResidualPot, "residual_pot", info.str());
        acmeOutput.writeGFToGnuplot(dgfResidualCon, "residual_con", info.str());
      }
    }

  protected:

    /**
     * This method attempts to reorder a matrix and rhs when it has not the desired vertex-block
     * structure (meaning the concentrations are vertex-blocked, but the potential unknowns come
     * at the end of the unknown vector as a separate block)
     *
     * This is a crude hack!
     *
     * At the end, A2 will contain the reordered matrix with blocksize 1, and A3 will contain
     * the reordered matrix with blocksize 4 - same holds for r2 and r3.
     *
     * As the whole GFS setup would have to be changed to able to use A3/r3, A2/r2 are assigned
     * to A/r, repectively.
     */
    void reorderMatrixAndRHS(Matrix& A, TestVector& r)
    {
      // ================= Reorder matrix ========================== //

      // Some typedefs
      typedef typename GridOperator::Traits::TrialGridFunctionSpace MultiGFS;
      typedef typename MultiGFS::Traits::GridViewType MultiGV;
      typedef typename MultiGV::Traits::IndexSet IndexSet;

      int n_total = this->gridoperator.trialGridFunctionSpace().size();
      int n_con = this->gridoperator.trialGridFunctionSpace().template child<0>().size();
      int n_pot = this->gridoperator.trialGridFunctionSpace().template child<1>().size();

      assert(n_con + n_pot == n_total);
      assert(permute.size() == n_total && inv_permute.size() == n_total);

      if(n_total != A.N())
      {
        DUNE_THROW(Dune::Exception, "Blocksize>1 detected, reordering matrix does only work with blocksize==1! (multigfs.size()="
            << n_total << ", A.N()=" << A.N());
      }

      // iterator types
      typedef typename Dune::PDELab::istl::raw_type<Matrix>::type RawMatrix;
      typedef typename RawMatrix::RowIterator rowiterator;
      typedef typename RawMatrix::ColIterator coliterator;
      typedef typename RawMatrix::row_type row_type;
      typedef typename RawMatrix::block_type block;

      RawMatrix& A_raw = Dune::PDELab::istl::raw(A);

      RawMatrix A2(A.N(), A.N(), A_raw.nonzeroes(), RawMatrix::row_wise);
      //RawMatrix P(A.N(),A.N(),A.N(),RawMatrix::row_wise);
      //RawMatrix P_inv(A.N(),A.N(),A.N(),RawMatrix::row_wise);

      typedef typename Dune::PDELab::istl::raw_type<TestVector>::type RawTestVector;
      RawTestVector& r_raw(r);
      RawTestVector r2(r_raw);

      typedef typename RawMatrix::CreateIterator Iter;

      // Build up sparsity pattern
      for (Iter row = A2.createbegin(); row != A2.createend(); ++row)
      {
        // Index of permuted matrix
        int i = row.index();

        // Get corresponding row of the original matrix
        int i_orig = permute[i];
        row_type orig_row = A_raw[i_orig];

        //debug_jochen << "New row " << i << ": ";

        // Loop over column entries within this row
        for (coliterator orig_col = orig_row.begin(); orig_col != orig_row.end(); ++orig_col)
        {
          int j_orig = orig_col.index();
          int j = inv_permute[j_orig];

          //debug_jochen << j << " ";

          // Initialize entry
          row.insert(j);
        }
        //debug_jochen << std::endl;
      }

      // Insert values
      for (rowiterator row = A2.begin(); row != A2.end(); ++row)
      {
        int i = row.index();
        int i_orig = permute[i];
        for (coliterator col = (*row).begin(); col != (*row).end(); ++col)
        {
          int j = col.index();
          int j_orig = permute[j];

          //debug_jochen << "Trying to assign"
          //  << "A2[" << i << "][" << j << "] with A[" << i_orig << "][" << j_orig << "] = "
          //  << A_raw[i_orig][j_orig] << std::endl;

          // Assign entry
          A2[i][j] = A_raw[i_orig][j_orig];
        }
      }

      // Also permute rhs!
      for (int i = 0; i < A_raw.N(); ++i)
      {
        int i_orig = permute[i];
        r2[i] = r_raw[i_orig];
      }

      // ================= END reorder matrix ========================== //


      // ============= Created blocked matrix from reordered matrix ============== //
      // This is only possible for one membrane element setups; since the ordering functionality for this case is now built-in,
      // the following blocking does not make sence and is disabled. It is furthermore not possible to have a Newton templated
      // with a blocksize-1 vector and change this vector to a blockize-4 vector in this function, so this would actually never
      // really work. It is jsut a proof-of-concept how to add blocking to a non-blocked matrix/vector.
//      const int BLOCKSIZE = 4;
//
//      typedef Dune::BlockVector<Dune::FieldVector<double, BLOCKSIZE> > BlockVector;
//      BlockVector r3(r2.flatsize() / BLOCKSIZE, r2.flatsize() / BLOCKSIZE);
//
//      for (int i = 0; i < r2.flatsize(); i++)
//      {
//        r3[i / BLOCKSIZE][i % BLOCKSIZE] = r2[i];
//      }
//
//      // Copy blocksize-1 matrix to blocksize-4 matrix
//      typedef Dune::BCRSMatrix<Dune::FieldMatrix<double, BLOCKSIZE, BLOCKSIZE> > BlockMatrix;
//      // iterator types
//      typedef typename BlockMatrix::RowIterator block_rowiterator;
//      typedef typename BlockMatrix::ColIterator block_coliterator;
//      typedef typename BlockMatrix::row_type block_row_type;
//      typedef typename BlockMatrix::block_type block_block;
//      typedef typename BlockMatrix::size_type block_size_type;
//
//      debug_jochen << "A.nonzeroes: " << A.nonzeroes() << std::endl;
//
//      block_size_type N_BLOCK = A.N() / BLOCKSIZE;
//      block_size_type NNZ_BLOCK = A.nonzeroes(); // viel zu gross!
//
//      BlockMatrix A3(N_BLOCK, N_BLOCK, NNZ_BLOCK, BlockMatrix::row_wise);
//
//      typedef BlockMatrix::CreateIterator BIter;
//      // Build up sparsity pattern
//      for (BIter row = A3.createbegin(); row != A3.createend(); ++row)
//      {
//        // Index of permuted matrix
//        int i = row.index();
//
//        // Loop over rows of original matrix and check if nonzeroes exist in this block
//        for (int i_orig = BLOCKSIZE * i; i_orig < BLOCKSIZE * (i + 1); i_orig++)
//        {
//          row_type orig_row = A2[i_orig];
//          // Loop over column entries within this row
//          for (coliterator orig_col = orig_row.begin(); orig_col != orig_row.end(); ++orig_col)
//          {
//            int j_orig = orig_col.index();
//            int j = j_orig / BLOCKSIZE;
//
//            assert(A2.exists(i_orig,j_orig));
//            //debug_jochen << j << " ";
//
//            // Initialize entry
//            row.insert(j);
//          }
//        }
//      }
//
//      // Insert values
//      for (block_rowiterator row = A3.begin(); row != A3.end(); ++row)
//      {
//        int i = row.index();
//        for (block_coliterator col = (*row).begin(); col != (*row).end(); ++col)
//        {
//          int j = col.index();
//
//          //debug_jochen << "Block A3[" << i << "][" << j << "]" << std::endl;
//
//          // Loop over values of original matrix
//          for (int i_orig = BLOCKSIZE * i; i_orig < BLOCKSIZE * (i + 1); i_orig++)
//          {
//            for (int j_orig = BLOCKSIZE * j; j_orig < BLOCKSIZE * (j + 1); j_orig++)
//            {
//              // Check if value exists
//              if (A2.exists(i_orig, j_orig))
//              {
//                //debug_jochen << "  assign ["
//                //  << (i_orig % BLOCKSIZE) << "][" << (j_orig % BLOCKSIZE)
//                //  << "] with A2[" << i_orig << "][" << j_orig << "] = "
//                //  << A2[i_orig][j_orig] << std::endl;
//
//                // Fill nested matrix with original value
//                A3[i][j][i_orig % BLOCKSIZE][j_orig % BLOCKSIZE] = A2[i_orig][j_orig];
//
//              } else
//              {
//                debug_jochen << "Value A2[" << i_orig << "][" << j_orig << "] does not exist!"
//                    << std::endl;
//              }
//            }
//          }
//        }
//      }
//
//
//      // Better print this shit
//      if (printMatrix)
//      {
//        std::stringstream kack;
//        kack << "octave/block_mat_" << (this->res.iterations + 1) << ".dat";
//        debug_verb << "Writing matrix to file [" << kack.str() << "]" << std::endl;
//        Dune::writeMatrixToMatlab(A2, kack.str());
//        if (printRhs)
//        {
//          std::stringstream kack2;
//          kack2 << "octave/block_rhs_" << (this->res.iterations + 1) << ".dat";
//          std::ofstream filestream(kack2.str());
//          Dune::printvector(filestream, r2, std::string("%rhs"), std::string(""));
//        }
//      }
//
//
//      // This is an example of how the blocked, reordered system could be solved!
//      BlockVector z3(r3);
//      z3 = 0;
//
//      Dune::MatrixAdapter<BlockMatrix, BlockVector, BlockVector> opa(A3);
//      Dune::SeqILU0<BlockMatrix, BlockVector, BlockVector> ilu0(A3, 1.0);
//
//      Dune::BiCGSTABSolver<BlockVector> my_ls(opa, ilu0, 1e-5, 5000, 5);
//      Dune::InverseOperatorResult stat;
//      my_ls.apply(z3, r3, stat);
//
//      Dune::PDELab::LinearSolverResult<double> res;
//      res.converged = stat.converged;
//      res.iterations = stat.iterations;
//      res.elapsed = stat.elapsed;
//      res.reduction = stat.reduction;
//      res.conv_rate = stat.conv_rate;
//
//      std::stringstream kack3;
//      kack3 << "octave/block_sol_" << (this->res.iterations + 1) << ".dat";
//      std::ofstream filestream(kack3.str());
//      Dune::printvector(filestream, z3, std::string("%sol"), std::string(""));
      // ============= END created blocked matrix from reordered matrix ============== //

      // Use reordered matrix and rhs
      Dune::PDELab::istl::raw(A) = A2;
      Dune::PDELab::istl::raw(r) = r2;
    }


    void getFileName(std::string& filename)
    {
      int mpi_rank = this->gridoperator.trialGridFunctionSpace().gridView().comm().rank();
      std::stringstream temp;
      temp << filename;
      temp << "_" << std::setw(2) << std::setfill('0') << (this->res.iterations + 1);
      temp << "_rank" << std::setw(2) << std::setfill('0') << mpi_rank;
      temp << ".dat";
      filename = temp.str();
    }

    std::string getMatFileName()
    {
      if(this->res.iterations == 0)
      {
        fn_matrix.increment();
      }
      std::string filename = fn_matrix.getName();
      getFileName(filename);
      return filename;
    }

    std::string getRhsFileName()
    {
      if(this->res.iterations == 0)
      {
        fn_vector.increment();
      }
      std::string filename = fn_vector.getName();
      getFileName(filename);
      return filename;
    }

    std::string getLinearSolutionFileName()
    {
      // This is called after the first solve
      if(this->res.iterations == 0)
      {
        fn_vector2.increment();
      }
      std::string filename = fn_vector2.getName();
      getFileName(filename);
      return filename;
    }



  private:
    AcmeOutput& acmeOutput;
    LOP& lop;
    SolutionContainer& solutionContainer;
    Dune::PDELab::FilenameHelper fn_matrix;
    Dune::PDELab::FilenameHelper fn_vector;
    Dune::PDELab::FilenameHelper fn_vector2;
    bool row_preconditioner;
    bool printMatrix;
    bool printRhs;
    bool printResidual;
    int nDOFsCytosol;
    int nDOFsMembrane;
    double lastReduction;

    bool reorderMatrix;
    std::vector<int> permute;
    std::vector<int> inv_permute;

    TstV residual;
    typename AcmeOutput::GFS_POT subGfsPot;
    typename AcmeOutput::GFS_CON subGfsCon;
    typename AcmeOutput::Traits::DGF_POT dgfResidualPot;
    typename AcmeOutput::Traits::DGF_CON dgfResidualCon;

    bool fullyImplicit;

    LineSearchStrategy strategy;


};

#endif /* DUNE_AX1_NEWTON_HH */
