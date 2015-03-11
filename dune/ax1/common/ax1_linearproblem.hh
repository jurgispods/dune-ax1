/*
 * ax1_linearproblem.hh
 *
 *  Created on: Aug 24, 2011
 *      Author: jpods
 *      (stolen from <dune/pdelab/stationary/linearproblem.hh>)
 */

#ifndef AX1_LINEARPROBLEM_HH
#define AX1_LINEARPROBLEM_HH

#include<dune/common/timer.hh>
#include<dune/pdelab/backend/backendselector.hh>

#include<dune/ax1/common/constants.hh>
#include<dune/ax1/common/rownorm_preconditioner.hh>

//===============================================================
// A class for solving linear stationary problems.
// It assembles the matrix, computes the right hand side and
// solves the problem.
// This is only a first vanilla implementation which has to be improved.
//===============================================================

template<class GOS, class LS, class V>
class Ax1StationaryLinearProblemSolver
{
  typedef typename V::ElementType Real;
  typedef typename GOS::template MatrixContainer<Real>::Type M;
  typedef typename GOS::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
  typedef typename Dune::PDELab::BackendVectorSelector<TrialGridFunctionSpace,Real>::Type W;
  typedef typename GOS::Traits::TrialGridFunctionSpace::Traits::GridViewType GV;

  public:

    Ax1StationaryLinearProblemSolver (const GOS& gos_, V& x_, LS& ls_, typename V::ElementType reduction_,
        const char* basename_, typename V::ElementType mindefect_ = 1e-99)
      : gos(gos_), ls(ls_), x(&x_), reduction(reduction_), mindefect(mindefect_)
      , fn_matrix(std::string(basename_) + "_mat")
      , fn_vector(std::string(basename_) + "_rhs")
      , printMatrix(false)
      , printRhs(false)
      , row_preconditioner(false)
      , verbosity(0)
    {
    }

    Ax1StationaryLinearProblemSolver (const GOS& gos_, LS& ls_, V& x_, typename V::ElementType reduction_,
        const char* basename_, typename V::ElementType mindefect_ = 1e-99)
      : gos(gos_), ls(ls_), x(&x_), reduction(reduction_), mindefect(mindefect_)
      , fn_matrix(std::string(basename_) + "_mat")
      , fn_vector(std::string(basename_) + "_rhs")
      , printMatrix(false)
      , printRhs(false)
      , row_preconditioner(false)
      , verbosity(0)
    {
    }

    Ax1StationaryLinearProblemSolver (const GOS& gos_, LS& ls_, typename V::ElementType reduction_,
        const char* basename_, typename V::ElementType mindefect_ = 1e-99)
      : gos(gos_), ls(ls_), x(0), reduction(reduction_), mindefect(mindefect_)
      , fn_matrix(std::string(basename_) + "_mat")
      , fn_vector(std::string(basename_) + "_rhs")
      , printMatrix(false)
      , printRhs(false)
      , row_preconditioner(false)
      , verbosity(0)
    {
    }

    void apply (V& x_) {
      x = &x_;
      apply();
    }

    void apply ()
    {
      Dune::Timer watch;
      double timing;

      // assemble matrix; optional: assemble only on demand!
      watch.reset();

      M m(gos);

      timing = watch.elapsed();
      // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
      if (gos.trialGridFunctionSpace().gridView().comm().rank()==0)
      {
        if(verbosity > 1) {
          debug_verb << "=== matrix setup (max) " << timing << " s" << std::endl;
        }
      }
      watch.reset();

      m = 0.0;
      gos.jacobian(*x,m);

      timing = watch.elapsed();
      // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
      if (gos.trialGridFunctionSpace().gridView().comm().rank()==0)
      {
        if(verbosity > 1) {
          debug_verb << "=== matrix assembly (max) " << timing << " s" << std::endl;
        }
      }

      // assemble residual
      watch.reset();

      W r(gos.testGridFunctionSpace(),0.0);
      gos.residual(*x,r);  // residual is additive

      timing = watch.elapsed();
      // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
      if (gos.trialGridFunctionSpace().gridView().comm().rank()==0)
      {
        if(verbosity > 1) {
          debug_verb << "=== residual assembly (max) " << timing << " s" << std::endl;
          }
      }

      if (row_preconditioner)
      {
        RowNormPreconditioner<GV> preconditioner(m,r);
      }


      // Debug output of matrix and rhs
      if (printMatrix)
      {
        fn_matrix.increment();
        std::stringstream filename;
        filename << fn_matrix.getName();
        filename << ".dat";
        debug_verb << "Writing matrix to file [" << filename.str() << "]" << std::endl;
        Dune::writeMatrixToMatlab(m.base(), filename.str());
        if(printRhs)
        {
          fn_vector.increment();
          filename.str("");
          filename << fn_vector.getName();
          filename << ".dat";
          debug_verb << "Writing rhs to file [" << filename.str() << "]" << std::endl;
          std::ofstream filestream(filename.str());
          Dune::printvector(filestream, r.base(), std::string("#Poisson rhs"), std::string(""));
          filestream.close();
        }
      }

      typename V::ElementType defect = ls.norm(r);

      // compute correction
      watch.reset();
      V z(gos.trialGridFunctionSpace(),0.0);
      typename V::ElementType red = std::min(reduction,defect/mindefect);

      if(verbosity > 1) {
        debug_verb << "=== solving (reduction: " << red << ") ";
      }

      ls.apply(m,z,r,red); // solver makes right hand side consistent
      timing = watch.elapsed();
      // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);

      if(verbosity > 1) {
        debug_verb << timing << " s" << std::endl;
      }

      // damping
      //z *= 0.7;

      // and update
      *x -= z;
    }

    void setPrintMatrix(bool printMatrix_)
    {
        printMatrix = printMatrix_;
    }

    void setPrintRhs(bool printRhs_)
    {
        printRhs = printRhs_;
    }

    void setRowPreconditioner(bool row_preconditioner_)
    {
     row_preconditioner = row_preconditioner_;
    }

    void setVerbosity(int verb)
    {
      verbosity = verb;
    }

    const LS& getLinearSolver()
    {
      return ls;
    }

    const GOS& getGridOperator()
    {
      return gos;
    }

  private:
    const GOS& gos;
    LS& ls;
    V* x;
    typename V::ElementType reduction;
    typename V::ElementType mindefect;
    Dune::PDELab::FilenameHelper fn_matrix;
    Dune::PDELab::FilenameHelper fn_vector;
    bool printMatrix;
    bool printRhs;
    bool row_preconditioner;
    int verbosity;

};


#endif /* AX1_LINEARPROBLEM_HH */
