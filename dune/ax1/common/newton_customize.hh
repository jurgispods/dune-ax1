#include"rownorm_preconditioner.hh"
#include<dune/pdelab/newton/newton.hh>

namespace Dune
{
  namespace PDELab
  {

    template<class GOS, class S, class TrlV, class TstV = TrlV>
    class CustomizeNewton: public NewtonSolver<GOS, S, TrlV, TstV> , public NewtonTerminate<GOS,
        TrlV, TstV> , public NewtonLineSearch<GOS, TrlV, TstV> , public NewtonPrepareStep<GOS,
        TrlV, TstV>
    {
        typedef GOS GridOperator;
        typedef S Solver;
        typedef TrlV TrialVector;
        typedef TstV TestVector;

        typedef typename GridOperator::Traits::TrialGridFunctionSpace::Traits::GridViewType GV;

        typedef typename TestVector::ElementType RFType;
        typedef typename GOS::template MatrixContainer<RFType>::Type Matrix;

      public:
        CustomizeNewton(GridOperator& go, TrialVector& u_, Solver& solver_) :
          NewtonBase<GOS, TrlV, TstV> (go, u_), NewtonSolver<GOS, S, TrlV, TstV> (go, u_, solver_),
              NewtonTerminate<GOS, TrlV, TstV> (go, u_),
              NewtonLineSearch<GOS, TrlV, TstV> (go, u_),
              NewtonPrepareStep<GOS, TrlV, TstV> (go, u_), min_linear_reduction(1e-3),
              reassemble_threshold(0.0), row_preconditioner(true)

        {
        }
        CustomizeNewton(GridOperator& go, Solver& solver_) :
          NewtonBase<GOS, TrlV, TstV> (go), NewtonSolver<GOS, S, TrlV, TstV> (go, solver_),
              NewtonTerminate<GOS, TrlV, TstV> (go), NewtonLineSearch<GOS, TrlV, TstV> (go),
              NewtonPrepareStep<GOS, TrlV, TstV> (go), min_linear_reduction(1e-3),
              reassemble_threshold(0.0), row_preconditioner(true)
        {
        }

        void setMinLinearReduction(RFType min_linear_reduction_)
        {
          min_linear_reduction = min_linear_reduction_;
        }

        void setReassembleThreshold(RFType reassemble_threshold_)
        {
          reassemble_threshold = reassemble_threshold_;
        }

        void setRowPreconditioner(bool row_preconditioner_)
        {
          row_preconditioner = row_preconditioner_;
        }

        void prepare_step(Matrix& A, TestVector& r)
        {

          TrialVector z(this->gridoperator.trialGridFunctionSpace());
          z = 1.0;
          this->reassembled = false;

          if (this->res.defect / this->prev_defect > reassemble_threshold)
          {
            if (this->verbosity_level >= 3) std::cout << "      Reassembling matrix..."
                << std::endl;
            A = 0.0; // TODO: Matrix interface
            this->gridoperator.jacobian(*this->u, A);
            if (row_preconditioner) RowNormPreconditioner<GV> preconditioner(A, r);
            this->reassembled = true;
          }

          this->linear_reduction = std::min(min_linear_reduction,
              this->res.defect * this->res.defect / (this->prev_defect * this->prev_defect));

          this->prev_defect = this->res.defect;

          ios_base_all_saver restorer(std::cout); // store old ios flags

          if (this->verbosity_level >= 4) std::cout << "      requested linear reduction:       "
              << std::setw(12) << std::setprecision(4) << std::scientific << this->linear_reduction
              << std::endl;

        }
      private:
        RFType min_linear_reduction;
        RFType reassemble_threshold;
        bool row_preconditioner;

    };

  }
}

