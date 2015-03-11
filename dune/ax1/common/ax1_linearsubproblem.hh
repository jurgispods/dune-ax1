#ifndef DUNE_AX1_LINEASUBRPROBLEM_HH
#define DUNE_AX1_LINEASUBRPROBLEM_HH

#include <iostream>

#include <dune/common/timer.hh>
#include <dune/common/deprecated.hh>
#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/backend/solver.hh>

namespace Dune {
  namespace PDELab {

    //===============================================================
    // A class for solving linear stationary problems.
    // It assembles the matrix, computes the right hand side and
    // solves the problem.
    // This is only a first vanilla implementation which has to be improved.
    //===============================================================


    template<typename GO, typename GFS_SUB, typename LS, typename V, int subIndex = 0>
    class Ax1LinearSubProblemSolver
    {
      typedef typename V::ElementType Real;
      typedef typename GO::Traits::Jacobian M;
      typedef typename GO::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
      typedef typename Dune::PDELab::BackendVectorSelector<TrialGridFunctionSpace,Real>::Type W;
      typedef GO GridOperator;
      typedef typename GO::Traits::TrialGridFunctionSpace::Traits::GridViewType GV;

      typedef typename M::Container::block_type M_SUB;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_SUB,typename V::field_type>::Type U_SUB;

      const int precision = 32;

    public:
      typedef StationaryLinearProblemSolverResult<double> Result;

      Ax1LinearSubProblemSolver(const GO& go, const GFS_SUB& gfsSub, V& x, LS& ls, const char* basename_, typename V::ElementType reduction,
          typename V::ElementType min_defect = 1e-99, int verbose=1) DUNE_DEPRECATED_MSG("Use StationaryLinearProblemSolver(const GO&, LS&, V&, ...) instead.")
        : _go(go)
        , _gfsSub(gfsSub)
        , _ls(ls)
        , _x(&x)
        , _reduction(reduction)
        , _min_defect(min_defect)
        , _hanging_node_modifications(false)
        , _keep_matrix(true)
        , _verbose(verbose)
        , fn_matrix(std::string(basename_) + "_mat")
        , fn_vector(std::string(basename_) + "_rhs")
        , printMatrix(false)
        , printRhs(false)
        , row_preconditioner(false)
      {}

      Ax1LinearSubProblemSolver(const GO& go, const GFS_SUB& gfsSub, LS& ls, V& x, const char* basename_, typename V::ElementType reduction,
          typename V::ElementType min_defect = 1e-99, int verbose=1)
        : _go(go)
        , _gfsSub(gfsSub)
        , _ls(ls)
        , _x(&x)
        , _reduction(reduction)
        , _min_defect(min_defect)
        , _hanging_node_modifications(false)
        , _keep_matrix(true)
        , _verbose(verbose)
        , fn_matrix(std::string(basename_) + "_mat")
        , fn_vector(std::string(basename_) + "_rhs")
        , printMatrix(false)
        , printRhs(false)
        , row_preconditioner(false)
      {}

      Ax1LinearSubProblemSolver (const GO& go, const GFS_SUB& gfsSub, LS& ls, const char* basename_,
          typename V::ElementType reduction, typename V::ElementType min_defect = 1e-99, int verbose=1)
        : _go(go)
        , _gfsSub(gfsSub)
        , _ls(ls)
        , _x()
        , _reduction(reduction)
        , _min_defect(min_defect)
        , _hanging_node_modifications(false)
        , _keep_matrix(true)
        , _verbose(verbose)
        , fn_matrix(std::string(basename_) + "_mat")
        , fn_vector(std::string(basename_) + "_rhs")
        , printMatrix(false)
        , printRhs(false)
        , row_preconditioner(false)
      {}

      //! Construct a StationaryLinearProblemSolver for the given objects and read parameters from a ParameterTree.
      /**
       * This constructor reads the parameter controlling its operation from a passed-in ParameterTree
       * instead of requiring the user to specify all of them as individual constructor parameters.
       * Currently the following parameters are read:
       *
       * Name                       | Default Value | Explanation
       * -------------------------- | ------------- | -----------
       * reduction                  |               | Required relative defect reduction
       * min_defect                 | 1e-99         | minimum absolute defect at which to stop
       * hanging_node_modifications | false         | perform required transformations for hanging nodes
       * keep_matrix                | true          | keep matrix between calls to apply() (but reassemble values every time)
       * verbosity                  | 1             | control amount of debug output
       *
       * Apart from reduction, all parameters have a default value and are optional.
       * The actual reduction for a call to apply() is calculated as r = max(reduction,min_defect/start_defect),
       * where start defect is the norm of the residual of x.
       */
      Ax1LinearSubProblemSolver(const GO& go, const GFS_SUB& gfsSub, LS& ls, V& x, const char* basename_, const ParameterTree& params)
        : _go(go)
        , _gfsSub(gfsSub)
        , _ls(ls)
        , _x(&x)
        , _reduction(params.get<typename V::ElementType>("reduction"))
        , _min_defect(params.get<typename V::ElementType>("min_defect",1e-99))
        , _hanging_node_modifications(params.get<bool>("hanging_node_modifications",false))
        , _keep_matrix(params.get<bool>("keep_matrix",true))
        , _verbose(params.get<int>("verbosity",1))
        , fn_matrix(std::string(basename_) + "_mat")
        , fn_vector(std::string(basename_) + "_rhs")
        , printMatrix(false)
        , printRhs(false)
        , row_preconditioner(false)
      {}

      //! Construct a StationaryLinearProblemSolver for the given objects and read parameters from a ParameterTree.
      /**
       * This constructor reads the parameter controlling its operation from a passed-in ParameterTree
       * instead of requiring the user to specify all of them as individual constructor parameters.
       * Currently the following parameters are read:
       *
       * Name                       | Default Value | Explanation
       * -------------------------- | ------------- | -----------
       * reduction                  |               | Required relative defect reduction
       * min_defect                 | 1e-99         | minimum absolute defect at which to stop
       * hanging_node_modifications | false         | perform required transformations for hanging nodes
       * keep_matrix                | true          | keep matrix between calls to apply() (but reassemble values every time)
       * verbosity                  | 1             | control amount of debug output
       *
       * Apart from reduction, all parameters have a default value and are optional.
       * The actual reduction for a call to apply() is calculated as r = max(reduction,min_defect/start_defect),
       * where start defect is the norm of the residual of x.
       */
      Ax1LinearSubProblemSolver(const GO& go, const GFS_SUB& gfsSub, LS& ls, const char* basename_, const ParameterTree& params)
        : _go(go)
        , _gfsSub(gfsSub)
        , _ls(ls)
        , _x()
        , _reduction(params.get<typename V::ElementType>("reduction"))
        , _min_defect(params.get<typename V::ElementType>("min_defect",1e-99))
        , _hanging_node_modifications(params.get<bool>("hanging_node_modifications",false))
        , _keep_matrix(params.get<bool>("keep_matrix",true))
        , _verbose(params.get<int>("verbosity",1))
        , fn_matrix(std::string(basename_) + "_mat")
        , fn_vector(std::string(basename_) + "_rhs")
        , printMatrix(false)
        , printRhs(false)
        , row_preconditioner(false)
      {}

      //! Set whether the solver should apply the necessary transformations for calculations on hanging nodes.
      void setHangingNodeModifications(bool b)
      {
        _hanging_node_modifications=b;
      }

      //! Return whether the solver performs the necessary transformations for calculations on hanging nodes.
      bool hangingNodeModifications() const
      {
        return _hanging_node_modifications;
      }

      //! Set whether the jacobian matrix should be kept across calls to apply().
      void setKeepMatrix(bool b)
      {
        _keep_matrix = b;
      }

      //! Return whether the jacobian matrix is kept across calls to apply().
      bool keepMatrix() const
      {
        return _keep_matrix;
      }

      const Result& result() const
      {
        return _res;
      }

      void apply(V& x) {
        _x = &x;
        apply();
      }

      void apply ()
      {
        Dune::Timer watch;
        double timing,assembler_time=0;

        // assemble matrix; optional: assemble only on demand!
        watch.reset();

        if (!_jacobian)
          {
            _jacobian = make_shared<M>(_go);
            timing = watch.elapsed();
            if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
              std::cout << "=== matrix setup (max) " << timing << " s" << std::endl;
            watch.reset();
            assembler_time += timing;
          }
        else if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
          std::cout << "=== matrix setup skipped (matrix already allocated)" << std::endl;

        (*_jacobian) = Real(0.0);
        if (_hanging_node_modifications)
          {
            Dune::PDELab::set_shifted_dofs(_go.localAssembler().trialConstraints(),0.0,*_x); // set hanging node DOFs to zero
            _go.localAssembler().backtransform(*_x); // interpolate hanging nodes adjacent to Dirichlet nodes
          }
        _go.jacobian(*_x,*_jacobian);

        timing = watch.elapsed();
        // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
        if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
          std::cout << "=== matrix assembly (max) " << timing << " s" << std::endl;
        assembler_time += timing;

        // assemble residual
        watch.reset();

        W r(_go.testGridFunctionSpace(),0.0);
        _go.residual(*_x,r);  // residual is additive

        timing = watch.elapsed();
        // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
        if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
          std::cout << "=== residual assembly (max) " << timing << " s" << std::endl;
        assembler_time += timing;
        _res.assembler_time = assembler_time;

        /* Here comes the crucial part: Since we are solving a subsystem of the complete system, we need to copy
         * the relevant part from the vectors r, z into new ISTLVectorBackends, on which the linear solver
         * will operate!
         */
        // Reference to sub-block in shared pointer to ISTLMatrixBackend (...)
        M_SUB& mSub = Dune::PDELab::istl::raw(*_jacobian)[subIndex][subIndex];

        // TODO Why not create a reference here?
        U_SUB rSub(_gfsSub);
        rSub.base() = Dune::PDELab::istl::raw(r)[subIndex];

        if (row_preconditioner)
        {
          RowNormPreconditioner<GV> preconditioner(mSub,istl::raw(rSub));
        }

        // Debug output of matrix and rhs
        if (printMatrix)
        {
          fn_matrix.increment();
          std::stringstream filename;
          filename << fn_matrix.getName();
          filename << ".dat";
          debug_verb << "Writing matrix to file [" << filename.str() << "]" << std::endl;
          Dune::writeMatrixToMatlab(mSub, filename.str(), precision);

          std::stringstream full_filename;
          full_filename << fn_matrix.getName();
          full_filename << "_full.dat";

          Dune::writeMatrixToMatlab(istl::raw(*_jacobian), full_filename.str(), precision);

          if(printRhs)
          {
            fn_vector.increment();
            filename.str("");
            filename << fn_vector.getName();
            filename << ".dat";
            debug_verb << "Writing rhs to file [" << filename.str() << "]" << std::endl;
            std::ofstream filestream(filename.str());
            Dune::printvector(filestream, rSub.base(), std::string("% block #" + std::to_string(subIndex) + "rhs"),
                std::string(""));
            filestream.close();

            std::stringstream full_filename;
            full_filename << fn_vector.getName();
            full_filename << "_full.dat";

            std::ofstream filestream_full(full_filename.str());
            Dune::printvector(filestream_full, r.base(), std::string("% full rhs"),
                std::string(""));
            filestream_full.close();
          }
        }

        typename V::ElementType defect = _ls.norm(rSub);

        if(std::isnan(defect))
          DUNE_THROW(Dune::Exception, "Initial residual is NaN!");

        std::cout << "Initial defect: " << defect << std::endl;
        std::cout << "Size of subproblem residual: " << rSub.flatsize() << " / " << rSub.N() << std::endl;

        // compute correction
        watch.reset();

        //V z(_go.trialGridFunctionSpace(),0.0);
        U_SUB zSub(_gfsSub, 0.0);

        typename V::ElementType red = std::max(_reduction,_min_defect/defect);

        if (_go.trialGridFunctionSpace().gridView().comm().rank()==0)
          std::cout << "=== solving (reduction: " << red << ") ";

        // Save rhs to later calculate achieved linear reduction
        U_SUB rhs_save(rSub);

        //_ls.apply(*_jacobian,z,r,red); // solver makes right hand side consistent
        _ls.apply(mSub,zSub,rSub,red); // solver makes right hand side consistent

        _linear_solver_result = _ls.result();
        timing = watch.elapsed();
        // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
        if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
          std::cout << timing << " s" << std::endl;
        _res.linear_solver_time = timing;

        _res.converged = _linear_solver_result.converged;
        _res.iterations = _linear_solver_result.iterations;
        _res.elapsed = _linear_solver_result.elapsed;
        _res.reduction = _linear_solver_result.reduction;
        _res.conv_rate = _linear_solver_result.conv_rate;
        _res.first_defect = static_cast<double>(defect);
        _res.defect = static_cast<double>(defect)*_linear_solver_result.reduction;
        _res.linear_solver_iterations = _linear_solver_result.iterations;

        debug_info << "End defect: " << _res.defect << std::endl;
        debug_info << "Achieved reduction (solver info): " << _res.reduction << std::endl;

        U_SUB linearResidual(rSub);
        mSub.mv(istl::raw(zSub), istl::raw(linearResidual));
        istl::raw(linearResidual) -= istl::raw(rhs_save);  // linear residual after solve, linearResidual = Az-r

        if(_verbose)
        {
          debug_info << "          M Frobenius norm: " << mSub.frobenius_norm() << std::endl; // norm of system matrix

          typename V::ElementType zNorm = istl::raw(zSub).two_norm();
          debug_info << "          z 2-norm: " << zNorm << std::endl; // norm of solution update
          //debug_info << "r 2-norm: " << rhs.two_norm() << std::endl;

          typename V::ElementType rhsNorm = istl::raw(rhs_save).two_norm();
          // initial residual norm == initial linear defect
          debug_info << "          initial linear defect: " << rhsNorm << std::endl;
          //norm of linear residual after solve

          typename V::ElementType linResNorm = istl::raw(linearResidual).two_norm();
          debug_info << "          linear defect after solve: " << linResNorm << std::endl;
          // Linear reduction = linear defect / initial linear defect
          debug_info << "          Achieved linear reduction: "
            << (linResNorm / rhsNorm) << std::endl;

          if(std::isnan(zNorm) || std::isnan(rhsNorm) || std::isnan(linResNorm))
            DUNE_THROW(Dune::Exception, "One of {zNorm, rhsNorm, linResNorm} is NaN!");
        }

        // and update
        if (_hanging_node_modifications)
          Dune::PDELab::set_shifted_dofs(_go.localAssembler().trialConstraints(),0.0,*_x); // set hanging node DOFs to zero

        // Apply update only to subsystem coefficients!
        Dune::PDELab::istl::raw(*_x)[subIndex] -= zSub;

        if (_hanging_node_modifications)
          _go.localAssembler().backtransform(*_x); // interpolate hanging nodes adjacent to Dirichlet nodes

        if (!_keep_matrix)
          _jacobian.reset();
      }

      //! Discard the stored Jacobian matrix.
      void discardMatrix()
      {
        if(_jacobian)
          _jacobian.reset();
      }

      const Dune::PDELab::LinearSolverResult<double>& ls_result() const{
        return _linear_solver_result;
      }

      void setMatrix(shared_ptr<M> matrix)
      {
        _jacobian = matrix;
      }

      void setVerbosity(int verb)
      {
        _verbose = verb;
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


    private:
      const GO& _go;
      const GFS_SUB& _gfsSub;
      LS& _ls;
      V* _x;
      shared_ptr<M> _jacobian;
      typename V::ElementType _reduction;
      typename V::ElementType _min_defect;
      Dune::PDELab::LinearSolverResult<double> _linear_solver_result;
      Result _res;
      bool _hanging_node_modifications;
      bool _keep_matrix;
      int _verbose;

      Dune::PDELab::FilenameHelper fn_matrix;
      Dune::PDELab::FilenameHelper fn_vector;
      bool printMatrix;
      bool printRhs;
      bool row_preconditioner;

    };

  } // namespace PDELab
} // namespace Dune

#endif
