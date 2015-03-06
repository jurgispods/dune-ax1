/*
 * ax1_solverbackend.hh
 *
 *  Created on: Jun 25, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_SOLVERBACKEND_HH
#define DUNE_AX1_SOLVERBACKEND_HH


namespace Dune{
  namespace PDELab{

#if HAVE_SUPERLU
    /**
     * @brief Solver backend using SuperLU as a direct solver.
     */
    class SubBlock_ISTLBackend_SEQ_SuperLU
      : public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] verbose_ print messages if true
      */
      explicit SubBlock_ISTLBackend_SEQ_SuperLU (int verbose_=1)
        : verbose(verbose_)
      {}


      /*! \brief make a linear solver object

        \param[in] maxiter Maximum number of allowed steps (ignored)
        \param[in] verbose_ print messages if true
      */
      SubBlock_ISTLBackend_SEQ_SuperLU (int maxiter, int verbose_)
        : verbose(verbose_)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W, class E>
      void apply(M& A, V& z, W& r, E reduction)
      {
        Dune::SuperLU<M> solver(A, verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(z, r, stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

      template<class V>
      typename Dune::template FieldTraits<typename V::field_type >::real_type norm(const V& v) const
      {
        return v.two_norm();
      }

    private:
      int verbose;
    };
#endif // HAVE_SUPERLU


    template<class GO, int s, template<class,class,class,int> class Preconditioner,
             //Adrian: template<class> class Solver>
             template<class,class,class> class Solver>
    class ISTLBackend_GMRES_AMG : public LinearResultStorage
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
      typedef istl::ParallelHelper<GFS> PHELPER;
      typedef typename GO::Traits::Jacobian M;
      typedef typename M::BaseT MatrixType;
      typedef typename GO::Traits::Domain V;
      typedef typename V::BaseT VectorType;
      typedef typename istl::CommSelector<s,Dune::MPIHelper::isFake>::type Comm;
#if HAVE_MPI
      typedef Preconditioner<MatrixType,VectorType,VectorType,1> Smoother;
      typedef Dune::BlockPreconditioner<VectorType,VectorType,Comm,Smoother> ParSmoother;
      typedef Dune::OverlappingSchwarzOperator<MatrixType,VectorType,VectorType,Comm> Operator;
#else
      typedef Preconditioner<MatrixType,VectorType,VectorType,1> ParSmoother;
      typedef Dune::MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#endif
      typedef typename Dune::Amg::SmootherTraits<ParSmoother>::Arguments SmootherArgs;
      typedef Dune::Amg::AMG<Operator,VectorType,ParSmoother,Comm> AMG;

    public:

      /**
       * @brief Parameters object to customize matrix hierachy building.
       */
      typedef Dune::Amg::Parameters Parameters;

      ISTLBackend_GMRES_AMG(const GFS& gfs_, unsigned maxiter_=5000,
                      int verbose_=1, bool reuse_=false,
                      bool usesuperlu_=true)
        : gfs(gfs_), phelper(gfs,verbose_), maxiter(maxiter_), params(15,2000),
          verbose(verbose_), reuse(reuse_), firstapply(true),
          usesuperlu(usesuperlu_)
      {
				params.setDefaultValuesIsotropic(GFS::Traits::GridViewType::Traits::Grid::dimension);
        params.setDebugLevel(verbose_);
#if !HAVE_SUPERLU
        if (gfs.gridView().comm().rank() == 0 && usesuperlu == true)
          {
            std::cout << "WARNING: You are using AMG without SuperLU!"
                      << " Please consider installing SuperLU,"
                      << " or set the usesuperlu flag to false"
                      << " to suppress this warning." << std::endl;
          }
#endif
      }

       /*! \brief set AMG parameters

        \param[in] params_ a parameter object of Type Dune::Amg::Parameters
      */
      void setParameters(const Parameters& params_)
      {
        params = params_;
      }

      void setparams(Parameters params_) DUNE_DEPRECATED_MSG("setparams() is deprecated, use setParameters() instead")
      {
        params = params_;
      }

      /**
       * @brief Get the parameters describing the behaviuour of AMG.
       *
       * The returned object can be adjusted to ones needs and then can be
       * reset using setParameters.
       * @return The object holding the parameters of AMG.
       */
      const Parameters& parameters() const
      {
        return params;
      }

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      typename V::ElementType norm (const V& v) const
      {
        typedef OverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        return psp.norm(v);
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      void apply(M& A, V& z, V& r, typename V::ElementType reduction)
      {
        Timer watch;
        Comm oocc(gfs.gridView().comm());
        MatrixType& mat=A.base();
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MatrixType,
          Dune::Amg::FirstDiagonal> > Criterion;
#if HAVE_MPI
        phelper.createIndexSetAndProjectForAMG(A, oocc);
        Operator oop(mat, oocc);
Dune::OverlappingSchwarzScalarProduct<VectorType,Comm> sp(oocc);
#else
        Operator oop(mat);
        Dune::SeqScalarProduct<VectorType> sp;
#endif
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;
        Criterion criterion(params);
        stats.tprepare=watch.elapsed();
        watch.reset();

        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        //only construct a new AMG if the matrix changes
        if (reuse==false || firstapply==true){
          amg.reset(new AMG(oop, criterion, smootherArgs, oocc));
          firstapply = false;
          stats.tsetup = watch.elapsed();
          stats.levels = amg->maxlevels();
          stats.directCoarseLevelSolver=amg->usesDirectCoarseLevelSolver();
        }
        watch.reset();
        //Adrian: Solver<VectorType> solver(oop,sp,*amg,reduction,maxiter,verb);
        int restart=100;
        bool recalcDefect = false;

        Solver<VectorType,VectorType,VectorType> solver(oop,sp,*amg,reduction,restart,maxiter,verb,recalcDefect);
        Dune::InverseOperatorResult stat;

        solver.apply(istl::raw(z),istl::raw(r),stat);
        stats.tsolve= watch.elapsed();
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

      /**
       * @brief Get statistics of the AMG solver (no of levels, timings).
       * @return statistis of the AMG solver.
       */
      const ISTLAMGStatistics& statistics() const
      {
        return stats;
      }

    private:
      const GFS& gfs;
      PHELPER phelper;
      unsigned maxiter;
      Parameters params;
      int verbose;
      bool reuse;
      bool firstapply;
      bool usesuperlu;
      shared_ptr<AMG> amg;
      ISTLAMGStatistics stats;
    };


    /**
     * @brief Overlapping parallel BiCGStab solver preconditioned with AMG smoothed by SSOR.
     * @tparam GO The type of the grid operator
     * (or the fakeGOTraits class for the old grid operator space).
     * @tparam s The bits to use for the globale index.
     */
    template<class GO, int s=96>
    class ISTLBackend_GMRES_AMG_SSOR
      : public ISTLBackend_GMRES_AMG<GO,
                               s,
                               Dune::SeqSSOR,
                               Dune::RestartedGMResSolver> //Dune::BiCGSTABSolver>
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
    public:
      /**
       * @brief Constructor
       * @param gfs_ The grid function space used.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       * @param reuse_ Set true, if the Matrix to be used is always identical
       * (AMG aggregation is then only performed once).
       * @param usesuperlu_ Set false, to suppress the no SuperLU warning
       */
      ISTLBackend_GMRES_AMG_SSOR(const GFS& gfs_, unsigned maxiter_=5000,
                                int verbose_=1, bool reuse_=false,
                                bool usesuperlu_=true)
        : ISTLBackend_GMRES_AMG<GO, s, Dune::SeqSSOR, Dune::RestartedGMResSolver>//Dune::BiCGSTABSolver>
          (gfs_, maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };


    /**
     * @brief Overlapping parallel BiCGStab solver preconditioned with AMG smoothed by ILU0.
     * @tparam GO The type of the grid operator
     * (or the fakeGOTraits class for the old grid operator space).
     * @tparam s The bits to use for the globale index.
     */
    template<class GO, int s=96>
    class ISTLBackend_GMRES_AMG_ILU0
      : public ISTLBackend_GMRES_AMG<GO,
                               s,
                               Dune::SeqILU0,
                               Dune::RestartedGMResSolver> //Dune::BiCGSTABSolver>
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
    public:
      /**
       * @brief Constructor
       * @param gfs_ The grid function space used.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       * @param reuse_ Set true, if the Matrix to be used is always identical
       * (AMG aggregation is then only performed once).
       * @param usesuperlu_ Set false, to suppress the no SuperLU warning
       */
      ISTLBackend_GMRES_AMG_ILU0(const GFS& gfs_, unsigned maxiter_=5000,
                                int verbose_=1, bool reuse_=false,
                                bool usesuperlu_=true)
        : ISTLBackend_GMRES_AMG<GO, s, Dune::SeqILU0, Dune::RestartedGMResSolver>//Dune::BiCGSTABSolver>
          (gfs_, maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };


    // wrapped sequential preconditioner
    template<class CC, class GFS, class P>
    class Ax1OverlappingWrappedPreconditioner
      : public Dune::Preconditioner<typename Dune::PDELab::BackendVectorSelector<GFS,typename P::domain_type::field_type>::Type,
                                    typename Dune::PDELab::BackendVectorSelector<GFS,typename P::range_type::field_type>::Type>
    {
    public:
      //! \brief The domain type of the preconditioner.
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,typename P::domain_type::field_type>::Type
      domain_type;
      //! \brief The range type of the preconditioner.
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,typename P::range_type::field_type>::Type
      range_type;

      // define the category
      enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::overlapping
      };

      //! Constructor.
      Ax1OverlappingWrappedPreconditioner (const GFS& gfs_, P& prec_, const CC& cc_,
                                        const istl::ParallelHelper<GFS>& helper_,
                                        const std::vector<int>& perm_,
                                        const std::vector<int>& inv_perm_)
        : gfs(gfs_), prec(prec_), cc(cc_), helper(helper_),
          perm(perm_),
          inv_perm(inv_perm_)
        {}

      /*!
        \brief Prepare the preconditioner.
      */
      virtual void pre (domain_type& x, range_type& b)
      {
        prec.pre(x,b);
      }

      /*!
        \brief Apply the preconditioner.
      */
      virtual void apply (domain_type& v, const range_type& d)
      {
        bool doReorder = (perm.size() > 0);

        range_type dd(d);
        set_constrained_dofs(cc,0.0,dd);

        prec.apply(istl::raw(v),istl::raw(dd));

        // When using reordered DOF vector, restore original ordering for communication!
        if(doReorder)
        {
          // Permute v back to original ordering v_orig
          domain_type v_orig(v);
          for(int i=0; i<perm.size(); i++)
          {
            istl::raw(v_orig)[perm[i]] = istl::raw(v)[i];
          }

          Dune::PDELab::AddDataHandle<GFS,domain_type> adddh(gfs,v_orig);

          if (gfs.gridView().comm().size()>1)
            gfs.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);

          // Permute v_orig again to reordered v
          for(int i=0; i<perm.size(); i++)
          {
            istl::raw(v)[inv_perm[i]] = istl::raw(v_orig)[i];
          }
        } else {
          Dune::PDELab::AddDataHandle<GFS,domain_type> adddh(gfs,v);

          if (gfs.gridView().comm().size()>1)
            gfs.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
        }
      }

      /*!
        \brief Clean up.
      */
      virtual void post (domain_type& x)
      {
        prec.post(istl::raw(x));
      }

    private:
      const GFS& gfs;
      P& prec;
      const CC& cc;
      const istl::ParallelHelper<GFS>& helper;
      const std::vector<int>& perm;
      const std::vector<int>& inv_perm;
    };

    template<typename GFS, typename X>
    class Ax1OVLPScalarProduct
      : public ScalarProduct<X>
    {
    public:
      enum {category=Dune::SolverCategory::overlapping};
      Ax1OVLPScalarProduct(const OVLPScalarProductImplementation<GFS>& implementation_,
          const std::vector<int>& perm_, const std::vector<int>& inv_perm_)
        : implementation(implementation_),
          perm(perm_),
          inv_perm(inv_perm_)
      {}

      virtual typename X::BaseT::field_type dot(const X& x, const X& y)
      {
        bool doReorder = (perm.size() > 0);

        // When using reordered DOF vector, restore original ordering for communication!
        if(doReorder)
        {
          // Permute v back to original ordering v_orig
          X x_orig(x);
          X y_orig(y);
          for(int i=0; i<perm.size(); i++)
          {
            istl::raw(x_orig)[perm[i]] = istl::raw(x)[i];
            istl::raw(y_orig)[perm[i]] = istl::raw(y)[i];
          }
          return implementation.dot(x_orig,y_orig);
        }

        return implementation.dot(x,y);
      }

      virtual typename X::BaseT::field_type norm (const X& x)
      {
        return sqrt(static_cast<double>(this->dot(x,x)));
      }

    private:
      const OVLPScalarProductImplementation<GFS>& implementation;
      const std::vector<int>& perm;
      const std::vector<int>& inv_perm;
    };

    // Base class for ILU0 as preconditioner, using a reordered solution vector; The index permutations perm and inv_perm
    // are handed over in the constructor, they work only for a flat (non-blocked) solution vector
    template<class GFS, class C,
             template<class> class Solver>
    class Ax1ReorderedISTLBackend_OVLP_ILU0_Base
      : public OVLPScalarProductImplementation<GFS>, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] c_ a constraints object
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      Ax1ReorderedISTLBackend_OVLP_ILU0_Base (const GFS& gfs_, const C& c_,
          const std::vector<int>& perm_, const std::vector<int>& inv_perm_,
          unsigned maxiter_=5000, int verbose_=1)
      : OVLPScalarProductImplementation<GFS>(gfs_), gfs(gfs_), c(), maxiter(maxiter_), verbose(verbose_)
        , perm(perm_)
        , inv_perm(inv_perm_)
      {
        if(perm.size() > 0)
        {
          debug_jochen << "Permuting the constraints container entries according to matrix reordering!" << std::endl;
          for(typename C::const_iterator cit = c_.begin(); cit!=c_.end(); ++cit)
          {
            typename C::key_type key = cit->first;

            debug_jochen << "Permuting constraint index: " << key[0] << " -> " << inv_perm[key[0]] << std::endl;
            key[0] = inv_perm[key[0]];
            c.insert(typename C::value_type(key,cit->second));
          }

          // Print constraints indices as used in solver backends
          for(typename C::const_iterator cit = c.begin(); cit!=c.end(); ++cit)
          {
            debug_verb << "Constraint DOF: 'global col' #" << (cit->first) << std::endl;
          }
        } else {
          debug_warn << "[Ax1ReorderedISTLBackend_OVLP_ILU0_Base()] Permutation empty, using original constraints container!" << std::endl;
          c = c_;
        }
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef OverlappingOperator<C,M,V,W> POP;
        POP pop(c,A);

        // Custom scalar product gets index permutations
        typedef Ax1OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this,perm,inv_perm);
        typedef SeqILU0<typename M::BaseT,typename V::BaseT,typename W::BaseT,1> SeqPrec;
        SeqPrec seqprec(istl::raw(A),1.0);

        // Custom preconditioner wrapper get index permutations
        typedef Ax1OverlappingWrappedPreconditioner<C,GFS,SeqPrec> WPREC;
        WPREC wprec(gfs,seqprec,c,this->parallelHelper(),perm,inv_perm);
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Solver<V> solver(pop,psp,wprec,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }
    private:
      const GFS& gfs;
      C c;
      unsigned maxiter;
      int steps;
      int verbose;
      const std::vector<int>& perm;
      const std::vector<int>& inv_perm;
    };


    /**
     * @brief Overlapping parallel BiCGStab solver with ILU0 preconditioner
     * @tparam GFS The Type of the GridFunctionSpace.
     * @tparam CC The Type of the Constraints Container.
     */
    template<class GFS, class CC>
    class Ax1ReorderedISTLBackend_OVLP_BCGS_ILU0
      : public Ax1ReorderedISTLBackend_OVLP_ILU0_Base<GFS,CC,Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] cc a constraints container object
        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
        Ax1ReorderedISTLBackend_OVLP_BCGS_ILU0 (const GFS& gfs, const CC& cc,
            const std::vector<int>& perm, const std::vector<int>& inv_perm,
            unsigned maxiter=5000, int verbose=1)
        : Ax1ReorderedISTLBackend_OVLP_ILU0_Base<GFS,CC,Dune::BiCGSTABSolver>(gfs, cc, perm, inv_perm, maxiter, verbose)
      {}
    };
  }
}


#endif /* DUNE_AX1_SOLVERBACKEND_HH */
