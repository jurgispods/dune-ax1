#ifndef DUNE_AX1_VECTORDISCRETEGRIDFUNCTIONGRADIENT_HH
#define DUNE_AX1_VECTORDISCRETEGRIDFUNCTIONGRADIENT_HH

#include <vector>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/common/shared_ptr.hh>

namespace Dune {
  namespace PDELab {

  /** \brief DiscreteGridFunctionGradient for vector-valued functions
   *
   *
   * \tparam T Type of PowerGridFunctionSpace
   * \tparam X Type of coefficients vector
   */
  template<typename T, typename X>
  class Ax1VectorDiscreteGridFunctionGradient
    : public GridFunctionInterface<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            //typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType,
            typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
            // dimension of image = dimDomain * #children
            T::CHILDREN*T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimDomain,
            //Dune::FieldVector<
              //typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType,
              Dune::FieldVector<
                typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
                T::CHILDREN*T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimDomain>
            //  ,T::CHILDREN
            //  >
            >,
            Ax1VectorDiscreteGridFunctionGradient<T,X>
          >
  {
    typedef T GFS;

      typedef GridFunctionInterface<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          //typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType,
          typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
          T::CHILDREN*T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimDomain,
          //Dune::FieldVector<
            //typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType,
            Dune::FieldVector<
              typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
              T::CHILDREN*T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimDomain>
            //,T::CHILDREN
           //>
          >,
          Ax1VectorDiscreteGridFunctionGradient<T,X>
        > BaseT;

  public:
    typedef typename BaseT::Traits Traits;
    typedef typename T::template Child<0>::Type ChildType;
    typedef typename ChildType::Traits::FiniteElementType::Traits::
            LocalBasisType::Traits LBTraits;

    typedef typename LBTraits::RangeFieldType RF;
    typedef typename LBTraits::JacobianType JT;

    Ax1VectorDiscreteGridFunctionGradient (const GFS& gfs, const X& x_)
      : pgfs(stackobject_to_shared_ptr(gfs)),
        lfs(gfs),
        lfs_cache(lfs),
        x_view(x_),
        xl(gfs.maxLocalSize()),
        J(gfs.maxLocalSize())
      {
      }

//    inline void evaluate (const typename Traits::ElementType& e,
//              const typename Traits::DomainType& x,
//              typename Traits::RangeType& y) const
//    {
//      myEvaluate(e,x,y);
//
//      lfs.bind(e);
//      lfs.vread(xg, xl);
//
//      // get Jacobian of geometry
//      const typename Traits::ElementType::Geometry::Jacobian& JgeoIT = e.geometry().jacobianInverseTransposed(x);
//
//      for(unsigned int k = 0; k != T::CHILDREN; ++k)
//      {
//        lfs.child(k).finiteElement().localBasis().evaluateJacobian(x, J);
//        y[k] = 0.0;
//
//        /*
//        typename Traits::RangeFieldType gradphi;
//        for(unsigned int i = 0; i < lfs.size(); ++i) {
//          // compute global gradient of shape function i
//          gradphi = 0;
//          JgeoIT.umv(J[i], gradphi);
//
//          // sum up global gradients, weighting them with the appropriate coeff
//          y[k].axpy(xl[lfs.child(k).localIndex(i)], gradphi);
//        }
//        */
//
//        for(unsigned i = 0; i != J.size() ; ++i)
//        {
//          // compute global gradient of shape function i
//          typename Traits::RangeFieldType gradphi(0.0);
//          JgeoIT.umv(J[i][0], gradphi[0]);
//
//          // sum up global gradients, weighting them with the appropriate coeff
//          y[k].axpy(xl[lfs.child(k).localIndex(i)], gradphi);
//        }
//
//      }
//
//  }

    inline void evaluate(const typename Traits::ElementType& e,
        const typename Traits::DomainType& x,
        typename Traits::RangeType& y) const
    {

      // evaluate shape function on the reference element as before
      lfs.bind(e);
      lfs_cache.update();
      x_view.bind(lfs_cache);
      x_view.read(xl);
      x_view.unbind();

      // get Jacobian of geometry
      const typename Traits::ElementType::Geometry::Jacobian
        JgeoIT = e.geometry().jacobianInverseTransposed(x);

      // Loop over PowerLFS and calculate gradient for each child separately
      for(unsigned int k = 0; k != T::CHILDREN; ++k)
      {
        // get local Jacobians/gradients of the shape functions
        std::vector<typename LBTraits::JacobianType> J(lfs.child(k).size());
        lfs.child(k).finiteElement().localBasis().evaluateJacobian(x,J);

        std::vector<Dune::FieldVector<RF,LBTraits::dimDomain> > gradphi(lfs.child(k).size());

        //typedef typename LFSU_POT::Traits::SizeType size_type;

        const int dimDomain =
            T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimDomain;
        Dune::FieldVector<typename Traits::RangeFieldType, dimDomain> childGradient(0.0);

        for (typename LFS::Traits::SizeType i=0; i<lfs.child(k).size(); i++)
        {
          JgeoIT.mv(J[i][0],gradphi[i]);
          //debug_jochen << "JgeoIT = " << JgeoIT << std::endl;
          //debug_jochen << "gradphi = " << Tools::getTypeName(gradphi) << std::endl;
          //debug_jochen << "y = " << Tools::getTypeName(y) << std::endl;
          //debug_jochen << "y[k] = " << Tools::getTypeName(y[k]) << std::endl;
          //debug_jochen << "xl = " << Tools::getTypeName(xl) << std::endl;

          childGradient.axpy(xl[lfs.child(k).localIndex(i)], gradphi[i]);
        }

        // Copy this child's gradient to the large vector containing all gradients
        for(int i=0; i<dimDomain; i++)
        {
          y[k*dimDomain + i] = childGradient[i];
          //debug_jochen << "grad y[" << (k*dimDomain + i) << "] = " << childGradient[i] << std::endl;
        }

        /*
        typename Traits::RangeType gradphi;
        for(unsigned int i = 0; i < lfs.size(); ++i) {
          // compute global gradient of shape function i
          gradphi = 0;
          JgeoIT.umv(J[i][0], gradphi);

          // sum up global gradients, weighting them with the appropriate coeff
          y[k].axpy(xl[lfs.child(k).localIndex(i)], gradphi);
        }*/
      }
    }




      //! get a reference to the GridView
    inline const typename Traits::GridViewType& getGridView () const
    {
      return pgfs->gridView();
    }

  private:
    typedef LocalFunctionSpace<GFS> LFS;
    typedef LFSIndexCache<LFS> LFSCache;
    typedef typename X::template ConstLocalView<LFSCache> XView;

    shared_ptr<GFS const> pgfs;
    mutable LFS lfs;
    mutable LFSCache lfs_cache;
    mutable XView x_view;
    mutable std::vector<RF> xl;
    mutable std::vector<JT> J;
  };



  }
}


#endif // DUNE_AX1_VECTORDISCRETEGRIDFUNCTIONGRADIENT_HH
