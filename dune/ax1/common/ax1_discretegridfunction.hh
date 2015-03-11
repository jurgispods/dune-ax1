/*
 * ax1_discretegridfunction.hh
 *
 *  Created on: Feb 9, 2012
 *      Author: jpods
 */
#ifndef DUNE_AX1_DISCRETEGRIDFUNCTION_HH
#define DUNE_AX1_DISCRETEGRIDFUNCTION_HH

#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

#include <dune/ax1/common/vectordiscretegridfunctiongradient.hh>



/** \brief Specialization of PDELab's VectorDiscreteGridFunction which
 * evaluates the function only on non-membrane elements and returns
 * 0 otherwise; use this only for concentration vectors!
 */
template<typename PHYSICS, typename T, typename X, std::size_t dimR = T::CHILDREN>
class Ax1VectorDiscreteGridFunction
  : public Dune::PDELab::VectorDiscreteGridFunction<T,X,dimR>
{
  typedef T GFS;
  typedef Dune::PDELab::VectorDiscreteGridFunction<T,X,dimR> BaseT;

public:
    typedef typename BaseT::Traits Traits;

    //! construct
    /**
     * \param gfs   GridFunctionSpace.
     * \param x_    Coefficient vector.
     * \param start Number of first child of gfs to use.
     */
    Ax1VectorDiscreteGridFunction(PHYSICS& physics_, const GFS& gfs, const X& x_, std::size_t start = 0)
    : BaseT(gfs, x_, start),
      physics(physics_)
    {
    }

    //! construct
    /**
     * \param gfs    GridFunctionSpace.
     * \param x_     Coefficient vector.
     * \param remap_ Subscriptable entity (i.e. a container, array, or
     *               pointer) with at least dimR entries.  The relevant
     *               entries are copied.
     *
     * \note If \c i denotes a component of the resulting grid function,
     *       then remap_[i] denotes the corresponding child of the
     *       gridfunctionspace.
     */
    template<class Remap>
    Ax1VectorDiscreteGridFunction(PHYSICS& physics_, const GFS& gfs, const X& x_, const Remap &remap_)
    : BaseT(gfs, x_, remap_),
      physics(physics_)
    {
    }

    // Evalute grid function; return 0 if the element is part of the membrane!
    inline void evaluate (const typename Traits::ElementType& e,
              const typename Traits::DomainType& x,
              typename Traits::RangeType& y) const
    {
      y = 0.0;
      if(not physics.isMembrane(e))
      {
        BaseT::evaluate(e,x,y);
      }
    }

private:
  PHYSICS& physics;
};


/** \brief Specialization of VectorDiscreteGridFunctionGradient which
 * evaluates the function only on non-membrane elements and returns
 * 0 otherwise; use this only for concentration vectors!
 */
template<typename PHYSICS, typename T, typename X>
class Ax1VectorDiscreteGridFunctionGradient
  :  public Dune::PDELab::VectorDiscreteGridFunctionGradient<T,X>
{
  typedef T GFS;
  typedef Dune::PDELab::VectorDiscreteGridFunctionGradient<T,X> BaseT;

public:
    typedef typename BaseT::Traits Traits;

    //! construct
    /**
     * \param gfs   GridFunctionSpace.
     * \param x_    Coefficient vector.
     * \param start Number of first child of gfs to use.
     */
    Ax1VectorDiscreteGridFunctionGradient(PHYSICS& physics_, const GFS& gfs, const X& x_)
    : BaseT(gfs, x_),
      physics(physics_)
    {
    }


    // Evalute grid function; return 0 if the element is part of the membrane!
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                                typename Traits::RangeType& y) const
    {
      y = 0.0;
      if(not physics.isMembrane(e))
      {
        BaseT::evaluate(e,x,y);
      }
    }

private:
  PHYSICS& physics;
};


#endif /* DUNE_AX1_DISCRETEGRIDFUNCTION_HH */
