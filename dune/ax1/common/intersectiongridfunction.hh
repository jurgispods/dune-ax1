/*
 * intersectiongridfunction.hh
 *
 *  Created on: Jun 10, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_INTERSECTIONGRIDFUNCTION_HH
#define DUNE_AX1_INTERSECTIONGRIDFUNCTION_HH

#include <dune/pdelab/common/function.hh>

namespace Dune
{
  namespace PDELab
  {


    //! \brief traits class holding function signature, same as in local function
    //! \tparam GV The type of the grid view the function lives on.
    //! \tparam RF The numeric type of the field representing the range.
    //! \tparam m The dimension of the range.
    //! \tparam R The type of the range.
    template<class GV, class RF, int m, class R>
    struct IntersectionGridFunctionTraits
      : public FunctionTraits<typename GV::Grid::ctype, GV::dimension-1,
                    Dune::FieldVector<typename GV::Grid::ctype,
                                                  GV::dimension-1>,
                  RF, m, R>
    {
      //! \brief Export grid view type in addition
      typedef GV GridViewType;
    };


    //! \brief A BoundaryGridFunction allows evaluation on boundary intersections
    // \tparam T The type of the BoundaryGridFunctionTraits.
    // \tparam Imp The type of the implementing class.
    template<class T, class Imp>
    class IntersectionGridFunctionInterface
    {
      public:
        //! \brief Export type traits of the boundary grid function.
        typedef T Traits;

        /** \brief Evaluate the GridFunction at given position

         Evaluates components of the grid function at the given position and
         returns these values in a vector.

             \param[in]  ig geometry of intersection with boundary
             \param[in]  x The position in entity-local coordinates
             \param[out] y The result of the evaluation
        */
         template<typename I>
        inline void evaluate (const I& ig,
                 const typename Traits::DomainType& x,
                 typename Traits::RangeType& y) const
        {
           asImp().evaluate(ig,x,y);
        }

         //! get a reference to the GridView
        inline const typename Traits::GridViewType& getGridView () const
        {
          return asImp().getGridView();
        }

      private:
        Imp& asImp () {return static_cast<Imp &> (*this);}
        const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };


    /** \brief leaf of a function tree
     *
     *  Classes derived from this class implement a \ref GridFunctionTree.
     *
     *  \tparam T   Traits class holding the functions signature
     *  \tparam Imp Class implementing the function.  Imp must be derived from
     *              GridFunctionBase in some way (Barton-Nackman-Trick).
     */
    template<class T, class Imp>
    class IntersectionGridFunctionBase
        : public IntersectionGridFunctionInterface<T,Imp>
        , public TypeTree::LeafNode
    {
    public:
        typedef GridFunctionTag ImplementationTag;
        //! Type of the GridView
      typedef typename T::GridViewType GridViewType;
    };


    //! Takes an IntersectionGridFunction and acts as a single component
    template<class T>
    class IntersectionGridFunctionSelectComponentAdapter
        : public IntersectionGridFunctionInterface<BoundaryGridFunctionTraits<typename T::Traits::GridViewType,
                                                                          typename T::Traits::RangeFieldType,1,
                                                                          Dune::FieldVector<typename T::Traits::RangeFieldType,1> > ,
                                              IntersectionGridFunctionSelectComponentAdapter<T> >
    {
      typedef IntersectionGridFunctionInterface<IntersectionGridFunctionTraits<typename T::Traits::GridViewType,
                                                                         typename T::Traits::RangeFieldType,1,
                                                                         Dune::FieldVector<typename T::Traits::RangeFieldType,1> > ,
                                              IntersectionGridFunctionSelectComponentAdapter<T> > BaseT;
      public:
        //! \brief Export type traits
        typedef typename BaseT::Traits Traits;

        IntersectionGridFunctionSelectComponentAdapter (const T& t_, int k_) : t(t_), k(k_) {}

        /** \brief Evaluate all basis function at given position

          Evaluates all shape functions at the given position and returns
          these values in a vector.
        */
        template<typename I>
        inline void evaluate (const I& is,
                  const typename Traits::DomainType& x,
                  typename Traits::RangeType& y) const
        {
          typename T::Traits::RangeType Y;
          t.evaluate(is,x,Y);
          y = Y[k];
        }

          //! get a reference to the GridView
        inline const typename Traits::GridViewType& getGridView () const
        {
          return t.getGridView();
        }


        //! set component to be selected
        void select (int k_)
        {
          k = k_;
        }

      private:
        const T& t;
        int k;
    };


  }
}



#endif /* DUNE_AX1_INTERSECTIONGRIDFUNCTION_HH */
