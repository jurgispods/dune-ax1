/*
 * ax1_processorgridfunction.hh
 *
 *  Created on: Apr 25, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_PROCESSORGRIDFUNCTION_HH
#define DUNE_AX1_PROCESSORGRIDFUNCTION_HH


template<typename GV>
class Ax1ProcessorGridFunction
  :  public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV, int, 1,
                                                                            Dune::FieldVector<int,1> >,
                                                                            Ax1ProcessorGridFunction<GV> >
{
  typedef Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV, int, 1,
                                         Dune::FieldVector<int,1> >,
                                         Ax1ProcessorGridFunction<GV> > BaseT;

public:
    typedef typename BaseT::Traits Traits;

    Ax1ProcessorGridFunction(const GV& gv_)
    : gv(gv_)
    {}


    // Evalute grid function; return 0 if the element is part of the membrane!
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                                typename Traits::RangeType& y) const
    {
      y = gv.comm().rank();
    }

    const typename Traits::GridViewType& getGridView() const
    {
      return gv;
    }

private:
    const GV& gv;
};


#endif /* DUNE_AX1_PROCESSORGRIDFUNCTION_HH */
