/*
 * error_norms.hh
 *
 *  Created on: Oct 26, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ERROR_NORMS_HH
#define DUNE_AX1_ERROR_NORMS_HH

#include <dune/pdelab/function/minus.hh>

/*! \brief Adapter returning ||f1(x)-f2(x)||^2 for two given grid functions

  \tparam T1  a grid function type
  \tparam T2  a grid function type
*/
template<typename T1, typename T2>
class DifferenceSquaredAdapter
  : public Dune::PDELab::GridFunctionBase<
            Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
            typename T1::Traits::RangeFieldType,
            1,
            Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >
            ,DifferenceSquaredAdapter<T1,T2> >
{
  public:
    typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                             typename T1::Traits::RangeFieldType,
                                             1,
                                             Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

    //! constructor
    DifferenceSquaredAdapter (const T1& t1_, const T2& t2_) : t1(t1_), t2(t2_) {}

    //! \copydoc GridFunctionBase::evaluate()
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                          typename Traits::RangeType& y) const
    {
      //debug_jochen << "[DifferenceSquaredAdapter::evaluate] @" << e.geometry().global(x) << std::endl;

      typename T1::Traits::RangeType y1;
      t1.evaluate(e,x,y1);
      //debug_jochen << "    y1: " << y1 << "  |";
      typename T2::Traits::RangeType y2;
      t2.evaluate(e,x,y2);
      //debug_jochen << "y2: " << y2 << std::endl;

      y1 -= y2;
      //debug_jochen << "  => diff = " << y1 << std::endl;

      y = y1.two_norm2();
    }

    inline const typename Traits::GridViewType& getGridView () const
    {
      return t1.getGridView();
    }

  private:
    const T1& t1;
    const T2& t2;
};

// FIXME Adapt this to parallel and cylinder cordinates!
class ErrorNorms
{
  public:

    template<typename GF1, typename GF2, bool useCylinderCoordinates=USE_CYLINDER_COORDINATES>
    static double l2Norm(GF1& gf1, GF2& gf2, int intorder=8)
    {
      typedef typename GF1::Traits::GridViewType GV;
      const GV& gv = gf1.getGridView();

      debug_jochen << "  Calculating L2 error using " << gv.size(0) << " elements..." << std::endl;

      const int dim = GV::dimension;

      DifferenceSquaredAdapter<GF1, GF2> diffSquared(gf1, gf2);
      //Dune::PDELab::MinusGridFunctionAdapter<GF1,GF2> diff(gf1, gf2);

      typedef typename GV::template Codim<0>::
        template Partition<Dune::Interior_Partition>::Iterator EIterator;
      typedef typename DifferenceSquaredAdapter<GF1, GF2>::Traits::DomainType DT;
      typedef typename DifferenceSquaredAdapter<GF1, GF2>::Traits::DomainFieldType DF;
      typedef typename DifferenceSquaredAdapter<GF1, GF2>::Traits::RangeType RT;
      typedef typename DifferenceSquaredAdapter<GF1, GF2>::Traits::RangeFieldType RF;

      typedef typename GV::template Codim<0>::Geometry GeometryOrig;
      // Use switch to choose original or cylinder Geometry type
      typedef typename Acme2CylGeometrySwitch::GeometrySwitch<GeometryOrig,useCylinderCoordinates>::type Geometry;

      double l2Error = 0.0;

      int count = 0;
      for (EIterator eit=gv.template begin<0,Dune::Interior_Partition>();
              eit!=gv.template end<0,Dune::Interior_Partition>(); ++eit)
      {
        const GeometryOrig& geoOrig = eit->geometry();
        const Geometry& geo(geoOrig);

        // select quadrature rule
        Dune::GeometryType gt = geo.type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        //debug_jochen << "[l2norm] Element #" << count++ << " @" << eit->geometry().center() << std::endl;
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        {
          RF weight = it->weight();
          RF slope = geo.integrationElement(it->position());

          RT error(0.0);
          //diffSquared.evaluate(*eit, it->position(), error);
          diffSquared.evaluate(*eit, it->position(), error);

//          RT err(0.0);
//          gf1.evaluate(*eit, it->position(), err);
//          debug_jochen << "gf1: " << err << std::endl;
//          gf2.evaluate(*eit, it->position(), err);
//          debug_jochen << "gf2: " << err << std::endl;
//
//          debug_jochen << "error: " << error << std::endl;
//          debug_jochen << "Adding value: " << error*error << std::endl;

          //l2Error += error.infinity_norm();
          l2Error += (error * weight * slope);
        }
      }

      l2Error = std::sqrt(l2Error);
      return l2Error;
    }

    template<typename GF1, typename GF2, bool useCylinderCoordinates=USE_CYLINDER_COORDINATES>
    static double maxNorm(GF1& gf1, GF2& gf2, int intorder=8)
    {
      typedef typename GF1::Traits::GridViewType GV;
      const GV& gv = gf1.getGridView();

      const int dim = GV::dimension;

      Dune::PDELab::MinusGridFunctionAdapter<GF1,GF2> diff(gf1, gf2);

      typedef typename GV::template Codim<0>::
        template Partition<Dune::Interior_Partition>::Iterator EIterator;
      typedef typename Dune::PDELab::MinusGridFunctionAdapter<GF1, GF2>::Traits::DomainType DT;
      typedef typename Dune::PDELab::MinusGridFunctionAdapter<GF1, GF2>::Traits::DomainFieldType DF;
      typedef typename Dune::PDELab::MinusGridFunctionAdapter<GF1, GF2>::Traits::RangeType RT;
      typedef typename Dune::PDELab::MinusGridFunctionAdapter<GF1, GF2>::Traits::RangeFieldType RF;

      DT maxPosition(0.0);

      typedef typename GV::template Codim<0>::Geometry GeometryOrig;
      // Use switch to choose original or cylinder Geometry type
      typedef typename Acme2CylGeometrySwitch::GeometrySwitch<GeometryOrig,useCylinderCoordinates>::type Geometry;

      double maxError = 0.0;
      for (EIterator eit=gv.template begin<0,Dune::Interior_Partition>();
          eit!=gv.template end<0,Dune::Interior_Partition>(); ++eit)
      {
        const GeometryOrig& geoOrig = eit->geometry();
        const Geometry& geo(geoOrig);

        // select quadrature rule
        Dune::GeometryType gt = geo.type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        {
          RF weight = it->weight();
          RF slope = geo.integrationElement(it->position());

          RT error(0.0);
          //diffSquared.evaluate(*eit, it->position(), error);
          diff.evaluate(*eit, it->position(), error);

          //l2Error += error.infinity_norm();
          if(error.infinity_norm() > maxError)
          {
            maxError = error.infinity_norm();
            maxPosition = eit->geometry().global(it->position());
          }
        }
      }
      debug_info << "  Max error encountered at position " << maxPosition << std::endl;
      return maxError;
    }
};


#endif /* DUNE_AX1_ERROR_NORMS_HH */
