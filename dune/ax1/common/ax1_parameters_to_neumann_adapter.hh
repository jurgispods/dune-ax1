/*
 * ax1_parameters_to_neumann_adapter.hh
 *
 *  Created on: Oct 13, 2014
 *      Author: jpods
 */

#ifndef DUNE_AX1_PARAMETERS_TO_NEUMANN_ADAPTER_HH
#define DUNE_AX1_PARAMETERS_TO_NEUMANN_ADAPTER_HH

template<typename PARAM, typename J_CON>
class Ax1ParametersToNeumannAdapter
  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::IntersectionGridFunctionTraits<typename PARAM::Traits::GridViewType,
                                       typename PARAM::Traits::RangeFieldType,
                                       1,
                                       Dune::FieldVector<typename PARAM::Traits::RangeFieldType,1> >,
                                       Ax1ParametersToNeumannAdapter<PARAM,J_CON> >
{
public:
  typedef Dune::PDELab::IntersectionGridFunctionTraits<typename PARAM::Traits::GridViewType,   // grid view type
                                           typename PARAM::Traits::RangeFieldType, // image field type (double)
                                           1,  // number of components of image
                                           Dune::FieldVector<typename PARAM::Traits::RangeFieldType,1> // image type (Dune::FieldVector<double, 1>)
                                           > Traits;

  Ax1ParametersToNeumannAdapter(const PARAM& param_, J_CON& jCon_)
  : param(param_),
    jCon(jCon_)
  {
  }



  template<typename I>
  inline void evaluate (const I &is,
                          const typename Traits::DomainType &x,
                          typename Traits::RangeType &y) const
  {
    typename J_CON::Traits::RangeType yCon;
    jCon.evaluate(is, x, yCon);

    std::vector<typename J_CON::Traits::RangeFieldType> yConVec(yCon.size(), 0.0);
    for(int j=0; j<yCon.size(); ++j)
    {
      yConVec[j] = yCon[j];
    }

    y = param.j(is, x, yConVec);
  }


  inline const typename Traits::GridViewType& getGridView () const
  {
    return param.getGridView();
  }

private:
  const PARAM& param;
  J_CON& jCon;
};

#endif /* DUNE_AX1_PARAMETERS_TO_NEUMANN_ADAPTER_HH */
