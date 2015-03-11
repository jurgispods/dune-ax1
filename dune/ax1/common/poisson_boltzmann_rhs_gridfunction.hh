/*
 * membranepotentialgridfunction.hh
 *
 *  Created on: Feb 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_POISSONBOLTZMANNRHSGRIDFUNCTION_HH
#define DUNE_AX1_POISSONBOLTZMANNRHSGRIDFUNCTION_HH

#include <dune/ax1/common/constants.hh>


template<typename DGF_POT, typename GF_INITIAL_CON, typename PHYSICS>
class PoissonBoltzmannRHSGridFunction
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<typename DGF_POT::Traits::GridViewType,
                                       typename DGF_POT::Traits::RangeFieldType,
                                       1,
                                       Dune::FieldVector<typename DGF_POT::Traits::RangeFieldType,1> >,
     PoissonBoltzmannRHSGridFunction<DGF_POT, GF_INITIAL_CON, PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename DGF_POT::Traits::GridViewType,   // grid view type
                                           typename DGF_POT::Traits::RangeFieldType, // image field type (double)
                                           1,                                        // number of components of image (1)
                                           Dune::FieldVector<typename DGF_POT::Traits::RangeFieldType,1> // image type (Dune::FieldVector<double, 1>)
                                           > Traits;

  //! constructor
  PoissonBoltzmannRHSGridFunction (DGF_POT& dgfPot_, GF_INITIAL_CON& gfInitialCon_, PHYSICS& physics_)
    : dgfPot(dgfPot_), gfInitialCon(gfInitialCon_), gridView(dgfPot_.getGridView()), physics(physics_),
      pot(0.0)
  {
    /*
    typename Traits::RangeType rt;
    typename Traits::RangeFieldType rf;
    typename Traits::DomainType dt;
    typename Traits::DomainFieldType df;
    typename DGF_CON::Traits::RangeType dgfConRT;

    debug_verb << "===== ChargeDensityGridFunction TYPEINFO ====" << std::endl;
    debug_verb << "Traits::RangeType = " << Tools::getTypeName(rt) << std::endl;
    debug_verb << "Traits::RangeFieldType = " << Tools::getTypeName(rf) << std::endl;
    debug_verb << "Traits::DomainType = " << Tools::getTypeName(dt) << std::endl;
    debug_verb << "Traits::DomainFieldType = " << Tools::getTypeName(df) << std::endl;
    debug_verb << "DGF_CON::Traits::RangeType = " << Tools::getTypeName(dgfConRT) << std::endl<< std::endl;
    */
  }

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    if(not physics.isMembrane(e))
    {
      //typename DGF_POT::Traits::RangeType pot;
      //dgfPot.evaluate(e,x,pot);

      typename GF_INITIAL_CON::Traits::RangeType initConc;
      gfInitialCon.evaluate(e, x, initConc);

      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        int valence = physics.getValence(j);
        //debug_verb << "init charge = " <<  (valence * initConc[j]) << ", exp-factor = " << std::exp(- valence * pot) << std::endl;
        y += valence * initConc[j] * std::exp(- valence * pot);
      }
      //debug_jochen << "pot: " << std::scientific << std::setprecision(16)<< pot << std::endl;
    }
  }

  inline void setPot(typename Traits::RangeFieldType& pot_)
  {
    pot = pot_;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  DGF_POT& dgfPot;
  GF_INITIAL_CON& gfInitialCon;
  const typename Traits::GridViewType& gridView;
  PHYSICS& physics;
  typename Traits::RangeFieldType pot;
};

#endif /* DUNE_AX1_POISSONBOLTZMANNRHSGRIDFUNCTION_HH */
