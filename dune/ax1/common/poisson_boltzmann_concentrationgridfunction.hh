/*
 * poisson_boltzmann_concentrationgridfunction.hh
 *
 *  Created on: Feb 24, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_POISSON_BOLTZMANN_CONCENTRATIONGRIDFUNCTION_HH
#define DUNE_AX1_POISSON_BOLTZMANN_CONCENTRATIONGRIDFUNCTION_HH

#include <dune/ax1/common/constants.hh>


template<typename INITIAL_CON, typename DGF_POT, typename SubGV, typename PHYSICS>
class PoissonBoltzmannConcentrationGridFunction
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<SubGV,
                                                                           typename INITIAL_CON::Traits::RangeFieldType,
                                                                           NUMBER_OF_SPECIES,
                                                                           typename INITIAL_CON::Traits::RangeType >,
           PoissonBoltzmannConcentrationGridFunction<INITIAL_CON,DGF_POT,SubGV,PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<SubGV,   // grid view type
                                           typename INITIAL_CON::Traits::RangeFieldType,  // range field type (double)
                                           NUMBER_OF_SPECIES,                             // # of components
                                           // Dune::FieldVector<Dune::FieldMatrix<double, dim, dim>, NUMBER_OF_SPECIES>
                                           typename INITIAL_CON::Traits::RangeType > Traits;

  //! constructor
  PoissonBoltzmannConcentrationGridFunction (INITIAL_CON& initialCon_, DGF_POT& dgfPot_, const SubGV& gridView_, PHYSICS& physics_)
    : initialCon(initialCon_),
      dgfPot(dgfPot_),
      gridView(gridView_),
      physics(physics_)
  {
    /*
    typename Traits::RangeType rt;
    typename Traits::RangeFieldType rf;
    typename Traits::DomainType dt;
    typename Traits::DomainFieldType df;
    typename INITIAL_CON_GRAD::Traits::RangeType initialConGradRT;

    debug_verb << "===== IonFluxGridFunction TYPEINFO ====" << std::endl;
    debug_verb << "Traits::RangeType = " << Tools::getTypeName(rt) << std::endl;
    debug_verb << "Traits::RangeFieldType = " << Tools::getTypeName(rf) << std::endl;
    debug_verb << "Traits::DomainType = " << Tools::getTypeName(dt) << std::endl;
    debug_verb << "Traits::DomainFieldType = " << Tools::getTypeName(df) << std::endl;
    debug_verb << "INITIAL_CON_GRAD::Traits::RangeType = " << Tools::getTypeName(initialConGradRT) << std::endl<< std::endl;
    */
  }

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    //int elemIndex = physics.getElementIndex(e);
    y = 0.0;

    if(not physics.isMembrane(e))
    {
      typename INITIAL_CON::Traits::RangeType conc0;
      typename DGF_POT::Traits::RangeType pot;

      // initialCon, dgfPot are grid functions living on the HOST grid! => Evaluate on host element
      const typename PHYSICS::Element& hostEntity = gridView.getHostElement(e);
      initialCon.evaluate(hostEntity,x,conc0);
      dgfPot.evaluate(hostEntity,x,pot);
      //pot = physics.convertTo_mV(pot);
      //debug_verb << "x=" << e.geometry().global(x) << " -> pot=" << pot << std::endl;

      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        int valence = physics.getValence(j);

        y[j] = conc0[j] * std::exp(- valence * pot);
      }
    }
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  INITIAL_CON& initialCon;
  DGF_POT& dgfPot;
  const typename Traits::GridViewType& gridView;
  PHYSICS& physics;
};


#endif /* DUNE_AX1_POISSON_BOLTZMANN_CONCENTRATIONGRIDFUNCTION_HH */
