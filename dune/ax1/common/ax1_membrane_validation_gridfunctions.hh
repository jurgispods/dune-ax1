/*
 * ax1_membrane_validation_gridfunctions.hh
 *
 *  Created on: Sep 5, 2012
 *      Author: jpods
 */

#ifndef AX1_MEMBRANE_VALIDATION_GRIDFUNCTIONS_HH_
#define AX1_MEMBRANE_VALIDATION_GRIDFUNCTIONS_HH_


#include <dune/ax1/common/constants.hh>


template<typename SubGV, typename DGF_CON, typename PHYSICS>
class Ax1NernstEquationGridFunction
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<
        SubGV,
        typename DGF_CON::Traits::RangeFieldType,
        1,
        Dune::FieldVector<typename DGF_CON::Traits::RangeFieldType,1> >,
      Ax1NernstEquationGridFunction<SubGV, DGF_CON, PHYSICS> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<SubGV,                           // grid view type
               typename DGF_CON::Traits::RangeFieldType,                     // range field type (double)
               1, // number of components of image
               Dune::FieldVector<typename DGF_CON::Traits::RangeFieldType,1>// image type (Dune::FieldVector<double, 1>)
               > Traits;

  //! constructor
  Ax1NernstEquationGridFunction (const SubGV& gv_, DGF_CON& dgfCon_, PHYSICS& physics_,
      int ionSpecies_ = 0)
    : gridView(gv_), dgfCon(dgfCon_), physics(physics_), ionSpecies(ionSpecies_)
  {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0.0;

    typename DGF_CON::Traits::RangeType conCytosol, conExtracellular;
    // Evaluate concentration jump on MD grid
    physics.getMembraneConcentrationJump(gridView.grid().multiDomainEntity(e), dgfCon, conCytosol, conExtracellular);

    typename PHYSICS::Real T = physics.getElectrolyte.getTemperature();
    int valence = physics.getValence(ionSpecies);

    // Now that we have all the information, evaluate Nernst equation!
    // E = (R*T)/(z*F) * ln(con_out/con_in) = (k_B * T) / (e*z) * ln(con_out/con_in)
    y = con_k * T / (con_e * valence) * std::log(conExtracellular[ionSpecies]/conCytosol[ionSpecies]);

    return;
   }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  const SubGV& gridView;
  DGF_CON& dgfCon;
  PHYSICS& physics;
  const int ionSpecies;
};


template<typename DGF_CON, typename PHYSICS>
class Ax1GoldmanEquationGridFunction
  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::BoundaryGridFunctionTraits<typename DGF_CON::Traits::GridViewType,
                                     typename DGF_CON::Traits::RangeFieldType,
                                     1,
                                     Dune::FieldVector<typename DGF_CON::Traits::RangeFieldType,1>
      >,
      Ax1GoldmanEquationGridFunction<DGF_CON, PHYSICS> >
{
public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<typename DGF_CON::Traits::GridViewType,   // grid view type
                                           typename DGF_CON::Traits::RangeFieldType, // image field type (double)
                                           1,  // number of components of image
                                           Dune::FieldVector<typename DGF_CON::Traits::RangeFieldType, 1> // image type (Dune::FieldVector<double, 1>)
                                           > Traits;

  typedef typename PHYSICS::ChannelSet CHANNELS;

  //! constructor
  Ax1GoldmanEquationGridFunction (DGF_CON& dgfCon_, PHYSICS& physics_)
    : gridView(dgfCon_.getGridView()), dgfCon(dgfCon_), physics(physics_)
  {
  }

  template<typename I>
  inline void evaluate (const I &is,
                          const typename Traits::DomainType &x,
                          typename Traits::RangeType &y) const
  {
    y = 0.0;

    const typename PHYSICS::ChannelSet& channels = physics.getChannelSet();
    int iIndex = physics.getIntersectionIndex(is);

    std::vector<typename Traits::RangeFieldType> leakConductances(3, 0.0);

    // Loop over channels and pick out leak channels
    for(int k=0; k<channels.size(); k++)
    {
      if(channels.getChannel(k).isLeakChannel())
      {
        // The total relative permeability is available via the 0-th gating particle of this leak channel
        leakConductances[channels.getChannel(k).getIonSpecies()] += channels.getGatingParticle(k, 0, iIndex);
      }
    }



    typename DGF_CON::Traits::RangeType conCytosol, conExtracellular;
    // Evaluate concentration jump on MD grid
    physics.getMembraneConcentrationJump(is, dgfCon, conCytosol, conExtracellular);

    typename Traits::RangeFieldType sumNumerator(0.0);
    typename Traits::RangeFieldType sumDenominator(0.0);
    typename Traits::RangeFieldType T = physics.getElectrolyte().getTemperature();

    // The maximum channel conductances are used for the permeabilities here. This is still valid as only
    // the _relative_ permeabilities matter!
    for(int j=0; j<NUMBER_OF_SPECIES; j++)
    {
      debug_verb << "leakConductance[" << ION_NAMES[j] << "] = " << leakConductances[j] << std::endl;

      int valence = physics.getValence(j);

      sumNumerator += (valence > 0 ? leakConductances[j] * conCytosol[j]
                                   : leakConductances[j] * conExtracellular[j]);
      sumDenominator += (valence > 0 ? leakConductances[j] * conExtracellular[j]
                                     : leakConductances[j] * conCytosol[j]);

    }

    // Now that we have all the information, evaluate Nernst equation! A_i:
    // E = - (k_B * T) / (e*z) * ln((sum[P_i * M_i_out]+sum[P_i * A_i_in])/(sum[P_i*M_i_in]+sum[P_i * A_i_out]))
    y = - con_k * T / (con_e) * std::log(sumNumerator/sumDenominator);

    // Bring to mV unit
    y *= 1000;

    debug_verb << "k*T/e = " << (con_k * T / (con_e)) << std::endl;
    debug_verb << "log(" << (sumNumerator/sumDenominator) << std::endl;

    return;
   }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gridView;
  }

private:
  const typename Traits::GridViewType& gridView;
  DGF_CON& dgfCon;
  PHYSICS& physics;

};

#endif /* AX1_MEMBRANE_VALIDATION_GRIDFUNCTIONS_HH_ */
