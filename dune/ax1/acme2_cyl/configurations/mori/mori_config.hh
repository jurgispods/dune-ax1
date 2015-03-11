/*
 * default.hh
 *
 *  Created on: Oct 13, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_MORI_CONFIG_HH
#define DUNE_AX1_ACME2CYL_MORI_CONFIG_HH

#include <dune/common/static_assert.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/electrolyte.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_initial.hh>
#include <dune/ax1/acme2_cyl/configurations/custom_potential.hh>
#include <dune/ax1/acme2_cyl/configurations/no_analytical_solution.hh>
#include <dune/ax1/acme2_cyl/configurations/zero_potential.hh>
#include <dune/ax1/acme2_cyl/configurations/mori/mori_poisson_parameters.hh>
#include <dune/ax1/acme2_cyl/configurations/mori/mori_nernst_planck_parameters.hh>

#include <dune/ax1/acme2_cyl/operator/acme2_cyl_mori_operator.hh>


template<typename T>
class MoriConfiguration
{
  public:
    static const bool hasAnalyticalSolution = false;
    static const bool useLogScaling = false;
    static const bool timeDependentBoundaryValues = false;

    static const int numberOfSpecies = 3;

    //! \brief global length scale scale [nm]
    static const double LENGTH_SCALE; // nano meters
    //! \brief global length time scale [ns]
    static const double TIME_SCALE; // micro seconds

    // diffusion constants
    static const double con_diffWaterNa; // m2 / s
    static const double con_diffWaterK; // m2 / s
    static const double con_diffWaterCl; // m2 / s

    static const double reduction;
    static const double absLimit;

    template<typename GV, typename RF, typename SubGV=GV>
    struct Traits
    {
      typedef InitialVoid<GV,RF,NUMBER_OF_SPECIES> INITIAL_CON;
      //typedef ZeroPotential<GV,RF,1> INITIAL_POT;
      typedef CustomBoundaryPotential<GV,RF,1> INITIAL_POT;
      typedef NoAnalyticalSolution<GV,RF,NUMBER_OF_SPECIES> SOLUTION_CON;
      typedef NoAnalyticalSolution<GV,RF,1> SOLUTION_POT;

      template<typename PHYSICS, typename GF_MEMB_FLUX, typename GF_MORI_FLUX>
      using NERNST_PLANCK_PARAMETERS = MoriNernstPlanckParameters<GV,RF,PHYSICS,INITIAL_CON,GF_MEMB_FLUX,GF_MORI_FLUX>;

      template<typename PHYSICS>
      using POISSON_PARAMETERS = MoriPoissonParameters<GV,RF,PHYSICS>;

      template<typename PARAM_CON, typename PARAM_POT, typename FEM_CON, typename FEM_POT, bool useMembContributions>
      using ELEC_OPERATOR = Acme2CylMoriOperator<PARAM_CON,PARAM_POT,FEM_CON,FEM_POT,useMembContributions>;

      static const bool useMembraneContributions = false;
      static const bool useImplicitMembraneFlux = true;
    };


    MoriConfiguration(const Acme2CylParameters& params)
    : name("mori"),
      electro(80.0, 279.45, con_mol, LENGTH_SCALE)
    {
      assert(NUMBER_OF_SPECIES == numberOfSpecies);

      // electrolyte definition ###########################

      Ion<T>    sodium(  1.0, "na" );
      electro.addIon(sodium);
      electro.setDiffConst(Na, con_diffWaterNa * TIME_SCALE / ( LENGTH_SCALE * LENGTH_SCALE ) );

      if(numberOfSpecies > 1)
      {
        Ion<T> potassium(  1.0, "k"  );
        electro.addIon(potassium);
        electro.setDiffConst(K, con_diffWaterK * TIME_SCALE / ( LENGTH_SCALE * LENGTH_SCALE ) );
      }
      if(numberOfSpecies > 2)
      {
        Ion<T>  chloride( -1.0, "cl" );
        electro.addIon(chloride);
        electro.setDiffConst(Cl, con_diffWaterCl * TIME_SCALE / ( LENGTH_SCALE * LENGTH_SCALE ) );
      }
    }

    Electrolyte<T> getElectrolyte()
    {
      return electro;
    }

    std::string& getConfigName()
    {
      return name;
    }

  private:
    Electrolyte<T> electro;
    std::string name;

};


template<typename T>
const double MoriConfiguration<T>::LENGTH_SCALE = 1.0e-9; // nano meters

template<typename T>
const double MoriConfiguration<T>::TIME_SCALE   = 1.0e-6; // micro seconds

template<typename T>
const double MoriConfiguration<T>::con_diffWaterNa = 1.33e-9; // m2 / s
template<typename T>
const double MoriConfiguration<T>::con_diffWaterK  = 1.96e-9; // m2 / s
template<typename T>
const double MoriConfiguration<T>::con_diffWaterCl = 2.03e-9; // m2 / s

template<typename T>
const double MoriConfiguration<T>::reduction = 1e-5;
template<typename T>
const double MoriConfiguration<T>::absLimit = 1e-5;


#endif /* DUNE_AX1_ACME2CYL_MORI_CONFIG_HH */
