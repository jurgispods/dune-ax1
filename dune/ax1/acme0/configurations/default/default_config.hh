/*
 * default.hh
 *
 *  Created on: Oct 13, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME0_DEFAULT_CONFIG_HH
#define DUNE_AX1_ACME0_DEFAULT_CONFIG_HH

#include <dune/common/static_assert.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/electrolyte.hh>
#include <dune/ax1/acme0/common/acme0_initial.hh>
#include <dune/ax1/acme0/configurations/no_analytical_solution.hh>
#include <dune/ax1/acme0/configurations/default/default_poisson_boundary.hh>
#include <dune/ax1/acme0/configurations/default/default_nernst_planck_boundary.hh>


template<typename T>
class DefaultConfiguration
{
  public:
    //! \brief global length scale scale [nm]
    static constexpr double LENGTH_SCALE = 1.0e-9; // nano meters
    //! \brief global length time scale [ns]
    static constexpr double TIME_SCALE   = 1.0e-9; // nano seconds

    // diffusion constants
    static constexpr double con_diffWaterNa = 1.33e-9; // m2 / s
    static constexpr double con_diffWaterK  = 1.96e-9; // m2 / s
    static constexpr double con_diffWaterCl = 2.03e-9; // m2 / s

    static const bool hasAnalyticalSolution = false;
    static const bool useLogScaling = false;

    static constexpr double reduction = 1e-6;
    static constexpr double absLimit = 1e-10;

    static const int numberOfSpecies = 3;

    template<typename GV, typename RF>
    struct Traits
    {
      typedef InitialCon<GV,RF,NUMBER_OF_SPECIES> INITIAL_CON;
      typedef NoAnalyticalSolution<GV,RF,NUMBER_OF_SPECIES> SOLUTION_CON;
      typedef NoAnalyticalSolution<GV,RF,1> SOLUTION_POT;
      typedef DefaultNernstPlanckBoundary<GV,RF,INITIAL_CON> NERNST_PLANCK_BOUNDARY;
      typedef DefaultPoissonBoundary<GV,RF> POISSON_BOUNDARY;
    };


    DefaultConfiguration()
    : name("default"),
      electro(Solvent<T>( 80.0 ), 279.45, con_mol, LENGTH_SCALE)
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


#endif /* DUNE_AX1_ACME0_DEFAULT_CONFIG_HH */
