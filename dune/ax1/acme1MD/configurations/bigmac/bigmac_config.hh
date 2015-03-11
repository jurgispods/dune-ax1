/*
 * bigmac.hh
 *
 *  Created on: Oct 13, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_BIGMAC_CONFIG_HH
#define DUNE_AX1_BIGMAC_CONFIG_HH

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/common/electrolyte.hh>
#include <dune/ax1/acme1MD/configurations/bigmac/bigmac_solution_con.hh>
#include <dune/ax1/acme1MD/configurations/bigmac/bigmac_solution_pot.hh>
#include <dune/ax1/acme1MD/configurations/bigmac/bigmac_initial.hh>
#include <dune/ax1/acme1MD/configurations/bigmac/bigmac_nernst_planck_boundary.hh>
#include <dune/ax1/acme1MD/configurations/bigmac/bigmac_poisson_boundary.hh>


template<typename T>
class BigmacConfiguration
{
  public:
    //! \brief global length scale scale [nm]
    static constexpr double LENGTH_SCALE = 1.0; // nano meters
    //! \brief global length time scale [ns]
    static constexpr double TIME_SCALE   = 1.0; // nano seconds

    // diffusion constants
    static constexpr double con_diffWaterNa = 1.0; // m2 / s

    static const bool hasAnalyticalSolution = true;
    static const bool useLogScaling = false;

    static constexpr double reduction = 1e-6;
    static constexpr double absLimit = 1e-10;

    static const int numberOfSpecies = 1;

    template<typename GV, typename RF, typename SubGV=GV>
    struct Traits
    {
      typedef BigmacCon<GV,RF,NUMBER_OF_SPECIES>     SOLUTION_CON;
      typedef BigmacPot<GV,RF,1>                     SOLUTION_POT;
      typedef BigmacInitialCon<GV,RF,NUMBER_OF_SPECIES,SOLUTION_CON> INITIAL_CON;
      typedef BigmacInitialPot<GV,RF,1,SOLUTION_POT> INITIAL_POT;
      typedef BigmacNernstPlanckBoundary<GV,RF,INITIAL_CON> NERNST_PLANCK_BOUNDARY;
      typedef BigmacPoissonBoundary<GV,RF> POISSON_BOUNDARY;
    };

    BigmacConfiguration()
    : name("bigmac"),
      electro(Solvent<T>( 1.0 ), temperatureToNormalizePoissonRHS(1.0), 1.0, LENGTH_SCALE)
    {
      assert(NUMBER_OF_SPECIES == numberOfSpecies);

      // electrolyte definition ###########################
      Ion<T>    sodium( 1.0, "na" );

      electro.addIon(sodium);

      electro.setDiffConst(Na, con_diffWaterNa * TIME_SCALE / ( LENGTH_SCALE * LENGTH_SCALE ) );
    }

    // rescale Temperature to give a RHS=1 for poisson equation
    T temperatureToNormalizePoissonRHS(const T & stdCon) const
    {
      return stdCon * con_e * con_e / ( con_k * con_eps0 );
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

#endif /* DUNE_AX1_BIGMAC_CONFIG_HH */
