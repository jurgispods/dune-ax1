/*
 * hamburger.hh
 *
 *  Created on: Oct 13, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_HAMBURGER_CONFIG_HH
#define DUNE_AX1_HAMBURGER_CONFIG_HH

#include <dune/pdelab/common/function.hh>

#include <dune/ax1/common/electrolyte.hh>
#include <dune/ax1/acme1/configurations/hamburger/hamburger_solution_con.hh>
#include <dune/ax1/acme1/configurations/hamburger/hamburger_solution_pot.hh>
#include <dune/ax1/acme1/configurations/hamburger/hamburger_initial.hh>
#include <dune/ax1/acme1/configurations/hamburger/hamburger_nernst_planck_boundary.hh>
#include <dune/ax1/acme1/configurations/hamburger/hamburger_poisson_boundary.hh>


template<typename T>
class HamburgerConfiguration
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
      typedef HamburgerCon<GV,RF,NUMBER_OF_SPECIES>     SOLUTION_CON;
      typedef HamburgerPot<GV,RF,1>                     SOLUTION_POT;
      typedef HamburgerInitial<GV,RF,NUMBER_OF_SPECIES,SOLUTION_CON> INITIAL_CON;
      typedef HamburgerNernstPlanckBoundary<SubGV,RF,INITIAL_CON> NERNST_PLANCK_BOUNDARY;
      typedef HamburgerPoissonBoundary<GV,RF> POISSON_BOUNDARY;
    };

    HamburgerConfiguration()
    : name("hamburger"),
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

#endif /* DUNE_AX1_HAMBURGER_CONFIG_HH */
