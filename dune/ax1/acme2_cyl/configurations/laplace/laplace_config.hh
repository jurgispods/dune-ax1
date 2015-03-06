/*
 * default.hh
 *
 *  Created on: Oct 13, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_LAPLACE_CONFIG_HH
#define DUNE_AX1_ACME2CYL_LAPLACE_CONFIG_HH

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/electrolyte.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_initial.hh>
#include <dune/ax1/acme2_cyl/configurations/custom_potential.hh>
#include <dune/ax1/acme2_cyl/configurations/no_analytical_solution.hh>
#include <dune/ax1/acme2_cyl/configurations/zero_potential.hh>
#include <dune/ax1/acme2_cyl/configurations/ES/ES_poisson_boundary.hh>
#include <dune/ax1/acme2_cyl/configurations/ES/ES_nernst_planck_boundary.hh>
#include <dune/ax1/acme2_cyl/configurations/ES/ES_initial.hh>


template<typename T>
class LaplaceConfiguration
{
  public:
    static const bool hasAnalyticalSolution = false;
    static const bool useLogScaling = false;
    static const bool timeDependentBoundaryValues = true;

    static const int numberOfSpecies = 0;

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
      typedef ESInitialPot<GV,RF,1> INITIAL_POT;
      typedef NoAnalyticalSolution<GV,RF,NUMBER_OF_SPECIES> SOLUTION_CON;
      typedef NoAnalyticalSolution<GV,RF,1> SOLUTION_POT;
      typedef ESNernstPlanckBoundary<GV,RF,INITIAL_CON> NERNST_PLANCK_BOUNDARY;
      typedef ESPoissonBoundary<GV,RF,INITIAL_POT> POISSON_BOUNDARY;
    };


    LaplaceConfiguration()
    : name("laplace"),
      electro(80.0, 279.45, con_mol, LENGTH_SCALE)
    {
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
const double LaplaceConfiguration<T>::LENGTH_SCALE = 1.0e-9; // nano meters

template<typename T>
const double LaplaceConfiguration<T>::TIME_SCALE   = 1.0e-6; // micro seconds

template<typename T>
const double LaplaceConfiguration<T>::con_diffWaterNa = 1.33e-9; // m2 / s
template<typename T>
const double LaplaceConfiguration<T>::con_diffWaterK  = 1.96e-9; // m2 / s
template<typename T>
const double LaplaceConfiguration<T>::con_diffWaterCl = 2.03e-9; // m2 / s

template<typename T>
const double LaplaceConfiguration<T>::reduction = 1e-5;
template<typename T>
const double LaplaceConfiguration<T>::absLimit = 1e-5;



#endif /* DUNE_AX1_ACME2CYL_DEFAULT_CONFIG_HH */
