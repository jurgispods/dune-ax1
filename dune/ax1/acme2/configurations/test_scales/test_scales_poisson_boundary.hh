/*
 * default_boundary.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2_TEST_SCALES_POISSON_BOUNDARY_HH
#define DUNE_AX1_ACME2_TEST_SCALES_POISSON_BOUNDARY_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/acme2/common/acme2_parametertree.hh>

template<typename GV, typename RF, typename GF_INITIAL_POT>
class TestScalesPoissonBoundary
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    TestScalesPoissonBoundary(const Acme2Parameters& params_, const GF_INITIAL_POT& gfInitialPot_)
    : params(params_),
      gfInitialPot(gfInitialPot_)
    {}

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, bool isMembraneInterface) const
    {
      if(isMembraneInterface)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::None;
      }

      //debug_jochen << "[P boundary] @ " << is.geometry().global(x) << ": ";

      //debug_jochen << (is.geometry().global(x)[1]+1e-6) << "<->" << params.yMax() << std::endl;
      bool top = is.geometry().global(x)[1]+1e-6 > params.yMax();
      bool bottom = is.geometry().global(x)[1]-1e-6 < params.yMin();
      bool left = is.geometry().global(x)[0]-1e-6 < params.xMin();
      bool right = is.geometry().global(x)[0]+1e-6 > params.xMax();

      // Dirichlet-0 on upper boundary (extracellular space)
      if(top)
      {
        //return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
        //debug_jochen << "Dirichlet" << std::endl;
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
      // Inside cell or extracellular side boundaries
      } else {
        // After equilibration: use Dirichlet for the side boundaries
        if(false /*(not params.doEquilibration() || time > params.tEquilibrium())
            && (left || right) */)
        {
          //debug_jochen << "Dirichlet" << std::endl;
          return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
        } else {
          //debug_jochen << "Neumann" << std::endl;
          return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
        }
      }
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, double time) const
    {
      typename Traits::DomainType xglobal = e.geometry().global(x);
      debug_jochen << e.geometry().center();

      // Return 0 at upper boundary
      if(xglobal[1]+1e-6 > params.yMax())
      {
        return 0.0;
      // Return equilibrium value at side boundaries
      } else {
        typename GF_INITIAL_POT::Traits::RangeType y;
        debug_jochen << ": g= " << y << std::endl;
        gfInitialPot.evaluate(e,x,y);
        return y;
      }

      return 0.0;
      /*
      if(xglobal > 0)
      {
        return 0.0;
      } else {
        return -2.7; // approx. -65mV
      }
      */
      //typename Traits::RangeFieldType norm = xglobal.two_norm2();

      //typename INITIAL_POT::Traits::RangeType pot = 0.0;
      //initialPot.evaluate(e,x,pot);
    }

    //! Neumann boundary condition
    typename Traits::RangeFieldType
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, double time) const
    {
      return 0.0;
    }

    //! outflow boundary condition
    typename Traits::RangeFieldType
    o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, double time) const
    {
      return 0.0;
    }

  private:
    const Acme2Parameters& params;
    const GF_INITIAL_POT& gfInitialPot;
};

#endif /* DUNE_AX1_ACME2_TEST_SCALES_POISSON_BOUNDARY_HH */
