/*
 * default_boundary.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2_DEFAULT_POISSON_BOUNDARY_HH
#define DUNE_AX1_ACME2_DEFAULT_POISSON_BOUNDARY_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/acme2/common/acme2_parametertree.hh>

template<typename GV, typename RF, typename GF_INITIAL_POT>
class DefaultPoissonBoundary
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    DefaultPoissonBoundary(const Acme2Parameters& params_, const GF_INITIAL_POT& gfInitialPot_)
    : params(params_),
      gfInitialPot(gfInitialPot_)
    {}

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, bool isMembraneInterface) const
    {
      // Fix potential at all boundaries when there is no membrane
      if(not params.useMembrane())
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
      }


      if(isMembraneInterface)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::None;
      }

      //debug_verb << "BOUNDARY @ x = " << is.geometry().global(x) << std::endl;

      //debug_jochen << (is.geometry().global(x)[1]+1e-6) << "<->" << params.yMax() << std::endl;
      bool top = is.geometry().global(x)[1]+1e-6 > params.yMax();
      bool bottom = is.geometry().global(x)[1]-1e-6 < params.yMin();
      bool left = is.geometry().global(x)[0]-1e-6 < params.xMin();
      bool right = is.geometry().global(x)[0]+1e-6 > params.xMax();

      // Dirichlet-0 on upper boundary (extracellular space)
      if(top)
      {
        //return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
      }

      // Dirichlet-0 on lower boundary (if it is extracellular space, like in the 2-membrane case)
      if(params.nMembranes() == 2 && bottom)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
      }

      // After equilibration: use Dirichlet for the side boundaries
      // EDIT: Deactivated Dirichlet BC for sides, as it is not necessary and hinders the AP to evolve!
      if(false /*(not params.doEquilibration() || time > params.tEquilibrium())
          && (left || right) */)
      {
        //debug_jochen << is.geometry().center() << ": Dirichlet" << std::endl;
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
      }

      //debug_jochen << is.geometry().center() << ": Neumann" << std::endl;
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;


    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, double time) const
    {
      //typename Traits::DomainType xglobal = e.geometry().global(x);
      bool top = e.geometry().global(x)[1]+1e-6 > params.yMax();
      bool bottom = e.geometry().global(x)[1]-1e-6 < params.yMin();
      bool left = e.geometry().global(x)[0]-1e-6 < params.xMin();
      bool right = e.geometry().global(x)[0]+1e-6 > params.xMax();

      //debug_jochen << e.geometry().center();

      // Return 0 at upper boundary
      if(top)
      {
        return 0.0;
      }
      if(bottom)
      {
        return 0.0;
      }

      // Return equilibrium value at side boundaries, deactivated in bctype method!
      if(left || right)
      {
        typename GF_INITIAL_POT::Traits::RangeType y;
        debug_jochen << ": g= " << y << std::endl;
        gfInitialPot.evaluate(e,x,y);
        return y;
      }

      return 0.0;
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

#endif /* DUNE_AX1_ACME2_DEFAULT_POISSON_BOUNDARY_HH */
