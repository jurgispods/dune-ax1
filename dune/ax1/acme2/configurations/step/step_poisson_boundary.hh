/*
 * default_boundary.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2_STEP_POISSON_BOUNDARY_HH
#define DUNE_AX1_ACME2_STEP_POISSON_BOUNDARY_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

template<typename GV, typename RF, typename GF_INITIAL_POT>
class StepPoissonBoundary
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    StepPoissonBoundary(const Acme2Parameters& params_, const GF_INITIAL_POT& gfInitialPot_)
    : params(params_)
    {}

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, bool isMembraneInterface) const
    {
      //debug_verb << "BOUNDARY @ x = " << is.geometry().global(x) << std::endl;
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, double time) const
    {
      //typename Traits::DomainType xglobal = e.geometry().global(x);
      //typename Traits::RangeFieldType norm = xglobal.two_norm2();

      //typename INITIAL_POT::Traits::RangeType pot = 0.0;
      //initialPot.evaluate(e,x,pot);

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
};

#endif /* DUNE_AX1_ACME2_STEP_POISSON_BOUNDARY_HH */
