/*
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_STEP_NERNST_PLANCK_BOUNDARY_HH
#define DUNE_AX1_ACME2CYL_STEP_NERNST_PLANCK_BOUNDARY_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

template<typename GV, typename RF, typename GF_INITIAL_CON>
class StepNernstPlanckBoundary
{


  public:
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    typedef GF_INITIAL_CON GFInitialCon;

    StepNernstPlanckBoundary(const Acme2CylParameters& params_, const GF_INITIAL_CON& gfInitialCon_)
    : params(params_),
      gfInitialCon(gfInitialCon_)
    {
    }

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, int ionSpecies, bool isMembraneInterface) const
    {
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        double time, int ionSpecies) const
    {
      //typename Traits::DomainType xglobal = e.geometry().global(x);
      //typename Traits::RangeFieldType norm = xglobal.two_norm2();
      return 0.0;
    }

    //! Neumann boundary condition
    typename Traits::RangeFieldType
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, int ionSpecies) const
    {
      return 0.0;
    }

    //! outflow boundary condition
    typename Traits::RangeFieldType
    o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, int ionSpecies) const
    {
      return 0.0;
    }

  private:
    const Acme2CylParameters& params;
    const GF_INITIAL_CON& gfInitialCon;
};

#endif /* DUNE_AX1_ACME2CYL_STEP_NERNST_PLANCK_BOUNDARY_HH */
