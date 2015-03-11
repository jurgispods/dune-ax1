/*
 * default_nernst_planck_boundary.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME0_DEFAULT_NERNST_PLANCK_BOUNDARY_HH
#define DUNE_AX1_ACME0_DEFAULT_NERNST_PLANCK_BOUNDARY_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>

template<typename GV, typename RF, typename GF_INITIAL_CON>
class DefaultNernstPlanckBoundary
{

  public:
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    DefaultNernstPlanckBoundary(const GF_INITIAL_CON& gfInitialCon_)
    : gfInitialCon(gfInitialCon_)
    {}

  //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, int ionSpecies) const
    {
      typename Traits::DomainType xglobal = is.geometry().global(x);

      if(xglobal > 0 && ionSpecies == K)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      } else {
        //return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      }
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        double time, int ionSpecies) const
    {
      typename Traits::DomainType xglobal = e.geometry().global(x);
      //typename Traits::RangeFieldType norm = xglobal.two_norm2();

      typename GF_INITIAL_CON::Traits::RangeType initConc(0.0);
      gfInitialCon.evaluate(e, x, initConc);

      return initConc[ionSpecies];
    }

    //! Neumann boundary condition
    typename Traits::RangeFieldType
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, int ionSpecies) const
    {


      typename Traits::DomainType xglobal = is.geometry().global(x);
      if(xglobal > 0 && ionSpecies == K)
      {
        debug_jochen << "t = " << time << ". [" << ION_NAMES[ionSpecies] << "] flux = " << std::abs(std::sin(time))
          << std::endl;
        return std::abs(std::sin(time));
      }
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
    const GF_INITIAL_CON& gfInitialCon;
};

#endif /* DUNE_AX1_ACME0_DEFAULT_NERNST_PLANCK_BOUNDARY_HH */
