/*
 * default_boundary.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_STEP_POISSON_BOUNDARY_HH
#define DUNE_AX1_ACME2CYL_STEP_POISSON_BOUNDARY_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

template<typename GV, typename RF, typename GF_INITIAL_POT>
class StepPoissonBoundary
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    StepPoissonBoundary(const Acme2CylParameters& params_, const GF_INITIAL_POT& gfInitialPot_)
    : params(params_)
    {}

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, bool isMembraneInterface) const
    {
      //debug_verb << "BOUNDARY @ x = " << is.geometry().global(x) << std::endl;
      //return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;

      bool top = is.geometry().global(x)[1]+1e-6 > params.yMax();
      bool bottom = is.geometry().global(x)[1]-1e-6 < params.yMin();
      bool left = is.geometry().global(x)[0]-1e-6 < params.xMin();
      bool right = is.geometry().global(x)[0]+1e-6 > params.xMax();

//      bool corner = (top && left) || (top && right) || (bottom && left) || (bottom && right);
//
//      for(int i=0; i<is.geometry().corners(); i++)
//      {
//        bool isCornerTop = is.geometry().corner(i)[1]+1e-6 > params.yMax();
//        bool isCornerBottom = is.geometry().corner(i)[1]-1e-6 < params.yMin();
//        bool isCornerLeft = is.geometry().corner(i)[0]-1e-6 < params.xMin();
//        bool isCornerRight = is.geometry().global(i)[0]+1e-6 > params.xMax();
//        corner = corner ||
//            (isCornerTop && isCornerLeft) || (isCornerTop && isCornerRight)
//            || (isCornerBottom && isCornerLeft) || (isCornerBottom && isCornerRight);
//      }
//
//      if(corner)
//      {
//        debug_jochen << "CORNER!! @ " << is.geometry().global(x) << std::endl;
//        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
//      }

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
    const Acme2CylParameters& params;
};

#endif /* DUNE_AX1_ACME2CYL_STEP_POISSON_BOUNDARY_HH */
