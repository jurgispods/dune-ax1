/*
 * default_boundary.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_ES_POISSON_BOUNDARY_HH
#define DUNE_AX1_ACME2CYL_ES_POISSON_BOUNDARY_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>

template<typename GV, typename RF, typename GF_INITIAL_POT>
class ESPoissonBoundary
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    enum { timeDependent = true};

    ESPoissonBoundary(const Acme2CylParameters& params_, const GF_INITIAL_POT& gfInitialPot_)
    : params(params_),
      gfInitialPot(gfInitialPot_),
      topDirichlet(params_.isTopBoundaryDirichlet_Potential()),
      bottomDirichlet(params_.isBottomBoundaryDirichlet_Potential()),
      leftCytosolDirichlet(params_.isLeftCytosolBoundaryDirichlet_Potential()),
      rightCytosolDirichlet(params_.isRightCytosolBoundaryDirichlet_Potential()),
      leftExtracellularDirichlet(params_.isLeftExtracellularBoundaryDirichlet_Potential()),
      rightExtracellularDirichlet(params_.isRightExtracellularBoundaryDirichlet_Potential()),
      membraneNeumann(params_.isMembraneBoundaryNeumann_Potential())
    {}

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, bool isMembraneInterface) const
    {
      BCType bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;

      if(isMembraneInterface)
      {
        if(membraneNeumann)
          return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
        else
          return Dune::PDELab::ConvectionDiffusionBoundaryConditions::None;
      }

      bool top = is.geometry().global(x)[1]+1e-6 > params.yMax();
      bool bottom = is.geometry().global(x)[1]-1e-6 < params.yMin();
      bool left = is.geometry().global(x)[0]-1e-6 < params.xMin();
      bool right = is.geometry().global(x)[0]+1e-6 > params.xMax();

      // TOP
      if(top)
      {
        bctype = (topDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                               : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        return bctype;
      }
      // BOTTOM
      if(bottom)
      {
        bctype = (bottomDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                  : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        return bctype;
      }
      // LEFT (do not distinguish between cytosol/extracellular boundaries for now)
      if(left)
      {
        bctype = (leftCytosolDirichlet && leftExtracellularDirichlet
            ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
            : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        return bctype;
      }
      // RIGHT (do not distinguish between cytosol/extracellular boundaries for now)
      if(right)
      {
        bctype = (rightCytosolDirichlet && rightExtracellularDirichlet
            ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
            : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        return bctype;
      }

      debug_jochen << "x = " << is.geometry().global(x) << " => bctype = " << bctype << std::endl;
      //return bctype;
      DUNE_THROW(Dune::Exception, "Error checking this boundary's position!");
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
        typename GF_INITIAL_POT::Traits::RangeType y;
        debug_jochen << ": g= " << y << std::endl;
        gfInitialPot.evaluate(e,x,y);
        return y;
      }

      // Return equilibrium value at side boundaries, deactivated in bctype method!
      if(left || right)
      {
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

    void setTime (double t)
    {
      debug_jochen << "set time to " << t << std::endl;
      gfInitialPot.setTime(t);
    }

  private:
    const Acme2CylParameters& params;
    const GF_INITIAL_POT& gfInitialPot;

    const bool topDirichlet;
    const bool bottomDirichlet;
    const bool leftCytosolDirichlet;
    const bool rightCytosolDirichlet;
    const bool leftExtracellularDirichlet;
    const bool rightExtracellularDirichlet;
    const bool membraneNeumann;
};

#endif /* DUNE_AX1_ACME2CYL_DEFAULT_POISSON_BOUNDARY_HH */
