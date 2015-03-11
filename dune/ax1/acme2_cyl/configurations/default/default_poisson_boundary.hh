/*
 * default_boundary.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_DEFAULT_POISSON_BOUNDARY_HH
#define DUNE_AX1_ACME2CYL_DEFAULT_POISSON_BOUNDARY_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>

template<typename GV, typename RF, typename GF_INITIAL_POT>
class DefaultPoissonBoundary
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    enum SideBoundary { Cytosol = 0, Debye = 1, Extracellular = 2};

    DefaultPoissonBoundary(const Acme2CylParameters& params_, const GF_INITIAL_POT& gfInitialPot_)
    : params(params_),
      gfInitialPot(gfInitialPot_),
      yMemb(params.yMemb()),
      topDirichlet(params_.isTopBoundaryDirichlet_Potential()),
      bottomDirichlet(params_.isBottomBoundaryDirichlet_Potential()),
      leftCytosolDirichlet(params_.isLeftCytosolBoundaryDirichlet_Potential()),
      rightCytosolDirichlet(params_.isRightCytosolBoundaryDirichlet_Potential()),
      leftDebyeDirichlet(params_.isLeftDebyeBoundaryDirichlet_Potential()),
      rightDebyeDirichlet(params_.isRightDebyeBoundaryDirichlet_Potential()),
      leftExtracellularDirichlet(params_.isLeftExtracellularBoundaryDirichlet_Potential()),
      rightExtracellularDirichlet(params_.isRightExtracellularBoundaryDirichlet_Potential()),
      debyeLayerWidth(params_.boundary.get("debyeLayerWidth",0.)),
      dMemb(params.dMemb()),
      useMembrane(params.useMembrane()),
      xmin(params.xMin()),
      xmax(params.xMax()),
      ymin(params.yMin()),
      ymax(params.yMax()),
      nMembranes(params.nMembranes())
    {}

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, bool isMembraneInterface) const
    {
      BCType bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;

      // Fix potential at all boundaries when there is no membrane
      if(not useMembrane)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
      }

      if(isMembraneInterface)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::None;
      }

      bool top = is.geometry().global(x)[1]+1e-6 > ymax;
      bool bottom = is.geometry().global(x)[1]-1e-6 < ymin;
      bool left = is.geometry().global(x)[0]-1e-6 < xmin;
      bool right = is.geometry().global(x)[0]+1e-6 > xmax;

      SideBoundary sideBoundary = SideBoundary::Extracellular;

      if(useMembrane && (left || right))
      {
        typename Traits::DomainFieldType yglobal = is.geometry().global(x)[1];
        switch(nMembranes)
        {
          case 1:
          {
            // Side boundaries below membrane are cytosol borders
            if (yglobal-1e-6 < yMemb[0]-debyeLayerWidth)
            {
              sideBoundary = SideBoundary::Cytosol;
              break;
            }
          }
          case 2:
          {
            // Side boundaries between membranes are cytosol borders
            if (yglobal+1e-6 > yMemb[0]+dMemb+debyeLayerWidth && yglobal-1e-6 < yMemb[1]-debyeLayerWidth)
            {
              sideBoundary = SideBoundary::Cytosol;
              break;
            }
          }
          case 3:
          {
            // Side boundaries below membrane 1 and between membranes 2 and 3 are cytosol borders
            if(yglobal-1e-6 < yMemb[0]-debyeLayerWidth ||
                  (yglobal+1e-6 > yMemb[1]+dMemb+debyeLayerWidth && yglobal-1e-6 < yMemb[2]-debyeLayerWidth))
            {
              sideBoundary = SideBoundary::Cytosol;
              break;
            }
          }
          // Check here for Debye layer boundary
          default:
          {
            if(nMembranes < 1 || nMembranes > 3)
              DUNE_THROW(Dune::NotImplemented, "Boundary types for more than 3 membranes not implemented!");

            for(int m = 0; m<nMembranes; m++)
            {
              // Point within the range of one debyeLayerWidth around the membrane are debye borders
              if (yglobal+1e-6 > yMemb[m]-debyeLayerWidth && yglobal-1e-6 < yMemb[m]+dMemb+debyeLayerWidth)
              {
                sideBoundary = SideBoundary::Debye;
                break;
              }
            }
          }
        }
      }

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
      // *** Handle LEFT boundary ***
      if(left)
      {
        if(sideBoundary == SideBoundary::Cytosol)
        {
          bctype = (leftCytosolDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                         : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        } else if (sideBoundary == SideBoundary::Debye) {
          bctype = (leftDebyeDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                       : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        } else {
          bctype = (leftExtracellularDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                               : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        }
        return bctype;
      }

      // *** Handle RIGHT boundary ***
      if(right)
      {
        if(sideBoundary == SideBoundary::Cytosol)
        {
          bctype = (rightCytosolDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                         : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        } else if (sideBoundary == SideBoundary::Debye) {
          bctype = (rightDebyeDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                       : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        } else {
          bctype = (rightExtracellularDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                               : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        }
        return bctype;
      }

      debug_jochen << "x = " << is.geometry().global(x) << " => POT bctype = " << bctype << std::endl;
      //return bctype;
      DUNE_THROW(Dune::Exception, "Error checking this boundary's position!");
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, double time) const
    {
      //typename Traits::DomainType xglobal = e.geometry().global(x);
      bool top = e.geometry().global(x)[1]+1e-6 > ymax;
      bool bottom = e.geometry().global(x)[1]-1e-6 < ymin;
      bool left = e.geometry().global(x)[0]-1e-6 < xmin;
      bool right = e.geometry().global(x)[0]+1e-6 > xmax;

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
    const Acme2CylParameters& params;
    const GF_INITIAL_POT& gfInitialPot;
    const std::vector<double> yMemb;

    const bool topDirichlet;
    const bool bottomDirichlet;
    const bool leftCytosolDirichlet;
    const bool rightCytosolDirichlet;
    const bool leftDebyeDirichlet;
    const bool rightDebyeDirichlet;
    const bool leftExtracellularDirichlet;
    const bool rightExtracellularDirichlet;

    const double debyeLayerWidth;
    const double dMemb;

    const bool useMembrane;
    const double xmin;
    const double xmax;
    const double ymin;
    const double ymax;
    const int nMembranes;
};

#endif /* DUNE_AX1_ACME2CYL_DEFAULT_POISSON_BOUNDARY_HH */
