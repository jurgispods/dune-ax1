/*
 * default_nernst_planck_boundary.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_DEFAULT_NERNST_PLANCK_BOUNDARY_HH
#define DUNE_AX1_ACME2CYL_DEFAULT_NERNST_PLANCK_BOUNDARY_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>

template<typename GV, typename RF, typename GF_INITIAL_CON>
class DefaultNernstPlanckBoundary
{

  public:
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    typedef GF_INITIAL_CON GFInitialCon;

    enum SideBoundary { Cytosol = 0, Debye = 1, Extracellular = 2};

    DefaultNernstPlanckBoundary(const Acme2CylParameters& params_, const GF_INITIAL_CON& gfInitialCon_)
      : params(params_),
        gfInitialCon(gfInitialCon_),
        yMemb(params.yMemb()),
        topDirichlet(params_.isTopBoundaryDirichlet_Concentration()),
        bottomDirichlet(params_.isBottomBoundaryDirichlet_Concentration()),
        leftCytosolDirichlet(params_.isLeftCytosolBoundaryDirichlet_Concentration()),
        rightCytosolDirichlet(params_.isRightCytosolBoundaryDirichlet_Concentration()),
        leftDebyeDirichlet(params_.isLeftDebyeBoundaryDirichlet_Concentration()),
        rightDebyeDirichlet(params_.isRightDebyeBoundaryDirichlet_Concentration()),
        leftExtracellularDirichlet(params_.isLeftExtracellularBoundaryDirichlet_Concentration()),
        rightExtracellularDirichlet(params_.isRightExtracellularBoundaryDirichlet_Concentration()),
        debyeLayerWidth(params_.boundary.get("debyeLayerWidth",0.)),
        dMemb(params.dMemb()),
        useMembrane(params.useMembrane()),
        xmin(params.xMin()),
        xmax(params.xMax()),
        ymin(params.yMin()),
        ymax(params.yMax()),
        nMembranes(params.nMembranes())
    {
      if(params.nMembranes() == 1 && bottomDirichlet)
        DUNE_THROW(Dune::Exception, "Lower boundary should be zero when using only one membrane!");
    }

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, int ionSpecies, bool isMembraneInterface) const
    {
      // Remember: BC types cannot be time-dependent, so the type set at the beginning will last
      // through the entire simulation run!

      // default: Dirichlet-bulk
      BCType bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;

      // Neumann-0 when there is no membrane
      if(not useMembrane)
      {
        bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
        return bctype;
      } else {

        // *** Handle MEMBRANE boundaries inside domain ***
        // Neumann-HH on the membrane interface
        if(isMembraneInterface)
        {
          bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
//          debug_jochen << "[membInterface] "<< "x = " << is.geometry().global(x)
//              << " => " << ION_NAMES[ionSpecies] << " bctype = " << bctype << std::endl;
          return bctype;
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

              for(int m = 0; m<params.nMembranes(); m++)
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

        // *** Handle TOP boundary ***
        if(top)
        {
          bctype = (topDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                 : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
//          debug_jochen << "[top] "<< "x = " << is.geometry().global(x)
//            << " => " << ION_NAMES[ionSpecies] << " bctype = " << bctype << std::endl;
          return bctype;
        }
        // *** Handle BOTTOM boundary ***
        if(bottom)
        {
          bctype = (bottomDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                    : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
//          debug_jochen << "[bottom] "<< "x = " << is.geometry().global(x)
//                        << " => " << ION_NAMES[ionSpecies] << " bctype = " << bctype << std::endl;
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
//          debug_jochen << "[left] "<< "x = " << is.geometry().global(x)
//                        << " => " << ION_NAMES[ionSpecies] << " bctype = " << bctype << std::endl;
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
//          debug_jochen << "[right] "<< "x = " << is.geometry().global(x)
//                        << " => " << ION_NAMES[ionSpecies] << " bctype = " << bctype << std::endl;
          return bctype;
        }
      }

      debug_jochen << "x = " << is.geometry().global(x) << " => CON bctype = " << bctype << std::endl;
      return bctype;

      DUNE_THROW(Dune::Exception, "Error checking this boundary's position!");
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        double time, int ionSpecies) const
    {
      typename Traits::DomainType xglobal = e.geometry().global(x);
      //typename Traits::RangeFieldType norm = xglobal.two_norm2();
      typename GF_INITIAL_CON::Traits::RangeType conc(0.0);
      gfInitialCon.evaluate(e, x, conc);

      //debug_jochen << xglobal << "gconc = " << conc[ionSpecies] << std::endl;

      return conc[ionSpecies];
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

#endif /* DUNE_AX1_ACME2CYL_DEFAULT_NERNST_PLANCK_BOUNDARY_HH */
