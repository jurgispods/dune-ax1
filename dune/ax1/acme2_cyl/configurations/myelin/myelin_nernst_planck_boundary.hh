/*
 * myelin_nernst_planck_boundary.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_MYELIN_NERNST_PLANCK_BOUNDARY_HH
#define DUNE_AX1_ACME2CYL_MYELIN_NERNST_PLANCK_BOUNDARY_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>

template<typename GV, typename RF, typename GF_INITIAL_CON>
class MyelinNernstPlanckBoundary
{

  public:
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    typedef GF_INITIAL_CON GFInitialCon;


    MyelinNernstPlanckBoundary(const Acme2CylParameters& params_, const GF_INITIAL_CON& gfInitialCon_)
      : params(params_),
        gfInitialCon(gfInitialCon_),
        yMemb(params.yMemb()),
        topDirichlet(params_.isTopBoundaryDirichlet_Concentration()),
        bottomDirichlet(params_.isBottomBoundaryDirichlet_Concentration()),
        leftCytosolDirichlet(params_.isLeftCytosolBoundaryDirichlet_Concentration()),
        rightCytosolDirichlet(params_.isRightCytosolBoundaryDirichlet_Concentration()),
        leftExtracellularDirichlet(params_.isLeftExtracellularBoundaryDirichlet_Concentration()),
        rightExtracellularDirichlet(params_.isRightExtracellularBoundaryDirichlet_Concentration())
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

      // myelin: Dirichlet-bulk
      BCType bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;

      // Neumann-0 when there is no membrane
      if(not params.useMembrane())
      {
        bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
        return bctype;
      } else {

        // *** Handle MEMBRANE boundaries inside domain ***
        // Neumann-HH on the membrane interface
        if(isMembraneInterface)
        {
          bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
          return bctype;
        }

        bool top = is.geometry().global(x)[1]+1e-6 > params.yMax();
        bool bottom = is.geometry().global(x)[1]-1e-6 < params.yMin();
        bool left = is.geometry().global(x)[0]-1e-6 < params.xMin();
        bool right = is.geometry().global(x)[0]+1e-6 > params.xMax();

        bool cytosolBoundary = false;
        if(params.useMembrane())
        {
          typename Traits::DomainFieldType yglobal = is.geometry().global(x)[1];
          switch(params.nMembranes())
          {
            case 1:
            {
              // Side boundaries below membrane are cytosol borders
              cytosolBoundary = (left || right)
                && yglobal-1e-6 < yMemb[0]; //- params.dMemb();
              break;
            }
            case 2:
            {
              // Side boundaries between membranes are cytosol borders
              cytosolBoundary = (left || right)
                && yglobal+1e-6 > yMemb[0] //+ params.dMemb()
                && yglobal-1e-6 < yMemb[1]; //+ params.dMemb() ;
              break;
            }
            case 3:
            {
              // Side boundaries below membrane 1 and between membranes 2 and 3 are cytosol borders
              cytosolBoundary = (left || right)
                && (yglobal-1e-6 < yMemb[0] ||
                    (yglobal+1e-6 > yMemb[1] && yglobal-1e-6 < yMemb[2]));
              break;
            }
            default:
              DUNE_THROW(Dune::NotImplemented, "Boundary types for more than 3 membranes not implemented!");
          }
        }

        // *** Handle TOP boundary ***
        if(top)
        {
          bctype = (topDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                 : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
          return bctype;
        }
        // *** Handle BOTTOM boundary ***
        if(bottom)
        {
          bctype = (bottomDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                    : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
          return bctype;
        }
        // *** Handle LEFT boundary ***
        if(left)
        {
          if(cytosolBoundary)
          {
            bctype = (leftCytosolDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
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
          if(cytosolBoundary)
          {
            bctype = (rightCytosolDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                            : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
          } else {
            bctype = (rightExtracellularDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                                  : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
          }
          return bctype;
        }
      }

      debug_jochen << "x = " << is.geometry().global(x) << " => CON bctype = " << bctype << std::endl;
      //return bctype;

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
    const bool leftExtracellularDirichlet;
    const bool rightExtracellularDirichlet;
};

#endif /* DUNE_AX1_ACME2CYL_MYELIN_NERNST_PLANCK_BOUNDARY_HH */
