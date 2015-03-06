/*
 * default_nernst_planck_boundary.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2_DEFAULT_NERNST_PLANCK_BOUNDARY_HH
#define DUNE_AX1_ACME2_DEFAULT_NERNST_PLANCK_BOUNDARY_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>

template<typename GV, typename RF, typename GF_INITIAL_CON>
class DefaultNernstPlanckBoundary
{

  public:
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    typedef GF_INITIAL_CON GFInitialCon;


    DefaultNernstPlanckBoundary(const Acme2Parameters& params_, const GF_INITIAL_CON& gfInitialCon_)
      : params(params_),
        gfInitialCon(gfInitialCon_)
    {}

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, int ionSpecies, bool isMembraneInterface) const
    {
      // default: Dirichlet-bulk
      BCType bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;

      // Neumann-0 when there is no membrane
      if(not params.useMembrane())
      {
        bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      } else {

        bool top = is.geometry().global(x)[1]+1e-6 > params.yMax();
        bool bottom = is.geometry().global(x)[1]-1e-6 < params.yMin();
        bool left = is.geometry().global(x)[0]-1e-6 < params.xMin();
        bool right = is.geometry().global(x)[0]+1e-6 > params.xMax();

        bool cytosolBoundary = false;
        if(params.nMembranes() == 1 && params.useMembrane())
        {
          // Side boundaries below membrane are cytosol borders
          cytosolBoundary = (left || right)
            && is.geometry().global(x)[1]-1e-6 < params.yMemb() - params.dMemb();
        }
        if(params.nMembranes() == 2 && params.useMembrane())
        {
          // Side boundaries between membranes are cytosol borders
          cytosolBoundary = (left || right)
            && is.geometry().global(x)[1]+1e-6 > params.yMemb() + params.dMemb()
            && is.geometry().global(x)[1]-1e-6 < params.yMemb() + params.dMemb() + params.getCellWidth();
        }

        // *** Handle LEFT/RIGHT boundaries ***
        // Neumann-0 on the left and right boundaries
        if(left || right)
        {
          bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;

          // New try: If config file flag 'closedCell' is true at the start of the simulation,
          // we completely close the cell border.
          // Left and right cytosol boundaries can be opened when 'closedCell' is false.
          // Remember: BC types cannot be time-dependent, so the type set at the beginning will last
          // through the entire simulation run!
          if(cytosolBoundary && not params.closedCell())
          {
            bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
          }
        }

        // *** Handle TOP boundary (default: Dirichlet-bulk) ***
        if(top && params.nMembranes() == 1)
        {
          std::string bla = params.general.get("1memb_upperBoundary","Dirichlet");
          if(bla == "Neumann")
          {
            bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
          }
        }

        // *** Handle BOTTOM boundary ***
        //Neumann-0 on lower boundary (only for the 1-membrane case as the lower boundary
        //is then representing the inner-cell symmetry axis!
        if(bottom && params.nMembranes() == 1)
        {
          bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;

          std::string bla = params.general.get("1memb_lowerBoundary","Neumann");
          if(bla == "Dirichlet")
          {
            bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
          }
        }

        // *** Handle MEMBRANE boundaries inside domain ***
        // Neumann-HH on the membrane interface
        if(isMembraneInterface)
        {
          bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
        }
      }

      //debug_jochen << "x = " << is.geometry().global(x) << " => bctype = " << bctype << std::endl;
      return bctype;
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
    const Acme2Parameters& params;
    const GF_INITIAL_CON& gfInitialCon;
};

#endif /* DUNE_AX1_ACME2_DEFAULT_NERNST_PLANCK_BOUNDARY_HH */
