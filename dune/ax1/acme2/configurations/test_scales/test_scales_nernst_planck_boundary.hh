/*
 * default_nernst_planck_boundary.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2_TEST_SCALES_NERNST_PLANCK_BOUNDARY_HH
#define DUNE_AX1_ACME2_TEST_SCALES_NERNST_PLANCK_BOUNDARY_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>

template<typename GV, typename RF, typename GF_INITIAL_CON>
class TestScalesNernstPlanckBoundary
{

  public:
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    typedef GF_INITIAL_CON GFInitialCon;


    TestScalesNernstPlanckBoundary(const Acme2Parameters& params_, const GF_INITIAL_CON& gfInitialCon_)
      : params(params_),
        gfInitialCon(gfInitialCon_)
    {}

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        double time, int ionSpecies, bool isMembraneInterface) const
    {
      bool top = is.geometry().global(x)[1]+1e-6 > params.yMax();
      bool bottom = is.geometry().global(x)[1]-1e-6 < params.yMin();
      bool left = is.geometry().global(x)[0]-1e-6 < params.xMin();
      bool right = is.geometry().global(x)[0]+1e-6 > params.xMax();

      //debug_jochen << "[NP boundary] @ " << is.geometry().global(x) << ": ";

      // Neumann-0 on the left and right boundaries
      if(left || right)
      {
        //debug_jochen << "Neumann" << std::endl;
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      }
      // Neumann-0 on bottom boundary (cytosol)
      if( bottom && params.closedCell())
      {
        //debug_jochen << "Neumann" << std::endl;
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      }
      // Neumann-HH on the membrane interface
      if(isMembraneInterface)
      {
        //debug_jochen << "Neumann" << std::endl;
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      // Dirichlet-bulk on the remaining boundaries (top and bottom)
      }

      //return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      // default: Dirichlet-bulk
      //debug_jochen << "Dirichlet" << std::endl;
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
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

#endif /* DUNE_AX1_ACME2_TEST_SCALES_NERNST_PLANCK_BOUNDARY_HH */
