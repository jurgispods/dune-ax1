/*
 * default_boundary.hh
 *
 *  Created on: Jan 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_BIGMAC_POISSON_BOUNDARY_HH
#define DUNE_AX1_BIGMAC_POISSON_BOUNDARY_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

template<typename GV, typename RF>
class BigmacPoissonBoundary
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

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
      typename Traits::DomainType xglobal = e.geometry().global(x);
      //typename Traits::RangeFieldType norm = xglobal.two_norm2();

      //typename INITIAL_POT::Traits::RangeType pot = 0.0;
      //initialPot.evaluate(e,x,pot);
      //debug_verb << "Gebäck" << std::endl;

    	if (xglobal > 0.0)
    	{
    		return 0.0;
    	}
    	else
    	{
    		double A = 1.0;
       	double v = 1.0;

       	return   2.0 * log ( cosh( A * ( - 5.0 - v * time ) ) ) + 5.0 * v
    		       - 2.0 * log ( cosh( A * (   5.0 - v * time ) ) ) + 5.0 * v;
    	}
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
};

#endif /* DUNE_AX1_BIGMAC_POISSON_BOUNDARY_HH */
