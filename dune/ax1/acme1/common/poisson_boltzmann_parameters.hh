/*
 * poisson_parameters.hh
 *
 *  Created on: Sep 9, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_POISSON_BOLTZMANN_PARAMETERS_HH
#define DUNE_AX1_POISSON_BOLTZMANN_PARAMETERS_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

template<typename GV, typename RF, typename PHYSICS, typename GF_PB_RHS, typename POISSON_BOUNDARY>
class PoissonBoltzmannParameters
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    PoissonBoltzmannParameters(PHYSICS& physics_, GF_PB_RHS& gfPoissonBoltzmannRHS_,
        POISSON_BOUNDARY& poissonBoundary_)
      : physics(physics_),
        gfPoissonBoltzmannRHS(gfPoissonBoltzmannRHS_),
        poissonBoundary(poissonBoundary_),
        poissonConstant(physics.getPoissonConstant()),
        time(0.0)
    {}

    //! tensor diffusion coefficient
    typename Traits::PermTensorType
    A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      typename Traits::PermTensorType I;

      int elemIndex = physics.getElementIndex(e);

      for (std::size_t i=0; i<Traits::dimDomain; i++)
      {
        for (std::size_t j=0; j<Traits::dimDomain; j++)
        {
          I[i][j] = (i==j) ? -physics.getPermittivity(elemIndex) : 0;
          //debug_verb << "[PoissonParameters] I[" << i << "][" << j << "] = " << I[i][j] << std::endl;
        }
      }

      //debug_verb << "elem[" << elemIndex << "] D = " << I << std::endl;
      return I;
    }

    //! velocity field
    typename Traits::RangeType
    b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      return 0.0;
    }


    //! sink term
    typename Traits::RangeFieldType
    c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      return 0.0;
    }

    //! source term
    typename Traits::RangeFieldType
    f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      //typename Traits::DomainType xglobal = e.geometry().global(x);
      //typename Traits::RangeFieldType norm = xglobal.two_norm2();

      typename GF_PB_RHS::Traits::RangeType poissonBoltzmannRHS;
      gfPoissonBoltzmannRHS.evaluate(e, x, poissonBoltzmannRHS);

      /*
      if(std::abs(e.geometry().center() + 5.35156) < 1e-3)
      {
        Dune::ios_base_all_saver laden(std::cout);
        debug_jochen << "x = " << e.geometry().global(x)
            << std::scientific << std::setprecision(16) << ", cd = " << chargeDensity << std::endl;
      }
      */

      typename Traits::RangeFieldType rhs = - poissonConstant * poissonBoltzmannRHS;

      //debug_verb << "[PoissonParameters] chargeDensity = " << chargeDensity << std::endl;
      //debug_verb << "[PoissonParameters] rhs = " << rhs << std::endl;

      return rhs;
    }

    //! source term depending on solution
    typename Traits::RangeFieldType
    f (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeFieldType& pot) const
    {
      gfPoissonBoltzmannRHS.setPot(pot);
      return f(e,x);
    }


    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      //debug_verb << "BOUNDARY @ x = " << is.geometry().global(x) << std::endl;
      return poissonBoundary.bctype(is, x, time);
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      return poissonBoundary.g(e, x, time);
    }

    //! Neumann boundary condition
    typename Traits::RangeFieldType
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      return poissonBoundary.j(is, x, time);
    }

    //! outflow boundary condition
    typename Traits::RangeFieldType
    o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      return poissonBoundary.o(is, x, time);
    }

    //! set time
    void setTime (double t)
    {
      time = t;
    }

    const PHYSICS& getPhysics() const
    {
      return physics;
    }

    std::string getName() const
    {
      return std::string("PoissonParameters");
    }


  private:
    PHYSICS& physics;
    GF_PB_RHS& gfPoissonBoltzmannRHS;
    POISSON_BOUNDARY& poissonBoundary;
    RF poissonConstant;
    double time;
};

#endif /* DUNE_AX1_POISSON_BOLTZMANN_PARAMETERS_HH */
