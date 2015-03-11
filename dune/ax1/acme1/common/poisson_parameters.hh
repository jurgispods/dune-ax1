/*
 * poisson_parameters.hh
 *
 *  Created on: Sep 9, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_POISSON_PARAMETERS_HH
#define DUNE_AX1_POISSON_PARAMETERS_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

template<typename GV, typename RF, typename PHYSICS, typename GF_CD, typename POISSON_BOUNDARY>
class PoissonParameters
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    PoissonParameters(PHYSICS& physics_, GF_CD& gfChargeDensity_, POISSON_BOUNDARY& poissonBoundary_)
      : physics(physics_),
        gfChargeDensity(gfChargeDensity_),
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
    f (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        std::vector<typename Traits::RangeFieldType> conc = std::vector<typename Traits::RangeFieldType>()) const
    {
      //typename Traits::DomainType xglobal = e.geometry().global(x);
      //typename Traits::RangeFieldType norm = xglobal.two_norm2();

      typename GF_CD::Traits::RangeType chargeDensity(0.0);

      // Use the discrete grid function to evaluate the potential gradient.
      // This is totally valid for the operator-split case!
      if(conc.size() == 0)
      {
        gfChargeDensity.evaluate(e, x, chargeDensity);
      } else {
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          chargeDensity += physics.getValence(j) * conc[j];
        }
      }
      //debug_jochen << "cd = " << chargeDensity << std::endl;

      typename Traits::RangeFieldType rhs = -chargeDensity * poissonConstant;

      //debug_verb << "[PoissonParameters] chargeDensity = " << chargeDensity << std::endl;
      //debug_verb << "[PoissonParameters] rhs = " << rhs << std::endl;

      return rhs;
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
    GF_CD& gfChargeDensity;
    POISSON_BOUNDARY& poissonBoundary;
    RF poissonConstant;
    double time;
};

#endif /* DUNE_AX1_POISSON_PARAMETERS_HH */
