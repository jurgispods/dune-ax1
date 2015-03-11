/*
 * poisson_parameters.hh
 *
 *  Created on: Sep 9, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME1MD_POISSON_PARAMETERS_HH
#define DUNE_AX1_ACME1MD_POISSON_PARAMETERS_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

template<typename GV, typename RF, typename PHYSICS, typename POISSON_BOUNDARY>
class PoissonParameters
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    PoissonParameters(PHYSICS& physics_, POISSON_BOUNDARY& poissonBoundary_)
      : physics(physics_),
        poissonBoundary(poissonBoundary_),
        poissonConstant(physics.getPoissonConstant())
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

      //debug_verb << "[PoissonParameters] @ " << e.geometry().global(x) << ", A = " << I << std::endl;
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

      // Missing concentration parameter 'conc' is only allowed when we are on the membrane!
      if(conc.size() == 0)
      {
        assert(physics.isMembrane(e));
      }

      typename Traits::RangeFieldType chargeDensity(0.0);
      for(int j=0; j<conc.size(); j++)
      {
        chargeDensity += physics.getValence(j) * conc[j];
      }

      typename Traits::RangeFieldType rhs = -chargeDensity * poissonConstant;

      //debug_verb << "[PoissonParameters] chargeDensity = " << chargeDensity << std::endl;
      //debug_verb << "[PoissonParameters] @ " << e.geometry().global(x) << ", rhs = " << rhs << std::endl;

      return rhs;
    }

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      bool isMembraneInterface = physics.isMembraneInterface(is);

      //debug_verb << "BOUNDARY @ x = " << is.geometry().global(x) << std::endl;
      return poissonBoundary.bctype(is, x, time, isMembraneInterface);
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
    POISSON_BOUNDARY& poissonBoundary;
    RF poissonConstant;
    double time;
};

#endif /* DUNE_AX1_ACME1MD_POISSON_PARAMETERS_HH */
