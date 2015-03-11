/*
 * nernst_planck_parameters.hh
 *
 *  Created on: Sep 1, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME0_NERNST_PLANCK_PARAMETERS_HH
#define DUNE_AX1_ACME0_NERNST_PLANCK_PARAMETERS_HH

#include <limits>

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>

template<typename GV, typename RF, typename PHYSICS, typename DGF_POT_GRAD, typename NERNST_PLANCK_BOUNDARY>
class NernstPlanckParameters
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    NernstPlanckParameters(PHYSICS& physics_, DGF_POT_GRAD& dgfPotGrad_,
        NERNST_PLANCK_BOUNDARY& nernstPlanckBoundary_)
      : physics(physics_),
        dgfPotGrad(dgfPotGrad_),
        nernstPlanckBoundary(nernstPlanckBoundary_)
    {
      // Numer of local power function space must equal number of ion species
      assert(NUMBER_OF_SPECIES == physics.numOfSpecies());
    }

    //! tensor diffusion coefficient
    typename Traits::PermTensorType
    A (const typename Traits::ElementType& e, const typename Traits::DomainType& x, int ionSpecies_) const
    {
      typename Traits::PermTensorType I;

      int elemIndex = physics.getElementIndex(e);
      //int subdomainIndex = physics.getSubdomainIndex(elemIndex);

      for (std::size_t i=0; i<Traits::dimDomain; i++)
      {
        for (std::size_t j=0; j<Traits::dimDomain; j++)
        {
          I[i][j] = (i==j) ? physics.getDiffCoeff(ionSpecies_, elemIndex) : 0;
          //I[i][j] = 0;
          //debug_verb << "I[" << i << "][" << j << "] = " << I[i][j] << std::endl;
        }
      }

      //debug_verb << "elem[" << elemIndex << "] D = " << I << std::endl;
      return I;
    }


    //! tensor diffusion coefficient
    typename Traits::PermTensorType
    A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      return this->A(e, x, ionSpecies);
    }


    //! velocity field
    typename Traits::RangeType
    b (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        int ionSpecies_, typename Traits::RangeType potGrad_ =
            typename Traits::RangeType(std::numeric_limits<typename Traits::RangeFieldType>::max())) const
    {
      //typename Traits::DomainType xglobal = e.geometry().global(x);
      typename Traits::RangeType v(0.0);

      int elemIndex = physics.getElementIndex(e);
      //int subdomainIndex = physics.getSubdomainIndex(elemIndex);

      typename Traits::RangeType D = physics.getDiffCoeff(ionSpecies_, elemIndex);
      typename Traits::RangeType z = physics.getValence(ionSpecies_);
      typename Traits::RangeType potGrad = potGrad_;

      if(potGrad == std::numeric_limits<typename Traits::RangeFieldType>::max())
      {
#ifndef ACME0_FULLY_IMPLICIT
        dgfPotGrad.evaluate(e, x, potGrad);
#else
        debug_warn << "You must hand over the current potential gradient for the fully-implicit method to work!" << std::endl;
        DUNE_THROW(Dune::Exception, "Potential gradient was not handed over to NernstPlanckParameters::b");
#endif
      }

      v = -D * z * potGrad;
      //v = 1.0;
      //debug_verb << "elem[" << elemIndex << "] v = " << v << std::endl;
      return v;
    }


    //! velocity field
    typename Traits::RangeType
    b (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeType potGrad
          = typename Traits::RangeType(std::numeric_limits<typename Traits::RangeFieldType>::max())) const
    {
      return this->b(e, x, ionSpecies, potGrad);
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
      return 0.0;
    }

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      return nernstPlanckBoundary.bctype(is, x, time, ionSpecies);
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      //typename Traits::DomainType xglobal = e.geometry().global(x);
      //typename Traits::RangeFieldType norm = xglobal.two_norm2();
      return nernstPlanckBoundary.g(e, x, time, ionSpecies);
    }

    //! Neumann boundary condition
    typename Traits::RangeFieldType
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      return nernstPlanckBoundary.j(is, x, time, ionSpecies);
    }

    //! outflow boundary condition
    typename Traits::RangeFieldType
    o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      return nernstPlanckBoundary.o(is, x, time, ionSpecies);
    }

    //! set time
    void setTime (double t)
    {
      time = t;
    }

    void setIonSpecies(int ionSpecies_)
    {
      ionSpecies = ionSpecies_;
    }

    const PHYSICS& getPhysics() const
    {
      return physics;
    }

    std::string getName() const
    {
      return std::string("NernstPlanckParameters");
    }

  private:
    PHYSICS& physics;
    DGF_POT_GRAD& dgfPotGrad;
    NERNST_PLANCK_BOUNDARY& nernstPlanckBoundary;
    double time;
    int ionSpecies;
};

#endif /* DUNE_AX1_ACME0_NERNST_PLANCK_PARAMETERS_HH */
