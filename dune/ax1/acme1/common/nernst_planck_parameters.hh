/*
 * nernst_planck_parameters.hh
 *
 *  Created on: Sep 1, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_NERNST_PLANCK_PARAMETERS_HH
#define DUNE_AX1_NERNST_PLANCK_PARAMETERS_HH

#include <limits>

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>

template<typename GV, typename RF, typename PHYSICS,
         typename DGF_POT_GRAD, typename NERNST_PLANCK_BOUNDARY>
class NernstPlanckParameters
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    typedef typename NERNST_PLANCK_BOUNDARY::GFInitialCon::Traits::RangeType ConRangeType;

    NernstPlanckParameters(const GV& gv_, PHYSICS& physics_,
        DGF_POT_GRAD& dgfPotGrad_,
        NERNST_PLANCK_BOUNDARY& nernstPlanckBoundary_,
        double tEquilibrium_)
      : gv(gv_),
        physics(physics_),
        dgfPotGrad(dgfPotGrad_),
        nernstPlanckBoundary(nernstPlanckBoundary_),
        time(0.0),
        tEquilibrium(tEquilibrium_),
        ionSpecies(0),
        channels(physics.getMembrane().getChannelSet()),
        potJump(channels.getMembraneElements().size()),
        conJump(channels.getMembraneElements().size()),
        conUp(channels.getMembraneElements().size())
    {
      // Numer of local power function space must equal number of ion species
      assert(NUMBER_OF_SPECIES == physics.numOfSpecies());
    }

    //! tensor diffusion coefficient
    typename Traits::PermTensorType
    A (const typename Traits::ElementType& e, const typename Traits::DomainType& x, int ionSpecies_) const
    {
      typename Traits::PermTensorType I;

      int elemIndex = physics.getElementIndex(e,gv);
      //int subdomainIndex = physics.getSubdomainIndex(elemIndex);

      for (std::size_t i=0; i<Traits::dimDomain; i++)
      {
        for (std::size_t j=0; j<Traits::dimDomain; j++)
        {
          I[i][j] = (i==j) ? physics.getDiffCoeff(ionSpecies_, elemIndex) : 0;
          //I[i][j] = 0;
          //debug_jochen << "I[" << i << "][" << j << "] = " << I[i][j] << std::endl;
        }
      }

      //debug_jochen << "elem[" << elemIndex << "] D = " << I << std::endl;
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

      int elemIndex = physics.getElementIndex(e,gv);
      //int subdomainIndex = physics.getSubdomainIndex(elemIndex);

      typename Traits::RangeType D = physics.getDiffCoeff(ionSpecies_, elemIndex);
      typename Traits::RangeType z = physics.getValence(ionSpecies_);
      typename Traits::RangeType potGrad = potGrad_;

      // Use the discrete grid function to evaluate the potential gradient.
      // This is totally valid for the operator-split case!
      if(potGrad == std::numeric_limits<typename Traits::RangeFieldType>::max())
      {
        dgfPotGrad.evaluate(e, x, potGrad);
      }

      v = -D * z * potGrad;
      //v = 1.0;
      //debug_jochen << "elem[" << elemIndex << "] v = " << v << std::endl;
      return v;
    }

    //! velocity field
    typename Traits::RangeType
    b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      return this->b(e, x, ionSpecies);
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
      bool isMembraneInterface = physics.isMembraneInterface(is, gv);
      return nernstPlanckBoundary.bctype(is, x, time, ionSpecies, isMembraneInterface);
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
      // Constant test flux for potassium ions
      if(ionSpecies == K)
      {
        if(physics.isMembraneInterface(is, gv))
        {
          typename Traits::RangeFieldType totalFlux(0.0);

          int mGlobalIndex = physics.getMembraneElementIndex(is, gv);
          int mLocalIndex = physics.getLocalMembraneElementIndex(mGlobalIndex);

          debug_jochen << "Membrane interface @ " << is.geometry().center() << std::endl;
          debug_jochen << "Membrane element index " << mGlobalIndex << "[" << mLocalIndex << "]" << std::endl;

          // Fixed membrane potential, updated each time step
          typename Traits::RangeFieldType currentPotJump = potJump[mLocalIndex];
          ConRangeType currentConJump = conJump[mLocalIndex];
          ConRangeType currentConUp = conUp[mLocalIndex];
          //typename DGF_POT::Traits::RangeType potJump = physics.getMembranePotentialJump(is, gv, dgfPot);

          // Compute conductance and resulting flux for every channel type
          for(int k=0; k<channels.size(); k++)
          {
            debug_jochen << " channel #" << k;
            if(channels.getChannel(k).getIonSpecies() == ionSpecies)
            {
              typename Traits::RangeFieldType conductance = channels.getEffConductance(k, mGlobalIndex);


              // TODO Calculate reversal potential from Nernst equation in each time step (Acme1Setup)
              //typename Traits::RangeFieldType relPot = currentPotJump;
              //relPot -= channels.getChannel(k).getReversalPotential();

              //debug_jochen << ", conductance = " << conductance;
              debug_jochen << ", potJump = " << currentPotJump;

              typename Traits::RangeFieldType flux(0.0);

              // TODO Fixed conductance for testing purposes
              conductance = 36.;

              Tools::calcTransMembraneFlux(physics, conductance, currentConJump[ionSpecies],
                  currentConUp[ionSpecies], currentPotJump, ionSpecies, flux);

              debug_jochen << " => [" << ION_NAMES[ionSpecies] << "] flux = " << flux;

              flux *= is.centerUnitOuterNormal();
              debug_jochen << ", normalFlux = " << flux;

              totalFlux += flux;
            }
            debug_jochen << std::endl;
          }

          //typename Traits::RangeFieldType flux = 0.1;

          //debug_jochen << "Opposite element: " << eOpposite.geometry().center() << std::endl;
          //debug_jochen << "[" << ION_NAMES[ionSpecies] << "] flux: " << flux << std::endl;
          debug_jochen << std::endl;

          return totalFlux;
        }
      }
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

    void updatePotJump(const std::valarray<RF>& potJump_)
    {
      potJump = potJump_;
    }

    void updateConJump(const std::valarray<ConRangeType>& conJump_)
    {
      conJump = conJump_;
    }

    void updateConUp(const std::valarray<ConRangeType>& conUp_)
    {
      conUp = conUp_;
    }

  private:
    const GV& gv;
    PHYSICS& physics;
    DGF_POT_GRAD& dgfPotGrad;
    NERNST_PLANCK_BOUNDARY& nernstPlanckBoundary;
    double time;
    const double tEquilibrium;
    int ionSpecies;

    const typename PHYSICS::ChannelSet& channels;
    std::valarray<RF> potJump;
    std::valarray<ConRangeType> conJump;
    std::valarray<ConRangeType> conUp;
};

#endif /* DUNE_AX1_NERNST_PLANCK_PARAMETERS_HH */
