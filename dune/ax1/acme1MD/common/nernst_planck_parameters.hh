/*
 * nernst_planck_parameters.hh
 *
 *  Created on: Sep 1, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME1MD_NERNST_PLANCK_PARAMETERS_HH
#define DUNE_AX1_ACME1MD_NERNST_PLANCK_PARAMETERS_HH

#include <limits>

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>

template<typename GV, typename RF, typename PHYSICS, typename NERNST_PLANCK_BOUNDARY, typename GF_MEMB_FLUX>
class NernstPlanckParameters
{
  public:

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    typedef typename NERNST_PLANCK_BOUNDARY::GFInitialCon::Traits::RangeType ConRangeType;

    NernstPlanckParameters(const GV& gv_, PHYSICS& physics_,
            NERNST_PLANCK_BOUNDARY& nernstPlanckBoundary_,
            GF_MEMB_FLUX& gfMembFlux_,
            double tEquilibrium_)
          : gv(gv_),
            physics(physics_),
            nernstPlanckBoundary(nernstPlanckBoundary_),
            gfMembFlux(gfMembFlux_),
            time(0.0),
            tEquilibrium(tEquilibrium_)
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
          //debug_jochen << "I[" << i << "][" << j << "] = " << I[i][j] << std::endl;
        }
      }

      //debug_verb << "[NernstPlanckParameters] @ " << e.geometry().global(x) << ", A = " << I << std::endl;
      return I;
    }

    //! velocity field
    typename Traits::RangeType
    b (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeType potGrad, int ionSpecies_) const
    {
      //typename Traits::DomainType xglobal = e.geometry().global(x);
      typename Traits::RangeType v(0.0);

      int elemIndex = physics.getElementIndex(e);
      //int subdomainIndex = physics.getSubdomainIndex(elemIndex);

      typename Traits::RangeType D = physics.getDiffCoeff(ionSpecies_, elemIndex);

      //debug_jochen << ION_NAMES[ionSpecies_] << " diff coeff = " << D << std::endl;

      typename Traits::RangeType z = physics.getValence(ionSpecies_);
      v = -D * z * potGrad;
      //v = 1.0;
      //debug_verb << "[NernstPlanckParameters] @ " << e.geometry().global(x) << ", b = " << v << std::endl;
      return v;
    }



    //! sink term
    typename Traits::RangeFieldType
    c (const typename Traits::ElementType& e, const typename Traits::DomainType& x, int ionSpecies_) const
    {
      return 0.0;
    }

    //! source term
    typename Traits::RangeFieldType
    f (const typename Traits::ElementType& e, const typename Traits::DomainType& x, int ionSpecies_) const
    {
      typename Traits::DomainType xglobal = e.geometry().global(x);
      //typename Traits::RangeFieldType norm = xglobal.two_norm2();

      typename Traits::DomainFieldType cytosolCenter = - (physics.getParams().xMax()
          - 0.5 * physics.getParams().dMemb()) / 2.;
      //debug_jochen << "### " << physics.getElementIndex(e) << " " << cytosolCenterIndex << std::endl;

      // Hack hack!! (1D) - Inject Na into the middle of the cell
      if(physics.getParams().doStimulation()
          && (e.geometry().corner(0) <= cytosolCenter && e.geometry().corner(1) >= cytosolCenter)
          && ((time > physics.getParams().tInj_start() && time < physics.getParams().tInj_end())
              /*|| (time > 15e3+physics.getParams().tInj_start() && time < 15e3+physics.getParams().tInj_end())*/)
          && ionSpecies_ == Na)
      {
        double con_inj = physics.getParams().stimulation.get("I_inj", 0.0);
        con_inj /= e.geometry().volume();

        //debug_jochen << "### x = "<< xglobal << ", increasing " << ION_NAMES[ionSpecies_]
        //    << " concentration by " << con_inj << std::endl;

        // Scale source term to match the units used by the operator
        return (physics.getTimeScale() / physics.getLengthScale()) * con_inj;
      }

      return 0.0;
    }

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, const int ionSpecies_) const
    {
      bool isMembraneInterface = physics.isMembraneInterface(is);
      return nernstPlanckBoundary.bctype(is, x, time, ionSpecies_, isMembraneInterface);
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, int ionSpecies_) const
    {
      //typename Traits::DomainType xglobal = e.geometry().global(x);
      //typename Traits::RangeFieldType norm = xglobal.two_norm2();
      return nernstPlanckBoundary.g(e, x, time, ionSpecies_);
    }

    //! Neumann boundary condition
    typename Traits::RangeFieldType
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, int ionSpecies_) const
    {
      // Constant test flux for potassium ions
      if(physics.isMembraneInterface(is))
      {
        Dune::PDELab::IntersectionGeometry<typename Traits::IntersectionType> ig(is, -1);
        typename GF_MEMB_FLUX::Traits::RangeType membFlux;
        gfMembFlux.evaluate(ig,x,membFlux);

        //debug_jochen << "BC @ [" << MD_SUBDOMAINS[physics.getSubDomainNumber(*is.inside())] <<
        //    "] element "<< is.inside()->geometry().center() << ": "
        //    << ION_NAMES[ionSpecies_] << " normal flux = " << membFlux[ionSpecies_] << std::endl;

        // Return trans-membrane flux for this ion species
        return membFlux[ionSpecies_];

        //if(ionSpecies_ == K) return physics.getParams().equilibration.get("dummyFlux", 0.0);
        //return 0.0;
      }


        /*
          // Add Na injection on intracellular domain boundary
          if(false && is.geometry().center() < 0.0 && ionSpecies_ == Na)
          {
            // Inject 10 pA for 0.5 ms
            typename Traits::RangeFieldType Na_inj = Tools::calculateNaInjectionFlux(physics, 10., 0.5 * 1e6);
            debug_jochen << "-- Injecting flux " << Na_inj << std::endl;
            return Na_inj;
          }
        }*/


      return nernstPlanckBoundary.j(is, x, time, ionSpecies_);
    }

    //! outflow boundary condition
    typename Traits::RangeFieldType
    o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, int ionSpecies_) const
    {
      return nernstPlanckBoundary.o(is, x, time, ionSpecies_);
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
     return std::string("NernstPlanckParameters");
    }

  private:
      const GV& gv;
      PHYSICS& physics;
      NERNST_PLANCK_BOUNDARY& nernstPlanckBoundary;
      GF_MEMB_FLUX& gfMembFlux;
      double time;
      const double tEquilibrium;
};

#endif /* DUNE_AX1_ACME1MD_NERNST_PLANCK_PARAMETERS_HH */
