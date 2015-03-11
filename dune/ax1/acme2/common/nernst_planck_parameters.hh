/*
 * nernst_planck_parameters.hh
 *
 *  Created on: Sep 1, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2_NERNST_PLANCK_PARAMETERS_HH
#define DUNE_AX1_ACME2_NERNST_PLANCK_PARAMETERS_HH

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
        typename Traits::RangeType& potGrad, int ionSpecies_) const
    {
      //typename Traits::DomainType xglobal = e.geometry().global(x);
      typename Traits::RangeType v(0.0);

      int elemIndex = physics.getElementIndex(e);
      //int subdomainIndex = physics.getSubdomainIndex(elemIndex);

      typename Traits::RangeFieldType D = physics.getDiffCoeff(ionSpecies_, elemIndex);

      //debug_jochen << ION_NAMES[ionSpecies_] << " diff coeff = " << D << std::endl;

      int z = physics.getValence(ionSpecies_);

      for(int i=0; i<v.size(); i++)
      {
        v[i] = -D * z * potGrad[i];
      }
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
      if(physics.getParams().doStimulation()
          && ((time > physics.getParams().tInj_start() && time < physics.getParams().tInj_end())
          /*|| (time > 15e3+physics.getParams().tInj_start() && time < 15e3+physics.getParams().tInj_end())*/))
      {
        typename Traits::DomainType xglobal = e.geometry().global(x);
        //typename Traits::RangeFieldType norm = xglobal.two_norm2();

        // Get position of injection form params file; default is (near) cytosol center
        typename Traits::DomainType injectionPosition;
        injectionPosition[0] = physics.getParams().stimulation.get("position_x",
            (physics.getParams().xMax() - 3.0) / 2.0);
        injectionPosition[1] = physics.getParams().stimulation.get("position_y",
            (physics.getParams().yMemb() - 3.0) / 2.0);

        typename Traits::DomainType firstCorner = e.geometry().corner(0);
        typename Traits::DomainType lastCorner = e.geometry().corner(e.geometry().corners()-1);

        bool elementContainsInjectionPosition = Tools::lessOrEqualThan(firstCorner, injectionPosition)
            && Tools::greaterOrEqualThan(lastCorner, injectionPosition);

        //debug_jochen << "### " << physics.getElementIndex(e) << " " << firstCorner
        //             << " <= " << injectionPosition << " <= " << lastCorner << std::endl;

        // Inject sodium into this element
        if(elementContainsInjectionPosition && ionSpecies_ == Na)
        {
          double con_inj = physics.getParams().stimulation.get("I_inj", 0.0);
          con_inj /= e.geometry().volume();

          // Bring to the right units (using microsecond time scale)
          con_inj *= physics.getTimeScale();

          //debug_jochen << "### x = "<< xglobal << ", increasing " << ION_NAMES[ionSpecies_]
          //    << " concentration by " << con_inj << std::endl;


          return con_inj;
        }
      }

      return 0.0;
    }

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, const int ionSpecies_) const
    {
      bool isMembraneInterface = physics.isMembraneInterface(is);

      //debug_verb << "Nernst-Planck BOUNDARY @ x = " << is.geometry().global(x) << ": "
      //    << nernstPlanckBoundary.bctype(is, x, time, ionSpecies_, isMembraneInterface) << std::endl;

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
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        int ionSpecies_) const
    {
      // Determine trans-membrane flux by evaluating membrane flux gridfunction
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
      }

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

#endif /* DUNE_AX1_ACME2_NERNST_PLANCK_PARAMETERS_HH */
