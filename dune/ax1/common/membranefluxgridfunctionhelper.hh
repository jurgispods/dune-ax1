/*
 * membranepotentialgridfunction.hh
 *
 *  Created on: Feb 17, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_MEMBRANEFLUXGRIDFUNCTIONHELPER_HH
#define DUNE_AX1_MEMBRANEFLUXGRIDFUNCTIONHELPER_HH

#include <dune/ax1/common/constants.hh>

template<typename GF_MEMB_FLUX>
class MembraneFluxGridFunction_DiffusionTerm
  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::IntersectionGridFunctionTraits<typename GF_MEMB_FLUX::Traits::GridViewType,
                                       typename GF_MEMB_FLUX::Traits::RangeFieldType,
                                       GF_MEMB_FLUX::Traits::RangeType::dimension,
                                       typename GF_MEMB_FLUX::Traits::RangeType>,
     MembraneFluxGridFunction_DiffusionTerm<GF_MEMB_FLUX> >
{
  public:
    typedef typename GF_MEMB_FLUX::Traits Traits;

    MembraneFluxGridFunction_DiffusionTerm(GF_MEMB_FLUX& gfMembFlux_)
      : gfMembFlux(gfMembFlux_)
    {
    }

    template<typename I>
    inline void evaluate (const I &is,
                            const typename Traits::DomainType &x,
                            typename Traits::RangeType &y) const
    {
      return gfMembFlux.evaluateDiffusionTerm(is,x,y);
    }


    inline const typename Traits::GridViewType& getGridView () const
    {
      return gfMembFlux.getGridView();
    }


  private:
    GF_MEMB_FLUX& gfMembFlux;


};

template<typename GF_MEMB_FLUX>
class MembraneFluxGridFunction_DriftTerm
  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::IntersectionGridFunctionTraits<typename GF_MEMB_FLUX::Traits::GridViewType,
                                       typename GF_MEMB_FLUX::Traits::RangeFieldType,
                                       GF_MEMB_FLUX::Traits::RangeType::dimension,
                                       typename GF_MEMB_FLUX::Traits::RangeType>,
     MembraneFluxGridFunction_DiffusionTerm<GF_MEMB_FLUX> >
{
  public:
    typedef typename GF_MEMB_FLUX::Traits Traits;

    MembraneFluxGridFunction_DriftTerm(GF_MEMB_FLUX& gfMembFlux_)
      : gfMembFlux(gfMembFlux_)
    {}

    template<typename I>
    inline void evaluate (const I &is,
                            const typename Traits::DomainType &x,
                            typename Traits::RangeType &y) const
    {
      return gfMembFlux.evaluateDriftTerm(is,x,y);
    }


    inline const typename Traits::GridViewType& getGridView () const
    {
      return gfMembFlux.getGridView();
    }


  private:
    GF_MEMB_FLUX& gfMembFlux;


};

template<typename GF_MEMB_FLUX>
class MembraneFluxGridFunction_LeakFlux
  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::IntersectionGridFunctionTraits<typename GF_MEMB_FLUX::Traits::GridViewType,
                                       typename GF_MEMB_FLUX::Traits::RangeFieldType,
                                       GF_MEMB_FLUX::Traits::RangeType::dimension,
                                       typename GF_MEMB_FLUX::Traits::RangeType>,
                                       MembraneFluxGridFunction_LeakFlux<GF_MEMB_FLUX> >
{
  public:
    typedef typename GF_MEMB_FLUX::Traits Traits;

    MembraneFluxGridFunction_LeakFlux(GF_MEMB_FLUX& gfMembFlux_)
      : gfMembFlux(gfMembFlux_)
    {}

    template<typename I>
    inline void evaluate (const I &is,
                            const typename Traits::DomainType &x,
                            typename Traits::RangeType &y) const
    {
      return gfMembFlux.evaluateLeakFlux(is,x,y);
    }


    inline const typename Traits::GridViewType& getGridView () const
    {
      return gfMembFlux.getGridView();
    }


  private:
    GF_MEMB_FLUX& gfMembFlux;

};


template<typename GF_MEMB_FLUX>
class MembraneFluxGridFunction_VoltageGatedFlux
  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::IntersectionGridFunctionTraits<typename GF_MEMB_FLUX::Traits::GridViewType,
                                       typename GF_MEMB_FLUX::Traits::RangeFieldType,
                                       GF_MEMB_FLUX::Traits::RangeType::dimension,
                                       typename GF_MEMB_FLUX::Traits::RangeType>,
                                       MembraneFluxGridFunction_VoltageGatedFlux<GF_MEMB_FLUX> >
{
  public:
    typedef typename GF_MEMB_FLUX::Traits Traits;

    MembraneFluxGridFunction_VoltageGatedFlux(GF_MEMB_FLUX& gfMembFlux_)
      : gfMembFlux(gfMembFlux_)
    {}

    template<typename I>
    inline void evaluate (const I &is,
                            const typename Traits::DomainType &x,
                            typename Traits::RangeType &y) const
    {
      return gfMembFlux.evaluateVoltageGatedFlux(is,x,y);
    }


    inline const typename Traits::GridViewType& getGridView () const
    {
      return gfMembFlux.getGridView();
    }


  private:
    GF_MEMB_FLUX& gfMembFlux;

};


template<typename GF_MORI_FLUX, typename DGF_POT, typename PHYSICS>
class MoriChargeGridFunction
  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::IntersectionGridFunctionTraits<typename GF_MORI_FLUX::Traits::GridViewType,
                                       typename GF_MORI_FLUX::Traits::RangeFieldType,
                                       GF_MORI_FLUX::Traits::RangeType::dimension,
                                       typename GF_MORI_FLUX::Traits::RangeType>,
                                       MoriChargeGridFunction<GF_MORI_FLUX,DGF_POT,PHYSICS> >
{
  public:
    typedef typename GF_MORI_FLUX::Traits Traits;

    MoriChargeGridFunction(GF_MORI_FLUX& gfMoriFlux_, const DGF_POT& dgfPot_, const PHYSICS& physics_)
      : gfMoriFlux(gfMoriFlux_),
        dgfPot(dgfPot_),
        physics(physics_)
    {}

    template<typename I>
    inline void evaluate (const I &is,
                            const typename Traits::DomainType &x,
                            typename Traits::RangeType &y) const
    {
      y = 0.0;

      if(physics.isMembraneInterface(is))
      {
        //debug_jochen << " --- x = " << is.geometry().global(x) << std::endl;

//        typename PHYSICS::FieldType perm = 0.0;
//        if(physics.isMembrane(*is.inside()))
//        {
//          perm = physics.getPermittivity(*is.inside());
//        } else {
//          perm = physics.getPermittivity(*is.outside());
//        }
//        int iIndex = physics.getIntersectionIndex(is);
//
//        debug_jochen << "         iIndex = " << iIndex << ", perm = " << perm << std::endl;

        const int subDomain = physics.getSubdomainIndex(*is.inside());

        typename PHYSICS::FieldType perm = physics.getMembranePermittivity(is);
        typename DGF_POT::Traits::RangeType membPot;
        physics.getMembranePotential(is, dgfPot, membPot);

        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          if(subDomain == ES)
            y[j] = gfMoriFlux.getNewSurfaceChargeDensity_Ext(is, j, perm, membPot);
          else if(subDomain == CYTOSOL)
            y[j] = gfMoriFlux.getNewSurfaceChargeDensity_Int(is, j, perm, membPot);
          else
            DUNE_THROW(Dune::Exception, "Expected inside element to be either ES or CYTOSOL!");
        }

        //debug_jochen << "         y = " << y << std::endl;
      }
    }


    inline const typename Traits::GridViewType& getGridView () const
    {
      return gfMoriFlux.getGridView();
    }


  private:
    GF_MORI_FLUX& gfMoriFlux;
    const DGF_POT& dgfPot;
    const PHYSICS& physics;

};

#endif /* DUNE_AX1_MEMBRANEFLUXGRIDFUNCTIONHELPER_HH */
