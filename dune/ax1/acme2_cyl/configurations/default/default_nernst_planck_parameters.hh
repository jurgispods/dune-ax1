/*
 * nernst_planck_parameters.hh
 *
 *  Created on: Sep 1, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_DEFAULT_NERNST_PLANCK_PARAMETERS_HH
#define DUNE_AX1_ACME2CYL_DEFAULT_NERNST_PLANCK_PARAMETERS_HH

#include <limits>

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_geometrywrapper.hh>

template<typename GV, typename RF, typename PHYSICS, typename INITIAL_GF, typename GF_MEMB_FLUX, typename GF_MORI_FLUX>
class DefaultNernstPlanckParameters
{
  public:
    typedef PHYSICS Physics;

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    // Hardcoded range field type for potential and concentrations
    typedef Dune::FieldVector<typename Traits::RangeFieldType,1> PotRangeType;
    typedef typename INITIAL_GF::Traits::RangeType ConRangeType;

    enum { numSpecies = NUMBER_OF_SPECIES };

    enum SideBoundary { Cytosol = 0, Debye = 1, Extracellular = 2};

    DefaultNernstPlanckParameters(const int ionSpecies_, const GV& gv_, PHYSICS& physics_,
            INITIAL_GF& initialGF_,
            GF_MEMB_FLUX& gfMembFlux_,
            GF_MORI_FLUX& gfMoriFlux_, // dummy, not used
            double tEquilibrium_)
    : ionSpecies(ionSpecies_),
      gv(gv_),
      physics(physics_),
      params(physics_.getParams()),
      initialGF(initialGF_),
      gfMembFlux(gfMembFlux_),
      topDirichlet(params.isTopBoundaryDirichlet_Concentration()),
      bottomDirichlet(params.isBottomBoundaryDirichlet_Concentration()),
      leftCytosolDirichlet(params.isLeftCytosolBoundaryDirichlet_Concentration()),
      rightCytosolDirichlet(params.isRightCytosolBoundaryDirichlet_Concentration()),
      leftDebyeDirichlet(params.isLeftDebyeBoundaryDirichlet_Concentration()),
      rightDebyeDirichlet(params.isRightDebyeBoundaryDirichlet_Concentration()),
      leftExtracellularDirichlet(params.isLeftExtracellularBoundaryDirichlet_Concentration()),
      rightExtracellularDirichlet(params.isRightExtracellularBoundaryDirichlet_Concentration()),
      yMemb(params.yMemb()),
      dMemb(params.dMemb()),
      debyeLayerWidth(params.boundary.get("debyeLayerWidth",0.)),
      useMembrane(params.useMembrane()),
      xmin(params.xMin()),
      xmax(params.xMax()),
      ymin(params.yMin()),
      ymax(params.yMax()),
      nMembranes(params.nMembranes()),
      time(0.0),
      tEquilibrium(tEquilibrium_),
      initConcEx(physics.numOfSpecies()),
      initConcIn(physics.numOfSpecies()),
      doStimulation(physics.getParams().doStimulation()),
      doMembFluxStimulation(physics.getParams().stimulation.get("memb_flux_stimulation", false)),
      membFluxStimulationIon(Cl),
      tinj_start(physics.getParams().tInj_start()),
      tinj_end(physics.getParams().tInj_end()),
      stim_x(physics.getParams().stimulation.get("position_x",
          (physics.getParams().xMax() - 3.0) / 2.0)),
      stim_y(physics.getParams().stimulation.get("position_y",
          (physics.getParams().yMemb()[0] - 3.0) / 2.0)),
      Iinj(physics.getParams().stimulation.get("I_inj", 0.0)),
      idStr(std::string("[" + ION_NAMES[ionSpecies] + "] "))
    {
      // Numer of local power function space must equal number of ion species
      assert(NUMBER_OF_SPECIES == physics.numOfSpecies());

      for(int i=0; i<initConcEx.size(); ++i)
      {
        initConcEx[i] = params.solution_ex.get(ION_NAMES[i], -1.);
        initConcIn[i] = params.solution_in.get(ION_NAMES[i], -1.);
      }
      //std::string initConcStr = params.initConc();

      std::string stimIon = physics.getParams().stimulation.get("memb_flux_ion", std::string("Cl"));
      if(stimIon == "Na")
        membFluxStimulationIon = Na;
      if(stimIon == "K")
        membFluxStimulationIon = K;
    }

    //! tensor diffusion coefficient
    typename Traits::PermTensorType
    A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      typename Traits::PermTensorType I;

      int elemIndex = physics.getElementIndex(e);
      //int subdomainIndex = physics.getSubdomainIndex(elemIndex);

      for (std::size_t i=0; i<Traits::dimDomain; i++)
      {
        for (std::size_t j=0; j<Traits::dimDomain; j++)
        {
          I[i][j] = (i==j) ? physics.getDiffCoeff(ionSpecies, elemIndex) : 0;
          //I[i][j] = 0;
          //debug_jochen << idStr << "I[" << i << "][" << j << "] = " << I[i][j] << std::endl;
        }
      }

      //debug_verb << idStr << "[NernstPlanckParameters] @ " << e.geometry().global(x) << ", A = " << I << std::endl;
      return I;
    }

    //! velocity field
    typename Traits::RangeType
    b (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeType& potGrad) const
    {
      //typename Traits::DomainType xglobal = e.geometry().global(x);
      typename Traits::RangeType v(0.0);

      int elemIndex = physics.getElementIndex(e);
      //int subdomainIndex = physics.getSubdomainIndex(elemIndex);

      typename Traits::RangeFieldType D = physics.getDiffCoeff(ionSpecies, elemIndex);

      //debug_jochen << idStr << " diff coeff = " << D << std::endl;

      int z = physics.getValence(ionSpecies);

      for(int i=0; i<v.size(); i++)
      {
        v[i] = -D * z * potGrad[i];
      }
      //v = 1.0;
      //debug_verb << idStr << "[NernstPlanckParameters] @ " << e.geometry().global(x) << ", b = " << v << std::endl;
      return v;
    }


    typename Traits::RangeType
    dbdu (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeType& phiGrad) const
    {
      int elemIndex = physics.getElementIndex(e);
      typename Traits::RangeFieldType D = physics.getDiffCoeff(ionSpecies, elemIndex);
      int z = physics.getValence(ionSpecies);

      typename Traits::RangeType v = phiGrad;
      v *= -D * z;

      return v;
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
      if(doStimulation
          && ((time > tinj_start && time < tinj_end)
          /*|| (time > 15e3+physics.getParams().tInj_start() && time < 15e3+physics.getParams().tInj_end())*/))
      {
        // TODO This is so horribly unperformant...

        // Wrap it up!
        typename Acme2CylGeometrySwitch::GeometrySwitch<typename Traits::ElementType::Geometry>::type
          geo(e.geometry());

        //typename Traits::DomainType xglobal = geo.global(x);
        //typename Traits::RangeFieldType norm = xglobal.two_norm2();

        // Get position of injection from params file; default is (near) cytosol center
        typename Traits::DomainType injectionPosition;
        injectionPosition[0] = stim_x;
        injectionPosition[1] = stim_y;

        typename Traits::DomainType firstCorner = geo.corner(0);
        typename Traits::DomainType lastCorner = geo.corner(geo.corners()-1);

        bool elementContainsInjectionPosition = Tools::lessOrEqualThan(firstCorner, injectionPosition)
            && Tools::greaterOrEqualThan(lastCorner, injectionPosition);

        //debug_jochen << idStr << "### " << physics.getElementIndex(e) << " " << firstCorner
        //             << " <= " << injectionPosition << " <= " << lastCorner << std::endl;

        // Inject sodium into this element
        if(elementContainsInjectionPosition && ionSpecies == Na)
        {
          double con_inj = Iinj;
          con_inj /= geo.volume();

          // Bring to the right units (using microsecond time scale)
          con_inj *= physics.getTimeScale();

          //debug_jochen << idStr << "### x = "<< xglobal << ", volume = " << geo.volume() << std::endl;
          //debug_jochen << idStr << "### Increasing " << ION_NAMES[ionSpecies_]
          //    << " concentration by " << con_inj << std::endl;

          return con_inj;
        }
      }

      return 0.0;
    }

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      bool isMembraneInterface = physics.isMembraneInterface(is);

      //debug_verb << idStr << "Nernst-Planck BOUNDARY @ x = " << is.geometry().global(x) << ": "
      //    << nernstPlanckBoundary.bctype(is, x, time, ionSpecies_, isMembraneInterface) << std::endl;

      // default: Dirichlet-bulk
      BCType bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;

      // Default: Neumann-0 when there is no membrane
      if(not useMembrane)
      {
        bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      } else {

        // *** Handle MEMBRANE boundaries inside domain ***
        // We currently only support Neumann bctypes on the membrane interface
        if(isMembraneInterface)
        {
          bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  //          debug_jochen << idStr << "[membInterface] "<< "x = " << is.geometry().global(x)
  //              << " => " << ION_NAMES[ionSpecies] << " bctype = " << bctype << std::endl;
          return bctype;
        }
      }

      // Check boundary location
      bool top = is.geometry().global(x)[1]+1e-6 > ymax;
      bool bottom = is.geometry().global(x)[1]-1e-6 < ymin;
      bool left = is.geometry().global(x)[0]-1e-6 < xmin;
      bool right = is.geometry().global(x)[0]+1e-6 > xmax;

      SideBoundary sideBoundary = SideBoundary::Extracellular;

      if(useMembrane && (left || right))
      {
        typename Traits::DomainFieldType yglobal = is.geometry().global(x)[1];
        switch(nMembranes)
        {
          case 1:
          {
            // Side boundaries below membrane are cytosol borders
            if (yglobal-1e-6 < yMemb[0]-debyeLayerWidth)
            {
              sideBoundary = SideBoundary::Cytosol;
              break;
            }
          }
          case 2:
          {
            // Side boundaries between membranes are cytosol borders
            if (yglobal+1e-6 > yMemb[0]+dMemb+debyeLayerWidth && yglobal-1e-6 < yMemb[1]-debyeLayerWidth)
            {
              sideBoundary = SideBoundary::Cytosol;
              break;
            }
          }
          case 3:
          {
            // Side boundaries below membrane 1 and between membranes 2 and 3 are cytosol borders
            if(yglobal-1e-6 < yMemb[0]-debyeLayerWidth ||
                  (yglobal+1e-6 > yMemb[1]+dMemb+debyeLayerWidth && yglobal-1e-6 < yMemb[2]-debyeLayerWidth))
            {
              sideBoundary = SideBoundary::Cytosol;
              break;
            }
          }
          // Check here for Debye layer boundary
          default:
          {
            if(nMembranes < 1 || nMembranes > 3)
              DUNE_THROW(Dune::NotImplemented, "Boundary types for more than 3 membranes not implemented!");

            for(int m = 0; m<params.nMembranes(); m++)
            {
              // Point within the range of one debyeLayerWidth around the membrane are debye borders
              if (yglobal+1e-6 > yMemb[m]-debyeLayerWidth && yglobal-1e-6 < yMemb[m]+dMemb+debyeLayerWidth)
              {
                sideBoundary = SideBoundary::Debye;
                break;
              }
            }
          }
        }
      }

      // *** Handle TOP boundary ***
      if(top)
      {
        bctype = (topDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                               : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
//          debug_jochen << idStr << "[top] "<< "x = " << is.geometry().global(x)
//            << " => " << ION_NAMES[ionSpecies] << " bctype = " << bctype << std::endl;
        return bctype;
      }
      // *** Handle BOTTOM boundary ***
      if(bottom)
      {
        bctype = (bottomDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                  : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
//          debug_jochen << idStr << "[bottom] "<< "x = " << is.geometry().global(x)
//                        << " => " << ION_NAMES[ionSpecies] << " bctype = " << bctype << std::endl;
        return bctype;
      }
      // *** Handle LEFT boundary ***
      if(left)
      {
        if(sideBoundary == SideBoundary::Cytosol)
        {
          bctype = (leftCytosolDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                         : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        } else if (sideBoundary == SideBoundary::Debye) {
          bctype = (leftDebyeDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                       : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        } else {
          bctype = (leftExtracellularDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                               : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        }
//          debug_jochen << idStr << "[left] "<< "x = " << is.geometry().global(x)
//                        << " => " << ION_NAMES[ionSpecies] << " bctype = " << bctype << std::endl;
        return bctype;
      }

      // *** Handle RIGHT boundary ***
      if(right)
      {
        if(sideBoundary == SideBoundary::Cytosol)
        {
          bctype = (rightCytosolDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                         : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        } else if (sideBoundary == SideBoundary::Debye) {
          bctype = (rightDebyeDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                       : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        } else {
          bctype = (rightExtracellularDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                               : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        }
//          debug_jochen << idStr << "[right] "<< "x = " << is.geometry().global(x)
//                        << " => " << ION_NAMES[ionSpecies] << " bctype = " << bctype << std::endl;
        return bctype;
      }


      debug_jochen << idStr << "x = " << is.geometry().global(x) << " => CON bctype = " << bctype << std::endl;

      DUNE_THROW(Dune::Exception, "Error checking this boundary's position!");
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      //int elemIndex = physics.getElementIndex(e);
      int subdomainIndex = physics.getSubdomainIndex(e);

      typename Traits::RangeFieldType y(0.0);

      // Default initialization: Homogeneous concentrations on each subdomain
      switch(subdomainIndex)
      {
        case CYTOSOL:
        {
          y = initConcIn[ionSpecies];
          break;
        }
        case ES:
        {
          y = initConcEx[ionSpecies];
          break;
        }
        case MEMBRANE:
        {
          if(time > 0.0)
          {
            DUNE_THROW(Dune::Exception, "Element is neither ES nor CYTOSOL!");
          }

          // No charge carriers inside the membrane!
          y = 0.0;
          // Enforce this, return immediately
          return y;
        }
        default:
          DUNE_THROW(Dune::Exception, "Element has an unknown subdomain index!");
      }

      // Override default behaviour by evaluating custom initial gridfunction
      // (default InitialVoid does nothing, so yCustom[ionSpecies_] will be equal to y
      typename INITIAL_GF::Traits::DomainType xglobal = e.geometry().global(x);
      typename INITIAL_GF::Traits::RangeType yCustom(y);
      initialGF.evaluateGlobal(xglobal,yCustom);

      return yCustom[ionSpecies];
    }

    //! Neumann boundary condition
    typename Traits::RangeFieldType
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        bool includeMoriFlux = false,
        bool implicit = false,
        const PotRangeType& uPot = PotRangeType(-1.0),
        const ConRangeType& uCon = ConRangeType(-1.0)) const
    {
      typename Traits::RangeFieldType y = 0.0;

      // Determine trans-membrane flux by evaluating membrane flux gridfunction
      if(physics.isMembraneInterface(is))
      {
        typename GF_MEMB_FLUX::Traits::RangeType membFlux;
        gfMembFlux.evaluate(is,x,membFlux,GF_MEMB_FLUX::FluxTerms::TOTAL_FLUX,ionSpecies,implicit,uPot,uCon);
        y = membFlux[ionSpecies];

        /* Now the intersection geometry can be wrapped, so that we have access to the
         * real cylinder surface
         */
        typedef typename Traits::IntersectionType::Geometry GEO_ORIG;
        typename Acme2CylGeometrySwitch::GeometrySwitch<GEO_ORIG>::type igeo(is.geometry());

        if(doMembFluxStimulation && (time > tinj_start && time < tinj_end) && ionSpecies == membFluxStimulationIon)
        {
          // Get position of injection from params file; default is (near) cytosol center
          typename Traits::DomainType injectionPosition;
          injectionPosition[0] = stim_x;
          injectionPosition[1] = stim_y;

          typename Traits::DomainType firstCorner = igeo.corner(0);
          typename Traits::DomainType lastCorner = igeo.corner(igeo.corners()-1);

          bool elementContainsInjectionPosition = (firstCorner[0]-1e-6 < stim_x && lastCorner[0]+1e-6 > stim_x);

          //debug_jochen << idStr << "### " << physics.getElementIndex(e) << " " << firstCorner
          //             << " <= " << injectionPosition << " <= " << lastCorner << std::endl;


          // Inject sodium into this element
          if(elementContainsInjectionPosition)
          {
            RF stimFlux = params.stimulation.get("memb_flux_inj", 0.0) / igeo.volume();

            // Smoothly vary stimulation by cosine function between t_inj_start and t_inj_end
            // (maximum at center of interval), taken from MoriPeskin2008 and adapted such that
            // the scaling function smoothly varies between 0 and 1 on the interval [0 1]
            // TODO note: stimulation is still discontinuous in space, maybe change this as well!
            RF t_rel = (time - tinj_start) / tinj_end;
            stimFlux *= 0.5 * (1 - std::cos(2*con_pi*t_rel));

            // Multiply by unit outer normal for correct orientation
            stimFlux *= is.centerUnitOuterNormal()[Traits::GridViewType::Grid::dimension - 1];

            debug_jochen << idStr << "  Adding value " << stimFlux << " to " << ION_NAMES[membFluxStimulationIon]
              << " flux" << " @ " << igeo.global(x) << std::endl;

            y += stimFlux;
          }
        }


        // Do surface scaling; membFlux should be scaled by factor (ext_surface / my_surface)
        // in order to always have the membrane flux multiplied by the extracellular membrane area,
        // even if we are on the intracellular side!
        RF my_surface = igeo.volume();
        RF ext_surface = 0.0;

        ext_surface = physics.getMembraneSurfaceArea(is);

        if(my_surface < is.geometry().volume())
        {
          debug_warn << idStr << "= My plain 2D surface: " << is.geometry().volume() << std::endl;
          debug_warn << idStr << "= My real surface: " << my_surface << std::endl;
          debug_warn << idStr << "= ES real surface: " << ext_surface << std::endl;

          DUNE_THROW(Dune::Exception, "Wrapped cylinder intersection volume is smaller than plain "
              << " 2D intersection volume!");
        }

        // Return trans-membrane flux for this ion species, scaled by membrane surfaces
        // => Integration over intersection in local operator will yield
        // flux = ext_surface * membFlux[ionSpecies_], the 'real' membrane current!
        RF surfaceScaleFactor = ext_surface / my_surface;
        y *= surfaceScaleFactor;

//        debug_jochen << idStr << "BC @ [" << MD_SUBDOMAINS[physics.getSubDomainNumber(*is.inside())]
//          << "] element "<< is.inside()->geometry().center() << ": "
//          << ION_NAMES[ionSpecies_] << " normal flux = " << y << std::endl;
      }

      // Default: Zero flux
      return y;
    }

    //! Neumann boundary condition j, scaled by membrane surface area
    typename Traits::RangeFieldType
    memb_flux (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {

      typedef typename Acme2CylGeometrySwitch::GeometrySwitch<typename Traits::IntersectionType::Geometry>::type GEO;
      GEO geo(is.geometry());

      return geo.volume() * j(is,x);
    }

    //! outflow boundary condition
    typename Traits::RangeFieldType
    o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      // Default: Zero outflow
      return 0.0;
    }

    //! set time
    void setTime (double t)
    {
      time = t;
      initialGF.setTime(time);
    }

    void preAssembly(bool implicit)
    {
      // Update channels for this Newton iteration
      debug_jochen << idStr << "# gfMembFlux.updateState" << std::endl;
      gfMembFlux.updateState();
      if(!implicit)
      {
        debug_jochen << idStr << "# gfMembFlux.updateFlux" << std::endl;
        gfMembFlux.updateFlux();
      }
    }

    const PHYSICS& getPhysics() const
    {
     return physics;
    }

    std::string getName() const
    {
     return std::string("NernstPlanckParameters");
    }

    const GV& getGridView() const
    {
      return gv;
    }

    // dummy method
    void prepareNextTimeStep(double nextTime)
    {
    }


private:
    const int ionSpecies;

    const GV& gv;
    PHYSICS& physics;
    const Acme2CylParameters& params;
    INITIAL_GF& initialGF;
    GF_MEMB_FLUX& gfMembFlux;

    const bool topDirichlet;
    const bool bottomDirichlet;
    const bool leftCytosolDirichlet;
    const bool rightCytosolDirichlet;
    const bool leftDebyeDirichlet;
    const bool rightDebyeDirichlet;
    const bool leftExtracellularDirichlet;
    const bool rightExtracellularDirichlet;
    const std::vector<double> yMemb;
    const double dMemb;
    const double debyeLayerWidth;
    const bool useMembrane;
    const double xmin;
    const double xmax;
    const double ymin;
    const double ymax;
    const int nMembranes;

    double time;
    const double tEquilibrium;

    std::vector<RF> initConcEx;
    std::vector<RF> initConcIn;

    const bool doStimulation;
    const bool doMembFluxStimulation;
    int membFluxStimulationIon;
    const double tinj_start;
    const double tinj_end;
    const double stim_x;
    const double stim_y;
    const double Iinj;
    const std::string idStr;
};


#endif /* DUNE_AX1_ACME2CYL_DEFAULT_NERNST_PLANCK_PARAMETERS_HH */
