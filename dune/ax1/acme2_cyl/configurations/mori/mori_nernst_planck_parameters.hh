/*
 * nernst_planck_parameters.hh
 *
 *  Created on: Sep 1, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_MORI_NERNST_PLANCK_PARAMETERS_HH
#define DUNE_AX1_ACME2CYL_MORI_NERNST_PLANCK_PARAMETERS_HH

#include <limits>

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/constants.hh>

template<typename GV, typename RF, typename PHYSICS, typename INITIAL_GF, typename GF_MEMB_FLUX,
typename GF_MORI_FLUX>
class MoriNernstPlanckParameters
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

    MoriNernstPlanckParameters(const int ionSpecies_, const GV& gv_, PHYSICS& physics_,
            INITIAL_GF& initialGF_,
            GF_MEMB_FLUX& gfMembFlux_,
            GF_MORI_FLUX& gfMoriFlux_,
            double tEquilibrium_)
    : ionSpecies(ionSpecies_),
      gv(gv_),
      physics(physics_),
      params(physics_.getParams()),
      initialGF(initialGF_),
      gfMembFlux(gfMembFlux_),
      gfMoriFlux(gfMoriFlux_),
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
      useLoadedData(params.boundary.get("useTimeDependentBoundaryValuesCon",false)),
      loadBoundaryLocation(params.boundary.get("loadBoundary",std::string("bottom"))),
      isLoadBoundaryDirichlet(params.isBoundaryDirichlet_Concentration(loadBoundaryLocation)),
      collapseToLineSource(params.boundary.get("collapseToLineSource",false)),
      filename(params.boundary.get("boundaryConLoadFilename",
          std::vector<std::string>(NUMBER_OF_SPECIES,"dummy_boundary_con_bottom.dat"))[ionSpecies_]),
      con_in(filename),
      nLinesPerTimeStep(0),
      boundaryValues_old(),
      boundaryValues_new(),
      timeOld(0.0),
      timeNew(0.0),
      startLine(0),
      isElementMapSetUp(false),
      epsTime(params.boundary.get("boundaryTimePrecision",1e-4)),
      idStr(std::string("[" + ION_NAMES[ionSpecies] + "] "))
    {
      // Number of local power function space must equal number of ion species
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

      init();
    }


    void init()
    {
      if(useLoadedData)
      {
        if(loadBoundaryLocation != "top" && loadBoundaryLocation != "bottom" && loadBoundaryLocation != "membrane")
          DUNE_THROW(Dune::NotImplemented, "Loading left and right boundaries currently not implemented!");

        if(! con_in.good())
          DUNE_THROW(Dune::Exception, "Simulation state file " << filename << " could not be found!");

        debug_info << idStr << "[MoriNernstPlanckParameters::init()] Using file '" << filename
           << "' for loading boundary data!" << std::endl;

        double t = -1;

        std::string line;
        std::getline(con_in, line);
        int nLine = 0;

        debug_jochen << idStr << "Searching for initial time value " << time << std::endl;

        // Search for time step
        while(t < time && std::abs(t-time) > epsTime && con_in.good())
        {
          std::getline(con_in, line);
          nLine++;

          size_t pos = line.find("# time:");
          if(pos == std::string::npos)
            continue;
          else
            debug_jochen << idStr << "Found '# time' token in line " << (nLine+1) << std::endl;

          std::stringstream line_str(line.substr(8));
          line_str >> t;
        }
        // We reached the end of the file while reading; assume there is a newline at the end of the file
        if(! con_in.good())
        {
          nLinesPerTimeStep = nLine-1;
        } else {
          // We found a line starting a new time value block; substract 3 (current line + 2 newlines)
          nLinesPerTimeStep = nLine-3;
        }
        debug_jochen << idStr << "Number of lines per time step: " << nLinesPerTimeStep << std::endl;
        boundaryValues_old.resize(nLinesPerTimeStep);
        boundaryValues_new.resize(nLinesPerTimeStep);

        // Go back to the beginning
        con_in.clear();
        con_in.seekg(0, std::ios::beg);
      }
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
          //debug_jochen << "I[" << i << "][" << j << "] = " << I[i][j] << std::endl;
        }
      }

      //debug_verb << "[NernstPlanckParameters] @ " << e.geometry().global(x) << ", A = " << I << std::endl;
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

      //debug_jochen << ION_NAMES[ionSpecies_] << " diff coeff = " << D << std::endl;

      int z = physics.getValence(ionSpecies);

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

        //debug_jochen << "### " << physics.getElementIndex(e) << " " << firstCorner
        //             << " <= " << injectionPosition << " <= " << lastCorner << std::endl;

        // Inject sodium into this element
        if(elementContainsInjectionPosition && ionSpecies == Na)
        {
          double con_inj = Iinj;
          con_inj /= geo.volume();

          // Bring to the right units (using microsecond time scale)
          con_inj *= physics.getTimeScale();

          //debug_jochen << "### x = "<< xglobal << ", volume = " << geo.volume() << std::endl;
          //debug_jochen << "### Increasing " << ION_NAMES[ionSpecies_]
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

      //debug_verb << "Nernst-Planck BOUNDARY @ x = " << is.geometry().global(x) << ": "
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
  //          debug_jochen << "[membInterface] "<< "x = " << is.geometry().global(x)
  //              << " => " << ION_NAMES[ionSpecies] << " bctype = " << bctype << std::endl;
          return bctype;
        }
      }

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
              DUNE_THROW(Dune::NotImplemented,  "Boundary types for more than 3 membranes not implemented!");

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
//          debug_jochen << "[top] "<< "x = " << is.geometry().global(x)
//            << " => " << ION_NAMES[ionSpecies] << " bctype = " << bctype << std::endl;
        return bctype;
      }
      // *** Handle BOTTOM boundary ***
      if(bottom)
      {
        bctype = (bottomDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                  : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
//          debug_jochen << "[bottom] "<< "x = " << is.geometry().global(x)
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
//          debug_jochen << "[left] "<< "x = " << is.geometry().global(x)
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
//          debug_jochen << "[right] "<< "x = " << is.geometry().global(x)
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
      typename Traits::RangeFieldType y(0.0);

      typename Traits::DomainType xglobal = e.geometry().global(x);
      bool top = xglobal[1] + 1e-6 > params.yMax();
      bool bottom = xglobal[1] - 1e-6 < params.yMin();
      bool left = xglobal[0] - 1e-6 < params.xMin();
      bool right = xglobal[0] + 1e-6 > params.xMax();
      // TODO Extend this to the case of multiple membranes / membrane element layers
      bool membrane = xglobal[1] + 1e-6 > params.yMemb()[0] && xglobal[1] - 1e-6 < params.yMemb()[0] + params.dMemb();

      bool boundary = top || bottom || left || right || membrane;

      if (useLoadedData && isLoadBoundaryDirichlet)
      {
        if ((bottom && loadBoundaryLocation == "bottom")
            || (membrane && loadBoundaryLocation == "membrane")
            || (top && loadBoundaryLocation == "top"))
        {
          // TODO Use elementMap in the same way as in j()

          double x_load = -1;
          double y_load = -1;
          double value = 0;

          bool isNewTime = std::abs(time-timeNew) < std::abs(time-timeOld);
          const std::vector<std::tuple<double,double,double> >& boundaryValues =
              (isNewTime ? boundaryValues_new : boundaryValues_old);

          // Search for x coordinate
          int i = 0;
          while (x_load < xglobal[0]) {
            x_load = std::get<0>(boundaryValues[i]);
            y_load = std::get<1>(boundaryValues[i]);
            value = std::get<2>(boundaryValues[i]);
            i++;
          }
          debug_jochen << idStr << "  At x=" << xglobal[0] << ", found old coordinate x=" << x_load
              << " -> value: " << value << std::endl;

          // Make sure we have found a matching x coordinate
          if (std::abs(x_load - xglobal[0]) > 1e-6)
            DUNE_THROW(Dune::Exception,
                "Could not find boundary value for x coordinate " << xglobal[0]
                 << ", found x=" << x_load << " instead (difference: " << std::abs(x_load-xglobal[0]) << ")!");

          y = value;

          // convert [mV] -> [1] (not necessary when loading from file)
          //y *= con_e / (1.0e3 * con_k * 279.45);

          return y;
        }
      }

      //int elemIndex = physics.getElementIndex(e);
      int subdomainIndex = physics.getSubdomainIndex(e);

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
      typename INITIAL_GF::Traits::RangeType yCustom(y);
      initialGF.evaluateGlobal(xglobal,yCustom);

      return yCustom[ionSpecies];
    }


    /*! \brief Neumann boundary condition
     */
    typename Traits::RangeFieldType
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        bool includeMoriFlux = true,
        bool implicit = false,
        const PotRangeType& uPot = PotRangeType(-1.0),
        const ConRangeType& uCon = ConRangeType(-1.0)) const
    {
      typedef typename Traits::IntersectionType::Geometry GEO_ORIG;

      typename Traits::RangeFieldType y(0.0);

      // Questionable condition to decide if we are at timeOld or timeNew: Check if difference from
      // time to timeNew is smaller than to timeOld; this surely only works if time is assumed to be
      // either timeOld or timeNew, not some time stage value in between!
      bool isNewTime = !(time > 0.0) || std::abs(time-timeNew) < std::abs(time-timeOld);

      if(! isNewTime)
      {
        DUNE_THROW(Dune::Exception, "Trying to evaluate Neumann boundary condition at old time "
            << time << ", this is not implemented for HH membrane flux! [timeOld = "
            << timeOld << ", timeNew = " << timeNew << "]");
      }

      if (useLoadedData && !isLoadBoundaryDirichlet)
      {
        typename Traits::DomainType xglobal = is.geometry().global(x);
        bool top = xglobal[1]+1e-6 > params.yMax();
        bool bottom = xglobal[1]-1e-6 < params.yMin();
        bool left = xglobal[0]-1e-6 < params.xMin();
        bool right = xglobal[0]+1e-6 > params.xMax();

        bool membrane = physics.isMembraneInterface(is);

        bool membraneActive = true;
        // Use config file flag 'isActive' to decide which of the (possible more than one) membranes is
        // assigned the loaded current density source. Membranes which are not active do not get any boundary source!
        if(membrane)
        {
          const std::vector<bool> isActive = params.isMembraneActive();
          int membNumber = physics.getMembraneNumber(is);
          membraneActive = isActive[membNumber];
        }

        if ((bottom && loadBoundaryLocation == "bottom")
            || (membrane && loadBoundaryLocation == "membrane" && membraneActive)
            || (top && loadBoundaryLocation == "top"))
        {
          double x_load = -1;
          double y_load = -1;
          double value = 0;

          // This is deactivated, see above: evaluating flux at old time is not currently supported
          //const std::vector<std::tuple<double,double,double> >& boundaryValues =
          //    (isNewTime ? boundaryValues_new : boundaryValues_old);
          const std::vector<std::tuple<double,double,double> >& boundaryValues = boundaryValues_new;

          typename PHYSICS::ElementPointer ep = is.inside();

//          if(physics.isMembrane(*ep))
//            DUNE_THROW(Dune::Exception, "Found membrane element although membrane contributions are switched off!");

          // Use inside element as identifier
          int elemIndex = physics.getElementIndex(*ep);
          if(! isElementMapSetUp)
          {
            if(elementMap.count(elemIndex) == 0)
            {
              int offset = elemIndex % physics.nElements(0);
              debug_jochen << idStr << "  Searching value for element #" << elemIndex << " (offset " << offset << ")" << std::endl;

              double x_load_prev = -1;
              double y_load_prev = -1;
              double value_prev = 0;

              // Search for x coordinate
              int i = 0;
              while (x_load < xglobal[0] && i<boundaryValues.size())
              {
                // Store old values
                x_load_prev = x_load;
                y_load_prev = y_load;
                value_prev = value;
                // Read new values
                x_load = std::get<0>(boundaryValues[i]);
                y_load = std::get<1>(boundaryValues[i]);
                value = std::get<2>(boundaryValues[i]);
                i++;
              }
              i--;
              // Since Neumann values are needed at quadrature points and not at the original (cell center)
              // coordinates, we need to find the nearest cell center for the given quadrature point at xglobal.
              // Put the values at the closest cell center into (x_load,y_load,value)
              if(x_load_prev > -1 && std::abs(x_load_prev-xglobal[0]) < std::abs(x_load-xglobal[0]))
              {
                x_load = x_load_prev;
                y_load = y_load_prev;
                value =  value_prev;
                i--;
              }

              debug_verb << idStr << "  At x=" << xglobal
                  << " (inside subdomain:" << physics.getSubDomainNumber(*is.inside()) << ", group '"
                  << physics.getGroupName(*is.inside())
                  << "'), found old x-coordinate " << x_load
                  << " @ index " << i << " -> value: " << value << std::endl;
              debug_verb << idStr << "  Corresponding flux in [mol / (m^2 s)]: "
                  << (value * physics.getLengthScale() / physics.getTimeScale()) << std::endl;

              // Store the index in the boundary vector corresponding to this element index
              elementMap[elemIndex] = i;
            } else {
              int i = elementMap.at(elemIndex);
              x_load = std::get<0>(boundaryValues[i]);
              y_load = std::get<1>(boundaryValues[i]);
              value = std::get<2>(boundaryValues[i]);
            }
          } else {
            // Throw exception if there is no entry in the elementMap for the current element
            if(elementMap.count(elemIndex) == 0)
              DUNE_THROW(Dune::Exception, "No entry in elementMap for elemIndex #" << elemIndex
                  << ", position " << xglobal);

            int i = elementMap.at(elemIndex);
            x_load = std::get<0>(boundaryValues[i]);
            y_load = std::get<1>(boundaryValues[i]);
            value = std::get<2>(boundaryValues[i]);
          }

          // Make sure we have found a matching x coordinate (distance should be <= 0.5 * x-diameter of cell
          const typename PHYSICS::Element& eInside = *is.inside();
          typename Traits::DomainType diam = eInside.geometry().corner(eInside.geometry().corners()-1);
          diam -= eInside.geometry().corner(0);
          if (x_load < 0.0 || std::abs(x_load - xglobal[0]) > diam[0])
            DUNE_THROW(Dune::Exception,
                "Could not find matching boundary value for x coordinate " << xglobal[0]
                 << ", found x=" << x_load << " instead (difference: " << std::abs(x_load-xglobal[0])
                 << ", cell diameter: " << diam[0] << ")!");

          //debug_jochen << idStr << " raw value: " << value << std::endl;

          y = value;

          // I'm not quite sure why the loaded flux has the wrong sign, because I didn't have to do that in Laplace
          // setup if I remember correctly, but then, wasn't the whole generated potential flipped and I had to
          // correct that in Matlab? Anyway, since the membrane flux output is generated from a membrane element,
          // but here we are coming from a cytosol or ES element, the factor -1 is certainly needed!
          y *= -1;

          /* Hack to use a line source instead of a cylinder mantle source: Divide by cylinder area volume
           * (to compensate for the integration over the intersection done in the operator) and multiply
           * by the plain 1D intersection volume (i.e., the length of the line)
           */
          if(collapseToLineSource)
          {
            typedef typename Traits::IntersectionType::Geometry GEO_ORIG;
            typedef typename Acme2CylGeometrySwitch::GeometrySwitch<GEO_ORIG>::type CYL_GEO;
            CYL_GEO cylGeo(is.geometry());

            typename Traits::DomainFieldType originalRadius = params.boundary.get("originalRadius",505.);

            // 2*pi*r * length = surface of the original cylinder segment
            typename Traits::DomainFieldType originalSurface = 2 * con_pi * originalRadius * is.geometry().volume();

            // Now scale this: Multiply by original surface (--> current!) and divide by
            // the membrane surface in the present geometry (--> again current density!)
            // (the calculated factor actually is originalRadius/currentRadius)
            //typename Traits::DomainFieldType surfaceScale = (originalSurface / cylGeo.volume());

            // New version; cylinder geometry trafo switched off in operator, so multiply density by original
            // surface (--> current) and divide by length of the line (i.e., intersection)
            // This enables us to use ymin = 0 to compare with LSA!
            typename Traits::DomainFieldType surfaceScale = (originalSurface / is.geometry().volume());

            if(time < 20.0)
              debug_jochen << "Collapsing to 1D line source by using factor " << surfaceScale << std::endl;

            y *= surfaceScale;
          }

          // Update: Don't forget the conductivity!
          //debug_verb << "   NOT Using conductivity " << physics.getPermittivity(*ep) << std::endl;
          //y /= physics.getPermittivity(*ep);

          // convert [mV] -> [1] (not necessary when loading from file)
          //y *= con_e / (1.0e3 * con_k * 279.45);

          // Multiply by -1 if we are on the opposite side of the membrane
          if(!collapseToLineSource && std::abs(xglobal[1]-y_load) > 1e-6)
          {
            //debug_jochen << "-- flip sign 1" << std::endl;
            y *= -1;
            // Use only the primary membrane side to impose Neumann BC, set other side to zero
            if(params.boundary.get("membraneOneSided",false))
            {
              y = 0.0;
            }

            // This is getting really confusing... surface-scaling is done for membrane interfaces later on,
            // but we also want to do this if we are on a load boundary, but the y-coordinate does not match,
            // which essentially means that we load a value from the other side of the membrane.
            if(!membrane)
            {
              /* Now the intersection geometry can be wrapped, so that we have access to the
               * real cylinder surface
               */
              typename Acme2CylGeometrySwitch::GeometrySwitch<GEO_ORIG>::type igeo(is.geometry());

              // Reconstruct original volume (2*pi*r_orig * plain 1D intersection volume)
              RF originalVolume = 2 * con_pi * y_load * is.geometry().volume();

              RF surfaceScaleFactor = originalVolume / igeo.volume();

              y *= surfaceScaleFactor;
            }
          }
          // Membrane normal vector is oriented differently, flip sign!
          if(physics.isMembrane(*is.inside()))
          {
            //debug_jochen << "-- flip sign 2" << std::endl;
            y *= -1;
          }

          //debug_verb << idStr << "   x = " << xglobal << ", loaded value y=" << y << std::endl;

        }
      } else {

        // Not loading flux data from file: Simply evaluate membrane flux GF
        if(physics.isMembraneInterface(is))
        {
          //        // FIXME Just a first try to enable fully-implicit solution including membrane flux update
          //        Dune::FieldVector<double,1> membPot(-2.7);
          //        Dune::FieldVector<double,3> conCytosol, conExtra;
          //        conCytosol[0] = 12; conCytosol[0] = 100;
          //        conCytosol[1] = 125; conCytosol[1] = 4;
          //        conCytosol[2] = 137; conCytosol[2] = 104;
          //        gfMembFlux.updateChannel(is,time,dt,membPot,conCytosol,conExtra);
          //        // End

          typename GF_MEMB_FLUX::Traits::RangeType membFlux;
          gfMembFlux.evaluate(is,x,membFlux,GF_MEMB_FLUX::FluxTerms::TOTAL_FLUX,ionSpecies,implicit,uPot,uCon);

          y = membFlux[ionSpecies];
        }

      }

      // Evaluation of Mori flux is independent of loading of boundary values in the default case
      if(physics.isMembraneInterface(is))
      {
        // Using values from an inactive Mori flux class is dangerous, since they will not necessarily be 0!
        includeMoriFlux = includeMoriFlux && gfMoriFlux.active();
        if(includeMoriFlux)
        {
          // Mori flux is calculated and added to total flux
          typename GF_MORI_FLUX::Traits::RangeType moriFlux;
          gfMoriFlux.evaluate(is,x,moriFlux,ionSpecies,implicit,uPot,uCon);

          y += moriFlux[ionSpecies];
        }

        /* Now the intersection geometry can be wrapped, so that we have access to the
         * real cylinder surface
         */
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

          //debug_jochen << "### " << physics.getElementIndex(e) << " " << firstCorner
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

        RF surfaceScaleFactor = ext_surface / my_surface;

//        debug_jochen << "BC @ [" << MD_SUBDOMAINS[physics.getSubDomainNumber(*is.inside())]
//          << "] element "<< is.inside()->geometry().center() << ": "
//          << ION_NAMES[ionSpecies] << " normal flux = " << y
//          << " (surface-scaled: " << (surfaceScaleFactor * y) << ")" << std::endl;

        y *= surfaceScaleFactor;

        // Return trans-membrane flux for this ion species, scaled by membrane surfaces
        // => Integration over intersection in local operator will yield
        // flux = ext_surface * membFlux[ionSpecies_], the 'real' membrane current!
        return y;
      }

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
      debug_jochen << "# gfMembFlux.updateState" << std::endl;
      gfMembFlux.updateState();
      if(!implicit)
      {
        debug_jochen << "# gfMembFlux.updateFlux" << std::endl;
        gfMembFlux.updateFlux();
      }

      // Update Mori charge layer for this Newton iteration
      debug_jochen << "# gfMoriFlux.updateState" << std::endl;
      gfMoriFlux.updateState();
      if(!implicit)
      {
        debug_jochen << "# gfMoriFlux.updateFlux" << std::endl;
        gfMoriFlux.updateFlux();
      }
    }

    const PHYSICS& getPhysics() const
    {
     return physics;
    }

    std::string getName() const
    {
     return std::string("MoriNernstPlanckParameters " + idStr);
    }

    const GV& getGridView() const
    {
      return gv;
    }

    void prepareNextTimeStep(double nextTime)
    {
      timeOld = time;
      timeNew = nextTime;

      if(useLoadedData)
      {
        boundaryValues_old = boundaryValues_new;

        // New method: Search for next time step in provided file
        assert(con_in.good());
        //debug_jochen << "Trying to read timeNew=" << timeNew << "@ x=" << x[0] << ", skipping first "
        //    << startLine << " lines!" << std::endl;

        // Set flag to indicate elementMap has been set up (happens after first time step)
        if(!isElementMapSetUp && time > 0.0 && elementMap.size() > 0)
        {
          debug_info << idStr << "=== [MoriNernstPlanckParameters::loadNextTimeStep(" << nextTime
              << ")] Setting flag 'isElementMapSetUp'"
              << " (for elementMap.size() = " << elementMap.size() << ", time = " << time << ")" << std::endl;
          isElementMapSetUp = true;
        }

        debug_info << idStr << " @ " << this << " Loading '" << loadBoundaryLocation
            << "' boundary data from file for requested time " << timeNew << "..." << std::endl;

        double t = -1;

        std::string line;

        // Search for time step
        while(t < timeNew && std::abs(t-timeNew) > epsTime && con_in.good())
        {
         std::getline(con_in, line);

         size_t pos = line.find("# time:");
         if(pos == std::string::npos)
           continue;

         std::stringstream line_str(line.substr(8));
         line_str >> t;
        }
        debug_jochen << idStr << "At requested time " << timeNew << ", found old timestep t=" << t << std::endl;

        // Make sure we have found a matching timestep
        if(std::abs(t-timeNew) > epsTime)
         DUNE_THROW(Dune::Exception, "Could not find boundary value for time "
             << timeNew  << ", found t=" << t << " instead!");

        // Cache all boundary entries belonging to the current timestep for quicker loading
        double x_load = -1;
        double y_load = -1;
        double value = 0;
        for(int i=0; i<nLinesPerTimeStep && con_in.good(); i++)
        {
          std::getline(con_in, line);

          std::stringstream line_str(line);
          line_str >> x_load >> y_load >> value;
          //debug_jochen << idStr  << " Processing line " << line << std::endl;
          //debug_verb << idStr << " Found old coordinate x=" << x_load << " (value " << value << ")" << std::endl;

          // TODO Optimization: Only push those coordinates that are actually present on this processor!
          boundaryValues_new[i] = std::tuple<double,double,double>(x_load,y_load,value);
        }
      }
    }

private:
    const int ionSpecies;

    const GV& gv;
    PHYSICS& physics;
    const Acme2CylParameters& params;
    INITIAL_GF& initialGF;
    GF_MEMB_FLUX& gfMembFlux;
    GF_MORI_FLUX& gfMoriFlux;

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

    const bool useLoadedData;
    std::string loadBoundaryLocation;
    const bool isLoadBoundaryDirichlet;
    const bool collapseToLineSource;
    std::string filename;
    std::ifstream con_in;
    int nLinesPerTimeStep;
    std::vector<std::tuple<double,double,double> > boundaryValues_old;
    std::vector<std::tuple<double,double,double> > boundaryValues_new;
    double timeOld;
    double timeNew;

    int startLine;
    mutable bool isElementMapSetUp;
    mutable std::map<int,int> elementMap;
    double epsTime;
    const std::string idStr;

};


#endif /* DUNE_AX1_ACME2CYL_MORI_NERNST_PLANCK_PARAMETERS_HH */
