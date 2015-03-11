/*
 * poisson_parameters.hh
 *
 *  Created on: Sep 9, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_LAPLACE_PARAMETERS_HH
#define DUNE_AX1_ACME2CYL_LAPLACE_PARAMETERS_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

template<typename GV, typename RF, typename PHYSICS>
class LaplaceParameters
{
  public:

    enum SideBoundary { Cytosol = 0, Debye = 1, Extracellular = 2};

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    LaplaceParameters(PHYSICS& physics_)
      : physics(physics_),
        params(physics_.getParams()),
        poissonConstant(physics_.getPoissonConstant()),
        topDirichlet(params.isTopBoundaryDirichlet_Potential()),
        bottomDirichlet(params.isBottomBoundaryDirichlet_Potential()),
        leftCytosolDirichlet(params.isLeftCytosolBoundaryDirichlet_Potential()),
        rightCytosolDirichlet(params.isRightCytosolBoundaryDirichlet_Potential()),
        leftDebyeDirichlet(params.isLeftDebyeBoundaryDirichlet_Potential()),
        rightDebyeDirichlet(params.isRightDebyeBoundaryDirichlet_Potential()),
        leftExtracellularDirichlet(params.isLeftExtracellularBoundaryDirichlet_Potential()),
        rightExtracellularDirichlet(params.isRightExtracellularBoundaryDirichlet_Potential()),
        membraneNeumann(params.isMembraneBoundaryNeumann_Potential()),
        yMemb(params.yMemb()),
        dMemb(params.dMemb()),
        debyeLayerWidth(params.boundary.get("debyeLayerWidth",0.)),
        useLoadedData(params.boundary.get("useTimeDependentBoundaryValues",false)),
        loadBoundaryLocation(params.boundary.get("loadBoundary",std::string("bottom"))),
        isLoadBoundaryDirichlet(params.isBoundaryDirichlet_Potential(loadBoundaryLocation)),
        collapseToLineSource(params.boundary.get("collapseToLineSource",false)),
        filename(params.boundary.get("boundaryLoadFilename","boundary_pot_bottom.dat")),
        pot_in(filename),
        nLinesPerTimeStep(0),
        boundaryValues(),
        startLine(0),
        isElementMapSetUp(false),
        time(0.0),
        epsTime(params.boundary.get("boundaryTimePrecision",1e-4))
    {
      init();
    }

    //! \brief New constructor for setting poissonConstant manually
    LaplaceParameters(PHYSICS& physics_, RF poissonConstant_)
      : physics(physics_),
        params(physics_.getParams()),
        poissonConstant(poissonConstant_),
        topDirichlet(params.isTopBoundaryDirichlet_Potential()),
        bottomDirichlet(params.isBottomBoundaryDirichlet_Potential()),
        leftCytosolDirichlet(params.isLeftCytosolBoundaryDirichlet_Potential()),
        rightCytosolDirichlet(params.isRightCytosolBoundaryDirichlet_Potential()),
        leftDebyeDirichlet(params.isLeftDebyeBoundaryDirichlet_Potential()),
        rightDebyeDirichlet(params.isRightDebyeBoundaryDirichlet_Potential()),
        leftExtracellularDirichlet(params.isLeftExtracellularBoundaryDirichlet_Potential()),
        rightExtracellularDirichlet(params.isRightExtracellularBoundaryDirichlet_Potential()),
        membraneNeumann(params.isMembraneBoundaryNeumann_Potential()),
        yMemb(params.yMemb()),
        dMemb(params.dMemb()),
        debyeLayerWidth(params.boundary.get("debyeLayerWidth",0.)),
        useLoadedData(params.boundary.get("useTimeDependentBoundaryValues",false)),
        loadBoundaryLocation(params.boundary.get("loadBoundary",std::string("bottom"))),
        isLoadBoundaryDirichlet(params.isBoundaryDirichlet_Potential(loadBoundaryLocation)),
        collapseToLineSource(params.boundary.get("collapseToLineSource",false)),
        filename(params.boundary.get("boundaryLoadFilename","boundary_pot_bottom.dat")),
        pot_in(filename),
        nLinesPerTimeStep(0),
        boundaryValues(),
        startLine(0),
        isElementMapSetUp(false),
        time(0.0),
        epsTime(params.boundary.get("boundaryTimePrecision",1e-4))
    {
      init();
      std::cout << "Collapse to line source: " << collapseToLineSource << std::endl;
    }

    void init()
    {
      if(useLoadedData)
      {
        if(loadBoundaryLocation != "top" && loadBoundaryLocation != "bottom" && loadBoundaryLocation != "membrane")
          DUNE_THROW(Dune::NotImplemented, "Loading left and right boundaries currently not implemented!");

        if(! pot_in.good())
          DUNE_THROW(Dune::Exception, "Simulation state file " << filename << " could not be found!");

        debug_info << "[LaplaceParameters::init()] Using file '" << filename
            << "' for loading boundary data!" << std::endl;

        double t = -1;

        std::string line;
        std::getline(pot_in, line);
        int nLine = 0;

        debug_jochen << "Searching for initial time value " << time << std::endl;

        // Search for time step
        while(t < time && std::abs(t-time) > epsTime && pot_in.good())
        {
         std::getline(pot_in, line);
         nLine++;

         size_t pos = line.find("# time:");
         if(pos == std::string::npos)
           continue;
         else
           debug_jochen << "Found '# time' token in line " << (nLine+1) << std::endl;

         std::stringstream line_str(line.substr(8));
         line_str >> t;
        }
        // We reached the end of the file while reading; assume there is a newline at the end of the file
        if(! pot_in.good())
        {
          nLinesPerTimeStep = nLine-1;
        } else {
          // We found a line starting a new time value block; substract 3 (current line + 2 newlines)
          nLinesPerTimeStep = nLine-3;
        }
        debug_jochen << "Number of lines per time step: " << nLinesPerTimeStep << std::endl;
        boundaryValues.resize(nLinesPerTimeStep);

        // Go back to the beginning
        pot_in.clear();
        pot_in.seekg(0, std::ios::beg);
      }
    }

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
          // Permittivity is conductivity => Scale with length scale to get units [S/nm]
          I[i][j] = (i==j) ? (-physics.getPermittivity(e) * physics.getLengthScale()) : 0;
          //I[i][j] = (i==j) ? -physics.getPermittivity(e) : 0;
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
      return typename Traits::RangeType(0.0);
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
      assert(conc.size() == 0);

      typename Traits::RangeFieldType rhs(0.0);
      if(conc.size() == 0) return rhs;

      typename Traits::RangeFieldType chargeDensity(0.0);
      for(int j=0; j<conc.size(); j++)
      {
        chargeDensity += physics.getValence(j) * conc[j];
      }

      rhs = -chargeDensity * poissonConstant;

      //debug_verb << "[PoissonParameters] chargeDensity = " << chargeDensity << std::endl;
      //debug_verb << "[PoissonParameters] @ " << e.geometry().global(x) << ", rhs = " << rhs << std::endl;

      return rhs;
    }

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      BCType bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;

      bool isMembraneInterface = physics.isMembraneInterface(is);

      //debug_verb << "Poisson BOUNDARY @ x = " << is.geometry().global(x) << ": "
      //    << poissonBoundary.bctype(is, x, time, isMembraneInterface) << std::endl;

      if(isMembraneInterface)
      {
        if(membraneNeumann)
          return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
        else
          return Dune::PDELab::ConvectionDiffusionBoundaryConditions::None;
      }

      bool top = is.geometry().global(x)[1]+1e-6 > params.yMax();
      bool bottom = is.geometry().global(x)[1]-1e-6 < params.yMin();
      bool left = is.geometry().global(x)[0]-1e-6 < params.xMin();
      bool right = is.geometry().global(x)[0]+1e-6 > params.xMax();

      SideBoundary sideBoundary = SideBoundary::Extracellular;

      if(params.useMembrane() && (left || right))
      {
        typename Traits::DomainFieldType yglobal = is.geometry().global(x)[1];
        switch(params.nMembranes())
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
            if(params.nMembranes() < 1 || params.nMembranes() > 3)
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

      // TOP
      if(top)
      {
        bctype = (topDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                               : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        return bctype;
      }
      // BOTTOM
      if(bottom)
      {
        bctype = (bottomDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                  : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        return bctype;
      }
      // *** Handle LEFT boundary ***
      if(left)
      {
        if(sideBoundary == SideBoundary::Cytosol)
        {
          bctype = (leftCytosolDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                         : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
          //debug_jochen << "Left cytosol @ " << is.geometry().global(x) << ", bctype: " << bctype << std::endl;
        } else if (sideBoundary == SideBoundary::Debye) {
          bctype = (leftDebyeDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                       : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
          //debug_jochen << "Left debye @ " << is.geometry().global(x) << ", bctype: " << bctype << std::endl;
        } else {
          bctype = (leftExtracellularDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                               : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
          //debug_jochen << "Left extra @ " << is.geometry().global(x) << ", bctype: " << bctype << std::endl;
        }

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
        return bctype;
      }

      debug_jochen << "x = " << is.geometry().global(x) << " => POT bctype = " << bctype << std::endl;
      //return bctype;
      DUNE_THROW(Dune::Exception, "Error checking this boundary's position!");
    }

    //! Dirichlet boundary condition value / initial value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      typename Traits::DomainType xglobal = e.geometry().global(x);
      bool top = xglobal[1] + 1e-6 > params.yMax();
      bool bottom = xglobal[1] - 1e-6 < params.yMin();
      bool left = xglobal[0] - 1e-6 < params.xMin();
      bool right = xglobal[0] + 1e-6 > params.xMax();
      // TODO Extend this to the case of multiple membranes / membrane element layers
      bool membrane = xglobal[1] + 1e-6 > params.yMemb()[0] && xglobal[1] - 1e-6 < params.yMemb()[0] + params.dMemb();

      bool boundary = top || bottom || left || right || membrane;

      typename Traits::RangeFieldType y(0.0);

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
          // Search for x coordinate
          int i = 0;
          while (x_load < xglobal[0]) {
            x_load = std::get<0>(boundaryValues[i]);
            y_load = std::get<1>(boundaryValues[i]);
            value = std::get<2>(boundaryValues[i]);
            i++;
          }
          debug_jochen << "  At x=" << xglobal[0] << ", found old coordinate x=" << x_load
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

      } else {

        // Deactivated
        if(false)
        {
          // Use a smooth Gaussian function
          const double baseline = -65;
          const double a = 110; // amplitude
          const double x0 = -2e6; // negative start coordinate of traveling signal
          const double v = 1e3; // velocity in nm/µs, here equal to 1m/s
          const double b = x0 + time * v; // center of the peak
          const double c = 1e6; // standard deviation (in nm)

          // evaluate Gaussian function
          double exponent = -(xglobal[0] - b) * (xglobal[0] - b) / (2 * c * c);
          y = baseline + a * std::exp(exponent);

          // convert [mV] -> [1]
          y *= con_e / (1.0e3 * con_k * 279.45);
          return y;
        }
      }

      // Default: return parameter file value
      if(physics.getGroupIndex(e) == ES)
      {
        y = physics.convertFrom_mV(params.solution_ex.get("pot", 0.0));
        return y;
      }
      if(physics.getGroupIndex(e) == CYTOSOL)
      {
        y = physics.convertFrom_mV(params.solution_in.get("pot", -65.0));
        //debug_jochen << "g= " << y << std::endl;
        return y;
      }
      if(time > 0.0)
      {
        DUNE_THROW(Dune::Exception, "Element is neither ES nor CYTOSOL!");
      }

      y = physics.convertFrom_mV(params.solution_in.get("pot", 0.0));
      return y;
    }

    //! Neumann boundary condition
    typename Traits::RangeFieldType
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      typename Traits::DomainType xglobal = is.geometry().global(x);
      bool top = xglobal[1]+1e-6 > params.yMax();
      bool bottom = xglobal[1]-1e-6 < params.yMin();
      bool left = xglobal[0]-1e-6 < params.xMin();
      bool right = xglobal[0]+1e-6 > params.xMax();

      bool membrane = physics.isMembraneInterface(is);
      // TODO Extend this to the case of multiple membranes / membrane element layers
      //bool membrane = xglobal[1] + 1e-6 > params.yMemb()[0] && xglobal[1] - 1e-6 < params.yMemb()[0] + params.dMemb();

//      debug_verb << "  xglobal=" << xglobal << ", top=" << top << ", bottom=" << bottom << ", left=" << left << ", right=" << right
//          << ", membrane=" << membrane << std::endl;

      typename Traits::RangeFieldType y(0.0);

      if (useLoadedData && !isLoadBoundaryDirichlet)
      {
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

          typename PHYSICS::ElementPointer ep = is.inside();

//          if(physics.isMembrane(*ep))
//            DUNE_THROW(Dune::Exception, "Found membrane element although membrane contributions are switched off!");



          // Use inside element as identifier
          int elemIndex = physics.getElementIndex(*ep);
          if(! isElementMapSetUp)
          {
            int offset = elemIndex % physics.nElements(0);
            //debug_jochen << "  Searching value for element #" << elemIndex << " (offset " << offset << ")" << std::endl;

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

            debug_verb << "  At x=" << xglobal
                << " (inside subdomain:" << physics.getSubDomainNumber(*is.inside()) << ", group '"
                << physics.getGroupName(*is.inside())
                << "'), found old x-coordinate " << x_load
                << " @ index " << i << " -> value: " << value << std::endl;
            debug_verb << "  Corresponding flux in [mol / (m^2 s)]: "
                << (value * physics.getLengthScale() / physics.getTimeScale()) << std::endl;

            // Store the index in the boundary vector corresponding to this element index
            elementMap[elemIndex] = i;
          } else {
            int i = elementMap.at(elemIndex);
            x_load = std::get<0>(boundaryValues[i]);
            y_load = std::get<1>(boundaryValues[i]);
            value = std::get<2>(boundaryValues[i]);
          }

          // Make sure we have found a matching x coordinate (distance should be <= 0.5 * x-diameter of cell
          typename Traits::DomainType diam = is.geometry().corner(is.geometry().corners()-1);
          diam -= is.geometry().corner(0);
          if (x_load < 0.0 || std::abs(x_load - xglobal[0]) > diam[0])
            DUNE_THROW(Dune::Exception,
                "Could not find matching boundary value for x coordinate " << xglobal[0]
                 << ", found x=" << x_load << " instead (difference: " << std::abs(x_load-xglobal[0]) << ")!");

          // Multiply by Poisson constant
          y = value * poissonConstant;

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
          debug_verb << "   NOT Using conductivity " << physics.getPermittivity(*ep) << std::endl;
          //y /= physics.getPermittivity(*ep);

          // convert [mV] -> [1] (not necessary when loading from file)
          //y *= con_e / (1.0e3 * con_k * 279.45);

          // Multiply by -1 if we are on the opposite side of the membrane
          if(!collapseToLineSource && std::abs(xglobal[1]-y_load) > 1e-6)
          {
            y *= -1;
            // Use only the primary membrane side to impose Neumann BC, set other side to zero
            if(params.boundary.get("membraneOneSided",false))
            {
              y = 0.0;
            }
          }
          // Membrane normal vector is oriented differently, flip sign!
          if(physics.isMembrane(*is.inside()))
          {
            y *= -1;
          }
          debug_verb << "   Returning value " << y << std::endl;

          // TODO Maybe do surface scaling as in NernstPlanckParameters
          return y;
        }
      }

      // Default: return zero flux condition
      return y;
    }

    //! outflow boundary condition
    typename Traits::RangeFieldType
    o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      return 0.0;
    }

    //! set time
    void setTime (double t)
    {
      time = t;

      // New method: Search for next time step in provided file
      if(useLoadedData)
      {
        assert(pot_in.good());
        //debug_jochen << "Trying to read time=" << time << "@ x=" << x[0] << ", skipping first "
        //    << startLine << " lines!" << std::endl;

        // Set flag to indicate elementMap has been set up (happens after first time step)
        if(!isElementMapSetUp && elementMap.size() > 0)
          isElementMapSetUp = true;

        debug_info << "= Loading '" << loadBoundaryLocation << "' boundary data from file for requested time "
            << time << "..." << std::endl;

        double t = -1;

        std::string line;

        // Search for time step
        while(t < time && std::abs(t-time) > epsTime && pot_in.good())
        {
         std::getline(pot_in, line);

         size_t pos = line.find("# time:");
         if(pos == std::string::npos)
           continue;

         std::stringstream line_str(line.substr(8));
         line_str >> t;
        }
        debug_jochen << "At requested time " << time << ", found old timestep t=" << t << std::endl;

        // Make sure we have found a matching timestep
        if(std::abs(t-time) > epsTime)
         DUNE_THROW(Dune::Exception, "Could not find boundary value for time "
             << time  << ", found t=" << t << " instead!");

        // Cache all boundary entries belonging to the current timestep for quicker loading
        double x_load = -1;
        double y_load = -1;
        double value = 0;
        for(int i=0; i<nLinesPerTimeStep && pot_in.good(); i++)
        {
          std::getline(pot_in, line);

          std::stringstream line_str(line);
          line_str >> x_load >> y_load >> value;
          //debug_jochen << " Processing line " << line << std::endl;
          debug_verb << " Found old coordinate x=" << x_load << std::endl;

          // TODO Optimization: Only push those coordinates that are actually present on this processor!
          boundaryValues[i] = std::tuple<double,double,double>(x_load,y_load,value);
        }
      }
    }

    const PHYSICS& getPhysics() const
    {
      return physics;
    }

    std::string getName() const
    {
      return std::string("LaplaceParameters");
    }

    double getTime() const
    {
      return time;
    }

    bool doCollapseToLineSource() const
    {
      return collapseToLineSource;
    }

    PHYSICS& physics;
    const Acme2CylParameters& params;
    RF poissonConstant;
    const bool topDirichlet;
    const bool bottomDirichlet;
    const bool leftCytosolDirichlet;
    const bool rightCytosolDirichlet;
    const bool leftDebyeDirichlet;
    const bool rightDebyeDirichlet;
    const bool leftExtracellularDirichlet;
    const bool rightExtracellularDirichlet;
    const bool membraneNeumann;
    const std::vector<double> yMemb;
    const double dMemb;
    const double debyeLayerWidth;

    const bool useLoadedData;
    std::string loadBoundaryLocation;
    const bool isLoadBoundaryDirichlet;
    const bool collapseToLineSource;
    std::string filename;
    std::ifstream pot_in;
    int nLinesPerTimeStep;
    std::vector<std::tuple<double,double,double> > boundaryValues;

    int startLine;
    mutable bool isElementMapSetUp;
    mutable std::map<int,int> elementMap;

    double time;
    double epsTime;
};

#endif /* DUNE_AX1_ACME2CYL_POISSON_PARAMETERS_HH */
