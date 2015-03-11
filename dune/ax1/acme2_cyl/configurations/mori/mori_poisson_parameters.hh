/*
 * poisson_parameters.hh
 *
 *  Created on: Sep 9, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_MORI_POISSON_PARAMETERS_HH
#define DUNE_AX1_ACME2CYL_MORI_POISSON_PARAMETERS_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

template<typename GV, typename RF, typename PHYSICS>
class MoriPoissonParameters
{
  public:
    typedef PHYSICS Physics;

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    // Hardcoded range field type for potential and concentrations
    typedef Dune::FieldVector<typename Traits::RangeFieldType,1> PotRangeType;
    typedef Dune::FieldVector<typename Traits::RangeFieldType,NUMBER_OF_SPECIES> ConRangeType;

    enum SideBoundary { Cytosol = 0, Debye = 1, Extracellular = 2};

    //! \brief New constructor for setting poissonConstant manually
    MoriPoissonParameters(const GV& gv_, PHYSICS& physics_, RF poissonConstant_)
      : gv(gv_),
        physics(physics_),
        params(physics_.getParams()),
        poissonConstant(poissonConstant_),
        lengthScale(physics_.getLengthScale()),
        timeScale(physics_.getTimeScale()),
        temp(physics.getElectrolyte().getTemperature()),
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
        useMembrane(params.useMembrane()),
        xmin(params.xMin()),
        xmax(params.xMax()),
        ymin(params.yMin()),
        ymax(params.yMax()),
        nMembranes(params.nMembranes()),
        initPotEx(physics.convertFrom_mV(params.solution_ex.get("pot", 0.0))),
        initPotIn(physics.convertFrom_mV(params.solution_in.get("pot", 0.0))),
        couplePotentialNeumannBoundaryToConcentration(
            params.boundary.get("couplePotentialNeumannBoundaryToConcentration",false)),
        useLoadedData(params.boundary.get("useTimeDependentBoundaryValuesPot",false)),
        loadBoundaryLocation(params.boundary.get("loadBoundary",std::string("bottom"))),
        isLoadBoundaryDirichlet(params.isBoundaryDirichlet_Potential(loadBoundaryLocation)),
        collapseToLineSource(params.boundary.get("collapseToLineSource",false)),
        filename(params.boundary.get("boundaryPotLoadFilename","boundary_pot_bottom.dat")),
        pot_in(filename),
        nLinesPerTimeStep(0),
        boundaryValues(),
        startLine(0),
        isElementMapSetUp(false),
        epsTime(params.boundary.get("boundaryTimePrecision",1e-4))
    {
      debug_verb << "MoriPoissonParameters(gv,physics,poissonConstant)" << std::endl;

      debug_jochen << "initPotIn = " << initPotIn << std::endl;
      debug_jochen << "initPotEx = " << initPotEx << std::endl;

      init();
    }

    void init()
    {
      if(useLoadedData)
      {
        if(loadBoundaryLocation != "top" && loadBoundaryLocation != "bottom" && loadBoundaryLocation != "membrane")
          DUNE_THROW(Dune::NotImplemented, "Loading left and right boundaries currently not implemented!");

        if(! pot_in.good())
          DUNE_THROW(Dune::Exception, "Simulation state file " << filename << " could not be found!");

        debug_info << "[MoriPoissonParameters::init()] Using file '" << filename
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
    A (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        const ConRangeType& conc = ConRangeType(-1.0)) const
    {
      if(conc[0] == -1.0)
        DUNE_THROW(Dune::Exception, "MoriPoissonParameters::A() was called without concentration vector parameter!");

      typename Traits::PermTensorType I;

      int elemIndex = physics.getElementIndex(e);

      typename Traits::RangeFieldType a(0.0);
      for(int j=0; j<conc.size(); j++)
      {
        // This definition of 'a' lacks the factor e/(kT) from Mori paper
        // (we already have a dimensionless potential); but is it correct?
        //a += con_e * con_mol * physics.getValence(j) * physics.getValence(j) *
        //    physics.getDiffCoeff(j, elemIndex) * conc[j];
        // Update: Need to scale by LS to get expressions in units of nanometres...
        a += physics.getLengthScale() * con_e * con_mol * physics.getValence(j) * physics.getValence(j) *
                    physics.getDiffCoeff(j, elemIndex) * conc[j];

        // Version with factor e/(kT) included
        //a += con_e * con_mol * physics.getValence(j) * physics.getValence(j) *
        //      physics.getDiffCoeff(j, elemIndex) * conc[j] * con_e / (con_k * temp);
      }

      for (std::size_t i=0; i<Traits::dimDomain; i++)
      {
        for (std::size_t j=0; j<Traits::dimDomain; j++)
        {
          // Positive: same sign as term 'bGrad' in operator
          I[i][j] = (i==j) ? a : 0;

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

    //! special term 'grad b' appearing in Mori potential equation (dependent on concentration gradient)
    typename Traits::RangeType
    bGrad (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        const std::vector<Dune::FieldVector<RF,GV::dimension> > AgraduCon) const
    {
      // Special bGrad for Mori electroneutral operator
      Dune::FieldVector<RF,GV::dimension> bGrad(0.0);

      // Mori model: Calculate weighted sum of concentration gradients
      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        // Start with reusing AgraduCon (= D_j * gradCon_j)
        Dune::FieldVector<RF,GV::dimension> bGradSingleCon = AgraduCon[j];

        // Multiply by e*N_a*z_j and scale by LS to get expressions in units of nanometres...
        bGradSingleCon *= lengthScale * con_e * con_mol * physics.getValence(j);

        bGrad += bGradSingleCon;
      }

      return bGrad;
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
        std::vector<typename Traits::RangeFieldType> fCon = std::vector<typename Traits::RangeFieldType>()) const
    {
      // Missing concentration parameter 'conc' is only allowed when we are on the membrane or when we
      // use the Laplace operator (no concentrations present!)
      assert(fCon.size() > 0 || (fCon.size() == 0 && physics.isMembrane(e)));

      typename Traits::RangeFieldType f(0.0);

      for(int j=0; j<fCon.size(); j++)
      {
        f += physics.getValence(j) * fCon[j];
      }

      f *= -1 * poissonConstant * (timeScale / lengthScale) * con_e / (con_k * temp);

      // TODO Remove this
      f *= params.general.get("arbFac_fPot",1.0);

      return f;
    }

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      bool isMembraneInterface = physics.isMembraneInterface(is);

      //debug_verb << "Poisson BOUNDARY @ x = " << is.geometry().global(x) << ": "
      //    << poissonBoundary.bctype(is, x, time, isMembraneInterface) << std::endl;

      BCType bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;

      // Do not allow Dirichlet boundary conditions on membrane for now!
      if(isMembraneInterface)
      {
        bctype = (membraneNeumann ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann
                                  : Dune::PDELab::ConvectionDiffusionBoundaryConditions::None);
        return bctype;
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
              DUNE_THROW(Dune::NotImplemented, "Boundary types for more than 3 membranes not implemented!");

            for(int m = 0; m<nMembranes; m++)
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
        } else if (sideBoundary == SideBoundary::Debye) {
          bctype = (leftDebyeDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                       : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
        } else {
          bctype = (leftExtracellularDirichlet ? Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                                               : Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann);
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
      }

      int subdomainIndex = physics.getSubdomainIndex(e);

      // Default initialization: Homogeneous concentrations on each subdomain
      switch(subdomainIndex)
      {
        case CYTOSOL:
        {
          y = initPotIn;
          break;
        }
        case ES:
        {
          y = initPotEx;
          break;
        }
        case MEMBRANE:
        {
          if(time > 0.0)
          {
            DUNE_THROW(Dune::Exception, "Element is neither ES nor CYTOSOL!");
          } else {

            typename Traits::DomainType xglobal = e.geometry().global(x);

            bool isOnMembraneInterface = false;
            double xMemb = -1;
            int membraneNumber = -1;

            for(int k=0; k<yMemb.size(); ++k)
            {
              // Position is located on membrane k
              if(xglobal[1] > yMemb[k]-1e-6 && xglobal[1] < yMemb[k]+dMemb+1e-6)
              {
                // Local coordinate between 0 and 1 with respect to membrane thickness
                xMemb = xglobal[1] - yMemb[k] / dMemb;
                membraneNumber = k;
                break;
              }
            }
            // Assign intrecellular potential to interior membrane positions
            if(xMemb < 0 || xMemb > 1)
              DUNE_THROW(Dune::Exception, "Membrane position (" << xglobal << "): Error assigning "
                  << " local membrane coordinate (" << xMemb << ")!");

            // For now, there is no generic way to determine which side of the membrane is cytosol/ES.
            // Use a hardcoded case-by-case test similar to that for tagging elements in acme2_cyl_par.cc!
            int numMembranes = yMemb.size();
            double y1, y2;
            // Case: Single membrane
            if(numMembranes == 1)
            {
              y1 = initPotIn;
              y2 = initPotEx;
            }
            // Two membranes enclosing cytosol
            if(numMembranes == 2)
            {
              // Tag this cell as extracellular if it is outside the two membrane layers
              if(membraneNumber==0)
              {
                y1 = initPotEx;
                y2 = initPotIn;
              } else {
                y1 = initPotIn;
                y2 = initPotEx;
              }
            }
            // 3 membranes: Extracellular space is between membranes 1 and 2 and above 3
            if(numMembranes == 3)
            {
              if(membraneNumber==0 || membraneNumber==2)
              {
                y1 = initPotIn;
                y2 = initPotEx;
              } else {
                y1 = initPotEx;
                y2 = initPotIn;
              }
            }

            // Now linearly interpolate between membrane interfaces
            y = y1 + xMemb * (y2-y1);

            break;
          }
        }
        default:
          DUNE_THROW(Dune::Exception, "Element has an unknown subdomain index!");
      }

      // TODO If boundary location-dependent Dirichlet values are desired, this has to be implemented!

//      bool top = e.geometry().global(x)[1]+1e-6 > ymax;
//      bool bottom = e.geometry().global(x)[1]-1e-6 < ymin;
//      bool left = e.geometry().global(x)[0]-1e-6 < xmin;
//      bool right = e.geometry().global(x)[0]+1e-6 > xmax;

      //debug_jochen << "g: " << e.geometry().global(x) << " -> " << y << std::endl;

      return y;
    }

    //! Neumann boundary condition
    typename Traits::RangeFieldType
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x,
        std::vector<typename Traits::RangeFieldType> jCon = std::vector<typename Traits::RangeFieldType>()) const
    {
      typename Traits::RangeFieldType y(0.0);

      bool membrane = physics.isMembraneInterface(is);

      if (useLoadedData && !isLoadBoundaryDirichlet)
      {
        typename Traits::DomainType xglobal = is.geometry().global(x);
        bool top = xglobal[1]+1e-6 > params.yMax();
        bool bottom = xglobal[1]-1e-6 < params.yMin();
        bool left = xglobal[0]-1e-6 < params.xMin();
        bool right = xglobal[0]+1e-6 > params.xMax();

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
          //y = value * poissonConstant;
          y = value;

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
          //debug_verb << "   Returning value " << y << std::endl;

          // TODO Maybe do surface scaling as in NernstPlanckParameters
          return y;
        }
      }

      // Couple potential Neumann BC to sum of concentration BC values
      if(membrane || couplePotentialNeumannBoundaryToConcentration)
      {
        for(int j=0; j<jCon.size(); j++)
        {
          // TODO Added weighing by valence here, does this give correct results? Without using Cl- ions,
          // this is equivalent to the previous implementation
          y += physics.getValence(j) * jCon[j];
        }

        // Mori: Multiply by poissonConstant = e * N_A * LS * LS * LS / TS;
        // an additional scaling factor comes from scaling up the diffusion coefficient with units (m^2 / s),
        // so the final scaling factor is equal to (e * N_A * LS)
        y *= poissonConstant * timeScale / (lengthScale * lengthScale);

        //const double arbFac = paramPot.getPhysics().getParams().general.get("arbitraryPotScalingFactor",1.0);
        //y *= arbFac;

        //debug_jochen << "  => coupled jPot value @ " << is.geometry().global(x) << ": " << y << std::endl;
      }

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
    }

    const PHYSICS& getPhysics() const
    {
      return physics;
    }

    const RF getPoissonConstant() const
    {
      return poissonConstant;
    }

    std::string getName() const
    {
      return std::string("MoriPoissonParameters");
    }

    void prepareNextTimeStep(double nextTime)
    {
      if(useLoadedData)
      {
        time = nextTime;

        // New method: Search for next time step in provided file
        assert(pot_in.good());
        //debug_jochen << "Trying to read time=" << time << "@ x=" << x[0] << ", skipping first "
        //    << startLine << " lines!" << std::endl;

        // Set flag to indicate elementMap has been set up (happens after first time step)
        if(!isElementMapSetUp && elementMap.size() > 0)
          isElementMapSetUp = true;

        debug_info << "= [pot] Loading '" << loadBoundaryLocation
           << "' boundary data from file for requested time " << time << "..." << std::endl;

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

    const GV& getGridView() const
    {
      return gv;
    }


  private:
    const GV& gv;
    PHYSICS& physics;
    const Acme2CylParameters& params;
    RF poissonConstant;
    RF lengthScale;
    RF timeScale;
    RF temp;
    double time;

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
    const bool useMembrane;
    const double xmin;
    const double xmax;
    const double ymin;
    const double ymax;
    const int nMembranes;

    RF initPotEx;
    RF initPotIn;

    const bool couplePotentialNeumannBoundaryToConcentration;
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
    double epsTime;

};

#endif /* DUNE_AX1_ACME2CYL_MORI_POISSON_PARAMETERS_HH */
