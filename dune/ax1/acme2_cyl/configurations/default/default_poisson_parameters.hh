/*
 * poisson_parameters.hh
 *
 *  Created on: Sep 9, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_DEFAULT_POISSON_PARAMETERS_HH
#define DUNE_AX1_ACME2CYL_DEFAULT_POISSON_PARAMETERS_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

template<typename GV, typename RF, typename PHYSICS>
class DefaultPoissonParameters
{
  public:
    typedef PHYSICS Physics;

    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    typedef Dune::FieldVector<typename Traits::RangeFieldType,1> PotRangeType;
    typedef Dune::FieldVector<typename Traits::RangeFieldType,NUMBER_OF_SPECIES> ConRangeType;

    enum SideBoundary { Cytosol = 0, Debye = 1, Extracellular = 2};

    DefaultPoissonParameters(const GV& gv_, PHYSICS& physics_)
      : gv(gv_),
        physics(physics_),
        params(physics_.getParams()),
        poissonConstant(physics.getPoissonConstant()),
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
        initPotIn(physics.convertFrom_mV(params.solution_in.get("pot", 0.0)))
    {}

    //! \brief New constructor for setting poissonConstant manually
    DefaultPoissonParameters(const GV& gv_, PHYSICS& physics_, RF poissonConstant_)
      : gv(gv_),
        physics(physics_),
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
        useMembrane(params.useMembrane()),
        xmin(params.xMin()),
        xmax(params.xMax()),
        ymin(params.yMin()),
        ymax(params.yMax()),
        nMembranes(params.nMembranes()),
        initPotEx(physics.convertFrom_mV(params.solution_ex.get("pot", 0.0))),
        initPotIn(physics.convertFrom_mV(params.solution_in.get("pot", 0.0)))
    {
      debug_jochen << "initPotIn = " << initPotIn << std::endl;
      debug_jochen << "initPotEx = " << initPotEx << std::endl;
    }

    //! tensor diffusion coefficient
    typename Traits::PermTensorType
    A (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        const ConRangeType& conc = ConRangeType(-1.0)) const
    {
      typename Traits::PermTensorType I;

      int elemIndex = physics.getElementIndex(e);

      for (std::size_t i=0; i<Traits::dimDomain; i++)
      {
        for (std::size_t j=0; j<Traits::dimDomain; j++)
        {
          I[i][j] = (i==j) ? -physics.getPermittivity(e) : 0;
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
        const ConRangeType& conc = ConRangeType(-1.0)) const
    {
      bool missingConc = (conc[0] == -1.0);

      // Missing concentration parameter 'conc' is only allowed when we are on the membrane or when we
      // use the Laplace operator (no concentrations present!)
      assert(!missingConc || (missingConc && physics.isMembrane(e)));

      typename Traits::RangeFieldType rhs(0.0);
      if(missingConc) return rhs;

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

    //! partial derivative of source term with respect to ion species j
    typename Traits::RangeFieldType
    dfdu (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        const int j, const PotRangeType& phi) const
    {
      typename Traits::RangeFieldType rhs = -physics.getValence(j) * phi * poissonConstant;

      return rhs;
    }

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
      bool isMembraneInterface = physics.isMembraneInterface(is);

      //debug_verb << "Poisson BOUNDARY @ x = " << is.geometry().global(x) << ": "
      //    << poissonBoundary.bctype(is, x, time, isMembraneInterface) << std::endl;

      BCType bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;

      // Fix potential at all boundaries when there is no membrane
      if(not useMembrane)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
      }

      if(isMembraneInterface)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::None;
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
      int subdomainIndex = physics.getSubdomainIndex(e);

      typename Traits::RangeFieldType y(0.0);

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
                xMemb = (xglobal[1] - yMemb[k]) / dMemb;
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
      return 0.0;
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
      return std::string("PoissonParameters");
    }

    // dummy method
    void prepareNextTimeStep(double nextTime)
    {
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


};

#endif /* DUNE_AX1_ACME2CYL_DEFAULT_POISSON_PARAMETERS_HH */
