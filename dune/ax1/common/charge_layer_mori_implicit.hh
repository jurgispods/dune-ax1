/*
 * charge_layer_mori_implicit.hh
 *
 *  Created on: Dec 3, 2014
 *      Author: jpods
 */

#ifndef DUNE_AX1_CHARGE_LAYER_MORI_IMPLICIT_HH
#define DUNE_AX1_CHARGE_LAYER_MORI_IMPLICIT_HH


template<typename DGF_CON, typename DGF_POT, typename PHYSICS, typename SolutionContainer>
class ChargeLayerMoriImplicit
  : public Dune::PDELab::IntersectionGridFunctionBase<
      Dune::PDELab::IntersectionGridFunctionTraits<typename DGF_CON::Traits::GridViewType,
                                       typename DGF_CON::Traits::RangeFieldType,
                                       DGF_CON::Traits::RangeType::dimension,
                                       typename DGF_CON::Traits::RangeType>,
     ChargeLayerMori<DGF_CON,DGF_POT,PHYSICS,SolutionContainer> >
{
public:
  typedef Dune::PDELab::IntersectionGridFunctionTraits<typename DGF_CON::Traits::GridViewType,   // grid view type
                                           typename DGF_CON::Traits::RangeFieldType, // image field type (double)
                                           DGF_CON::Traits::RangeType::dimension,  // number of components of image
                                           typename DGF_CON::Traits::RangeType // image type (Dune::FieldVector<double, 1>)
                                           > Traits;
    typedef typename PHYSICS::FieldType T;
    typedef std::vector<T> VPOT;
    typedef Dune::FieldVector<T, NUMBER_OF_SPECIES> VCON;
    typedef std::vector<VCON> V;

  public:

    enum MembraneInterface { INTERNAL = 1, EXTERNAL = 2};

    //! constructor
    ChargeLayerMoriImplicit (DGF_CON& dgfCon_, DGF_POT& dgfPot_, PHYSICS& physics_,
        SolutionContainer& solutionContainer_, double tEquilibrium_ = 20000)
    : physics(physics_),
      solutionContainer(solutionContainer_),
      tEquilibrium(tEquilibrium_),
      gridView(dgfCon_.getGridView()),
      isActive(physics.getParams().useMori()),
      useMoriCorrection(physics.getParams().useMoriCorrection()),
      isInitialized(false),
      partiallyInitialized(false),
      timeStepAccepted(true),
      membraneIndices(physics.getChannelSet().getMembraneIndices()),
      dMemb(physics.getParams().dMemb() * physics.getLengthScale()),
      temp(physics.getElectrolyte().getTemperature()),
      stdCon(physics.getElectrolyte().getStdCon()),
      epsElec(physics.getElectrolyte().getPermittivity()),
      time(0.0),
      dt(1.0),
      tau(1e-9), // Time is measured in units [s] => 1 ns
      gamma_ext_old(physics.getChannelSet().getMembraneIndices().size(), VCON(0.0)),
      gamma_ext_new(physics.getChannelSet().getMembraneIndices().size(), VCON(0.0)),
      gamma_int_old(physics.getChannelSet().getMembraneIndices().size(), VCON(0.0)),
      gamma_int_new(physics.getChannelSet().getMembraneIndices().size(), VCON(0.0)),
      useMoriOperatorSplit(physics.getParams().general.get("useMoriOperatorSplit",false)),
      iterateMembraneFlux(physics.getParams().boundary.get("fullyImplicitMembraneFlux",false)
          || useMoriOperatorSplit),
      surfaceChargeScale((physics.getTimeScale() / physics.getLengthScale()) / (con_e * con_mol)),
      verboseMembraneFluxClasses(physics.getParams().boundary.get("verboseMembraneFluxClasses",false))
    {
      assert(DGF_CON::Traits::RangeType::dimension == NUMBER_OF_SPECIES);
      if(!iterateMembraneFlux)
        DUNE_THROW(Dune::Exception, "This application uses the implicit version of membrane flux classes, but the "
            << "corresponding flag 'fullyImplicitMembraneFlux' in config file is not set correctly!");
    }

    template<typename I>
    inline void evaluate (const I &is,
                          const typename Traits::DomainType &x,
                          typename Traits::RangeType &y) const
    {
      evaluate(is, x, y, -1, false);
    }

    template<typename I>
    inline void evaluate (const I &is,
                          const typename Traits::DomainType &x,
                          typename Traits::RangeType &y,
                          int ionSpecies,
                          bool implicit,
                          const typename DGF_POT::Traits::RangeType& uPot = typename DGF_POT::Traits::RangeType(-1.0),
                          const typename DGF_CON::Traits::RangeType& uCon = typename DGF_CON::Traits::RangeType(-1.0)) const
    {
      y = 0.0;

      if(isActive)
      {
        if(implicit)
          assert(uCon[0] > -1.0);

        assert(time >= 0.0 && dt > 0.0);

        // Deactivate output
        if(! implicit || !verboseMembraneFluxClasses)
          debug_verb.push(false);

        const std::vector<bool> membraneActive = physics.getParams().isMembraneActive();

//        debug_verb << "membraneActive: [ ";
//        for(int i=0; i<membraneActive.size(); ++i)
//        {
//          debug_verb << membraneActive[i] << " ";
//        }
//        debug_verb << "]" << std::endl;

        const double real_dt = dt * physics.getTimeScale();

        // We can now use the newtonSolution's pointer to the current solution vector
        // instead of the solution from the last time step.
        DGF_POT dgfPotNew(solutionContainer.getGfsPot(), *solutionContainer.getSolutionPotNew());
        // We also need the old (last time step) potential to calculate the time derivative
        DGF_POT dgfPotOld(solutionContainer.getGfsPot(), *solutionContainer.getSolutionPotOld());

        //typename DGF_CON::BASE_DGF dgfConElecNew(solutionContainer.getGfsCon(), *solutionContainer.getSolutionCon());
        //DGF_CON dgfConNew(gridView, dgfConElecNew, physics);


        if(physics.isMembraneInterface(is))
        {
          // Determine on which side of a membrane interface we are!
          int membraneInterface = 0;
          // Now we need to distinguish between cytosol and extracellular elements!
          switch(physics.getSubdomainIndex(*is.inside()))
          {
            case CYTOSOL:
            {
              membraneInterface = MembraneInterface::INTERNAL;
              break;
            }
            case ES:
            {
              membraneInterface = MembraneInterface::EXTERNAL;
              break;
            }
            default:
              // We come from an inside membrane element; check outside element to determine if this is and interior
              // (cytosol) or exterioi (ES) interface!
              switch(physics.getSubdomainIndex(*is.outside()))
              {
                case CYTOSOL:
                {
                  membraneInterface = MembraneInterface::INTERNAL;
                  break;
                }
                case ES:
                {
                  membraneInterface = MembraneInterface::EXTERNAL;
                  break;
                }
                default:
                  DUNE_THROW(Dune::Exception, "Membrane flux can only be evaluated on a real membrane interface, not on"
                      << " membrane-membrane intersections!");
              }
          }

          int mIndex = physics.getMembraneIndex(is);
          int membNumber = physics.getMembraneNumber(is);

//          // TODO Debug, remove
          if(implicit && mIndex != 0)
            debug_verb.push(false);

          // Do nothing when membrane is passive
          if(membraneActive[membNumber])
          {
            // Attention! We are coming from a non-primary (not even normal-oriented) intersection;
            // get the primary intersection in order to get the correct indices!
            typename PHYSICS::ElementIntersectionIterator miit = physics.getPrimaryMembraneIntersection(is);
            int iIndex = physics.getIntersectionIndex(*miit);
            assert(checkIndex(iIndex));

            int mIndex = membraneIndices.at(iIndex);

            // Open channels in middle membrane element only
            if(true /*mIndex == (membGV.size(0)/2 - 1)*/)
            {
              debug_verb << std::endl << "====================";
              debug_verb << "Calculating Mori flux for membrane element #" << iIndex << "[" << mIndex << "]"
                    << " (memb #" << membNumber << ") ====================" << std::endl;

              // get membrane permittivity
              typename PHYSICS::FieldType perm = physics.getMembranePermittivity(is);

              // We need the membrane potential for calculation of Mori flux
              typename DGF_POT::Traits::RangeType potJump;
              // Now check if current unknown values have been handed over by operator; if not, we use the old values
              // from the solution container
              if(implicit && uCon[0] == -1)
                  DUNE_THROW(Dune::Exception, "When using the implicit version of evaluate, you must provide current values "
                      << "for both potential and concentrations, even though concentrations are not being used!");

              if(implicit)
              {
                // Implicit case: Use special version of method handing over current value of unknown
                physics.getMembranePotentialJump(is, dgfPotNew, potJump, uPot);
              } else {
                physics.getMembranePotentialJump(is, dgfPotNew, potJump);
              }
              typename DGF_POT::Traits::RangeType potJump_mV = physics.convertTo_mV(potJump);

              typename DGF_POT::Traits::RangeType potJumpOld;
              physics.getMembranePotentialJump(is, dgfPotOld, potJumpOld);
              typename DGF_POT::Traits::RangeType potJumpOld_mV = physics.convertTo_mV(potJumpOld);

              // Update flux
              typename Traits::RangeType totalMoriFlux_ext(0.0);
              typename Traits::RangeType totalMoriFlux_int(0.0);

              int jLow = ionSpecies;
              int jUp = ionSpecies;
              if(ionSpecies == -1)
              {
                jLow = 0;
                jUp = NUMBER_OF_SPECIES-1;
              }
              for(int j=jLow; j<=jUp; j++)
              {
                debug_verb << "-------------------------------------" << std::endl;

                if(useMoriCorrection)
                {
                  DUNE_THROW(Dune::NotImplemented, "Mori concentration profile correction not implemented for the implicit case!");
                }

                /* The surface charge density returned by charge layer class has units [C / m^2],
                 * the (approximated) time derivative has units [C / (s * m^2)]. To match our membrane flux
                 * units [(TS/LS) * mol / (m^2 * s)], we need to multiply by scaling factor
                 * surfaceChargeScale = (TS/LS) / (e*N_A)!
                 */
                debug_verb << " Calculating Mori flux for perm = " << perm
                    << ", real_dt = " << real_dt << " [implicit=" << implicit << "]" << std::endl;


                // capacitance per unit area (F / m^2)
                const T C_m = con_eps0 * perm / dMemb;
                const T membPotOld = 1e-3 * potJumpOld_mV;
                const T membPotNew = 1e-3 * potJump_mV;

                T sigma_old = 0.0;
                T sigma_new = 0.0;
                if(membraneInterface == MembraneInterface::INTERNAL)
                {
                  sigma_old = gamma_int_old[mIndex][j] * C_m * membPotOld;
                  sigma_new = gamma_int_new[mIndex][j] * C_m * membPotNew;

                  Dune::ios_base_all_saver haselnussSchnaps(std::cout);
                  debug_verb << "  [int] old sigma = " << gamma_int_old[mIndex][j] << " * " << C_m
                      << " * " << membPotOld <<" = " << sigma_old << std::endl;
                  debug_verb << "  [int] new sigma = " << gamma_int_new[mIndex][j] << " * " << C_m
                      << " * " << membPotNew <<" = " << sigma_new << std::endl;
                  debug_verb << std::scientific << "  gamma diff: "  << (gamma_int_new[mIndex][j] - gamma_int_old[mIndex][j])
                      << ", pot diff: " << std::scientific << (membPotNew - membPotOld) << std::endl;

                } else {
                  sigma_old = gamma_ext_old[mIndex][j] * C_m * membPotOld;
                  sigma_new = gamma_ext_new[mIndex][j] * C_m * membPotNew;

                  Dune::ios_base_all_saver haselnussSchnaps(std::cout);
                  debug_verb << "  [ext] old sigma = " << gamma_ext_old[mIndex][j] << " * " << C_m
                      << " * " << membPotOld <<" = " << sigma_old << std::endl;
                  debug_verb << "  [ext] new sigma = " << gamma_ext_new[mIndex][j] << " * " << C_m
                      << " * " << membPotNew <<" = " << sigma_new << std::endl;
                  debug_verb << std::scientific << "  gamma diff: "  << (gamma_ext_new[mIndex][j] - gamma_ext_old[mIndex][j])
                      << ", pot diff: " << std::scientific << (membPotNew - membPotOld) << std::endl;
                }

                const T moriFlux = surfaceChargeScale * (sigma_new - sigma_old) / real_dt;

                debug_verb << " => Mori [" << ION_NAMES[j] << "] flux = " << moriFlux
                    << (membraneInterface == MembraneInterface::INTERNAL ? " (CY)" : " (ES)") << std::endl;

                y[j] = moriFlux;

                debug_verb << "-------------------------------------" << std::endl;
              }
              debug_verb << "-------------------------------------" << std::endl;
            }
          }

//          // TODO Debug, remove
          if(implicit && mIndex != 0)
            debug_verb.pop();

          // 2D: Multiply by y-component of outer normal, as the membrane is oriented horizontally
          y *= is.centerUnitOuterNormal()[Traits::GridViewType::Grid::dimension - 1];
        }

        // Reactivate output
        if(! implicit || !verboseMembraneFluxClasses)
          debug_verb.pop();
      }
    }

    inline const typename Traits::GridViewType& getGridView () const
    {
      return gridView;
    }

    inline void updateState(double time_, double dt_)
    {
      time = time_;
      dt = dt_;

      updateState();
    }

    inline void updateState(const typename SolutionContainer::U& u)
    {
      DUNE_THROW(Dune::NotImplemented, "updateState(u) currently not supported. Use updateState() instead and make"
        << " sure that solutionContainer contains the desired solution vector(s)!");
      assert(u == *solutionContainer.getSolutionConNew());
      updateState();
    }

    // TODO Revamp output of this method (membrane info, make output optional, etc.)
    inline void updateState()
    {
      if(isActive)
      {
        assert(time >= 0.0 && dt > 0.0);

        // Prevent multiple channel initializations within one timestep!
        if(! isInitialized && partiallyInitialized)
          return;

        timeStepAccepted = false;

        // We can now use the newtonSolution's pointer to the current solution vector
        // instead of the solution from the last time step.
        typename DGF_CON::BASE_DGF dgfConElecNew(solutionContainer.getGfsCon(), *solutionContainer.getSolutionConNew());
        DGF_CON dgfConNew(gridView, dgfConElecNew, physics);

        const std::vector<bool> membraneActive = physics.getParams().isMembraneActive();
        for(typename PHYSICS::MIterator mit = physics.mBegin(); mit != physics.mEnd(); ++mit)
        {
          int membNumber = physics.getMembraneNumber(*mit);
          if(membraneActive[membNumber])
          {
            int iIndex = physics.getIntersectionIndex(*mit);

            assert(checkIndex(iIndex));
            int mIndex = membraneIndices.at(iIndex);

            debug_verb << std::endl << "====================";
            debug_verb << " Updating Mori states for membrane element #" << iIndex << "[" << mIndex << "]"
                  << " (memb #" << membNumber << ") ====================" << std::endl;

            typename DGF_CON::Traits::RangeType conCytosol, conExtra;
            physics.getMembraneConcentrationJump(*mit, dgfConNew, conCytosol, conExtra);

            // Update channels (i.e., gating particles)
            if(not isInitialized)
            {
              // Init active channels once at time == tEquilibrium
              if(not physics.getParams().doEquilibration() ||
                  (std::abs(time-tEquilibrium) < 1e-6 || time > tEquilibrium)) // time >= tEquilibrium
              {
                // Init channels
                debug_info << "Initializing charge layer to CY/ES concentrations of " << conCytosol
                    << " / " << conExtra << std::endl;
                init(mIndex, conCytosol, conExtra);

                partiallyInitialized = true;
              }
            } else {
              // Do one time step for charge layer
              double real_dt = dt * physics.getTimeScale();

              timeStep(mIndex, real_dt, conCytosol, conExtra);
            }
            debug_verb << "-------------------------------------" << std::endl;
          }
        }
      }
    }

    inline void updateFlux(double time_, double dt_)
    {
      DUNE_THROW(Dune::NotImplemented, "Method updateFlux() should never be called in implicit flux calculation!");
    }

    inline void updateFlux(const typename SolutionContainer::U& u)
    {
      DUNE_THROW(Dune::NotImplemented, "Method updateFlux() should never be called in implicit flux calculation!");
    }

    inline void updateFlux()
    {
      DUNE_THROW(Dune::NotImplemented, "Method updateFlux() should never be called in implicit flux calculation!");
    }

    bool isPartiallyInitialized() const
    {
      return partiallyInitialized;
    }

    bool isCompletelyInitialized() const
    {
      return isInitialized;
    }

    inline void acceptTimeStep()
    {
      timeStepAccepted = true;

      if(partiallyInitialized)
      {
        isInitialized = true;

        debug_info << std::endl;
        debug_info << "============================" << std::endl;
        debug_info << "^^ Mori charge layer initialized ^^" << std::endl;
        debug_info << "============================" << std::endl;
        debug_info << std::endl;

        partiallyInitialized = false;
      }

      // Assign current values of state variables to be used as the old state
      gamma_ext_old = gamma_ext_new;
      gamma_int_old = gamma_int_new;
    }

    inline void discardTimeStep()
    {
      timeStepAccepted = false;

      // Set time, dt to default values
      time = -1;
      dt = -1;
    }

    void setActive(bool active)
    {
      isActive = active;
    }

    bool active() const
    {
      return isActive;
    }

    //! \brief return exterior surface charge density with units [C / m^2]
    template<typename IS>
    T getNewSurfaceChargeDensity_Ext(const IS& is, const int j, const T perm, const T memb_pot) const
    {
      int m = physics.getMembraneIndex(is);

//      typename PHYSICS::FieldType perm = 0.0;
//      if(physics.isMembrane(*is.inside()))
//      {
//        perm = physics.getPermittivity(*is.inside());
//      } else {
//        perm = physics.getPermittivity(*is.outside());
//      }

      // capacitance per unit area (F / m^2) (Mori 2006: C*_m)
      T C_m = con_eps0 * perm / dMemb;

      // Bring membrane potential to units of [V]
      return (gamma_ext_new[m][j] * C_m * 1e-3 * memb_pot);
    }

    //! \brief return interior surface charge density with units [C / m^2]
    template<typename IS>
    T getNewSurfaceChargeDensity_Int(const IS& is, const int j, const T perm, const T memb_pot) const
    {
      int m = physics.getMembraneIndex(is);

//      typename PHYSICS::FieldType perm = 0.0;
//      if(physics.isMembrane(*is.inside()))
//      {
//        perm = physics.getPermittivity(*is.inside());
//      } else {
//        perm = physics.getPermittivity(*is.outside());
//      }

      // capacitance per unit area (F / m^2)
      T C_m = con_eps0 * perm / dMemb;
      // Bring membrane potential to units of [V]
      return (gamma_int_new[m][j] * C_m * 1e-3 * memb_pot);
    }

  private:
    // Calculate gamma_snake values and use as initial values
    void init(int m, const VCON& conCytosol, const VCON& conExtra)
    {
      // Calculate new gamma_snake values
      T sum_ext(0.0);
      T sum_int(0.0);

      for(int j=0; j<conCytosol.size(); j++)
      {
        const T z2 = physics.getValence(j)*physics.getValence(j);
        sum_ext += z2*conExtra[j];
        sum_int += z2*conCytosol[j];
      }

      VCON gamma_sum_ext(0.0);
      VCON gamma_sum_int(0.0);

      debug_info << "Initializing charge layer fractions to: " << std::endl;
      for(int j=0; j<conCytosol.size(); j++)
      {
        const T z2 = physics.getValence(j)*physics.getValence(j);
        gamma_ext_old[m][j] = z2 * conExtra[j] / sum_ext;
        gamma_int_old[m][j] = z2 * conCytosol[j] / sum_int;
        gamma_ext_new[m][j] = gamma_ext_old[m][j];
        gamma_int_new[m][j] = gamma_int_old[m][j];

        gamma_sum_ext += gamma_ext_old[m][j];
        gamma_sum_int += gamma_int_old[m][j];

        debug_info << "gamma_ext[" << j << "] = " << gamma_ext_old[m][j] << " / "
          << "gamma_int[" << j << "] = " << gamma_int_old[m][j] << std::endl;
      }
      debug_info << "SUM: " << gamma_sum_ext << " / " << gamma_sum_int << std::endl;
    }

    //! \brief Do one (BE) time step for membrane element i
    //! \param dt Real dt [s]
    //! \param v Potential [mV]
    void timeStep(int m, const T dt_, const VCON& conCytosol, const VCON& conExtra)
    {
      // Calculate new gamma_snake values
      VCON gamma_ext_snake(0.0);
      VCON gamma_int_snake(0.0);

      T sum_ext(0.0);
      T sum_int(0.0);

      for(int j=0; j<conCytosol.size(); j++)
      {
        const T z2 = physics.getValence(j)*physics.getValence(j);
        sum_ext += z2*conExtra[j];
        sum_int += z2*conCytosol[j];
      }

      for(int j=0; j<conCytosol.size(); j++)
      {
        const T z2 = physics.getValence(j)*physics.getValence(j);
        gamma_ext_snake[j] = z2 * conExtra[j] / sum_ext;
        gamma_int_snake[j] = z2 * conCytosol[j] / sum_int;
      }

      debug_verb << "  Charge layer timestep (dt =" << dt_ << ")" << std::endl;
      // Update gamma values
      for(int j=0; j<conCytosol.size(); j++)
      {
        gamma_ext_new[m][j] = (tau * gamma_ext_old[m][j] + dt_ * gamma_ext_snake[j]) / (tau + dt_);
        gamma_int_new[m][j] = (tau * gamma_int_old[m][j] + dt_ * gamma_int_snake[j]) / (tau + dt_);

        debug_verb << "   " << gamma_ext_old[m][j] << " -> " << gamma_ext_new[m][j] << std::endl;
        debug_verb << "   " << gamma_int_old[m][j] << " -> " << gamma_int_new[m][j] << std::endl;
      }

      // Check if sum equals one
      T gamma_sum_ext(0.0);
      T gamma_sum_int(0.0);
      for(int j=0; j<conCytosol.size(); j++)
      {
        gamma_sum_ext += gamma_ext_new[m][j];
        gamma_sum_int += gamma_int_new[m][j];
      }
      if(std::abs(gamma_sum_ext-1.0) > 1e-6 || std::abs(gamma_sum_int-1.0) > 1e-6)
      {
        DUNE_THROW(Dune::Exception, "ChargeLayerMori::timeStep(): Sum of new gamma's is not equal to 1.0! "
            << " (ES: " << gamma_sum_ext << " / CY: " << gamma_sum_int << ")");
      }
    }

    bool checkIndex(const int i) const
    {
      if(membraneIndices.find(i) == membraneIndices.end())
      {
        DUNE_THROW(Dune::Exception, "Given element #" << i << " is not a membrane element!");
      }
      return true;
    }

    PHYSICS& physics;
    SolutionContainer& solutionContainer;
    double tEquilibrium;

    const typename Traits::GridViewType& gridView;
    bool isActive;
    bool useMoriCorrection;

    bool isInitialized;
    bool partiallyInitialized;
    bool timeStepAccepted;

    const std::map<int, int>& membraneIndices;

    const T dMemb;
    const T temp;
    const T stdCon;
    const T epsElec;
    T time;
    T dt;
    const T tau;

    V gamma_ext_old;
    V gamma_ext_new;
    V gamma_int_old;
    V gamma_int_new;

    const bool useMoriOperatorSplit;
    const bool iterateMembraneFlux;

    const T surfaceChargeScale;

    const bool verboseMembraneFluxClasses;
};


#endif /* DUNE_AX1_CHARGE_LAYER_MORI_IMPLICIT_HH */
