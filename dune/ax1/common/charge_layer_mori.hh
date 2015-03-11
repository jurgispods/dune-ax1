/*
 * charge_layer_mori.hh
 *
 *  Created on: Jul 30, 2014
 *      Author: jpods
 */

#ifndef DUNE_AX1_CHARGE_LAYER_MORI_HH
#define DUNE_AX1_CHARGE_LAYER_MORI_HH

// TODO This implementation uses the intrinsic membrane capacitance per unit area (C_m)
// instead of the effective membrane capacitance C*_m from Mori and Peskin 2008, which
// also includes capacitance contributions by Debye layers. It still needs to be investigated
// if this makes a considerable difference!

template<typename DGF_CON, typename DGF_POT, typename PHYSICS, typename SolutionContainer>
class ChargeLayerMori
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

    //! constructor
    ChargeLayerMori (DGF_CON& dgfCon_, DGF_POT& dgfPot_, PHYSICS& physics_,
        SolutionContainer& solutionContainer_, double tEquilibrium_ = 20000)
    : dgfCon(dgfCon_),
      dgfPot(dgfPot_),
      physics(physics_),
      solutionContainer(solutionContainer_),
      tEquilibrium(tEquilibrium_),
      gridView(dgfCon_.getGridView()),
      isActive(physics.getParams().useMori()),
      useMoriCorrection(physics.getParams().useMoriCorrection()),
      isInitialized(false),
      partiallyInitialized(false),
      timeStepAccepted(true),
      membraneIndices(physics.getChannelSet().getMembraneIndices()),
      cachedMoriFlux_ext(membraneIndices.size(), typename Traits::RangeType(0.0)),
      cachedMoriFlux_int(membraneIndices.size(), typename Traits::RangeType(0.0)),
      dMemb(physics.getParams().dMemb()),
      temp(physics.getElectrolyte().getTemperature()),
      stdCon(physics.getElectrolyte().getStdCon()),
      epsElec(physics.getElectrolyte().getPermittivity()),
      time(-1),
      dt(-1),
      tau(1e-9), // Time is measured in units [s] => 1 ns
      memb_pot_old(physics.getChannelSet().getMembraneIndices().size(), 0.0),
      memb_pot_new(physics.getChannelSet().getMembraneIndices().size(), 0.0),
      gamma_ext_old(physics.getChannelSet().getMembraneIndices().size(), VCON(0.0)),
      gamma_ext_new(physics.getChannelSet().getMembraneIndices().size(), VCON(0.0)),
      gamma_int_old(physics.getChannelSet().getMembraneIndices().size(), VCON(0.0)),
      gamma_int_new(physics.getChannelSet().getMembraneIndices().size(), VCON(0.0)),
      iterateMembraneFlux(physics.getParams().boundary.get("fullyImplicitMembraneFlux",false)
          || physics.getParams().general.get("useMoriOperatorSplit",false))
    {
      assert(DGF_CON::Traits::RangeType::dimension == NUMBER_OF_SPECIES);
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
      if(! iterateMembraneFlux)
        assert(timeStepAccepted);

      assert(implicit == false);

      y = 0.0;

      if(physics.isMembraneInterface(is))
      {

        int mIndex = physics.getMembraneIndex(is);

        assert(mIndex < cachedMoriFlux_ext.size() && mIndex < cachedMoriFlux_int.size());

        // Now we need to distinguish between cytosol and extracellular elements!
        switch(physics.getSubdomainIndex(*is.inside()))
        {
          case CYTOSOL:
          {
            y = cachedMoriFlux_int[mIndex];
            break;
          }
          case ES:
          {
            y = cachedMoriFlux_ext[mIndex];
            break;
          }
          default:
            // We come from an inside membrane element; check outside element to determine if this is and interior
            // (cytosol) or exterioi (ES) interface!
            switch(physics.getSubdomainIndex(*is.outside()))
            {
              case CYTOSOL:
              {
                y = cachedMoriFlux_int[mIndex];
                break;
              }
              case ES:
              {
                y = cachedMoriFlux_ext[mIndex];
                break;
              }
              default:
                DUNE_THROW(Dune::Exception, "Membrane flux can only be evaluated on a real membrane interface, not on"
                    << " membrane-membrane intersections!");
            }
        }
        // 2D: Multiply by y-component of outer normal, as the membrane is oriented horizontally
        y *= is.centerUnitOuterNormal()[Traits::GridViewType::Grid::dimension - 1];
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
      assert(u == *solutionContainer.getSolutionPot());
      updateState();
    }


    inline void updateState()
    {
      if(isActive)
      {
        assert(time > -1 && dt > -1);

        timeStepAccepted = false;

        // Prevent multiple channel initializations within one timestep!
        if(! isInitialized && partiallyInitialized)
          return;

        // We can now use the newtonSolution's pointer to the current solution vector
        // instead of the solution from the last time step. Implement this for a truly
        // fully-implicit numerical method!
        DGF_POT dgfPotNew(solutionContainer.getGfsPot(), *solutionContainer.getSolutionPot());
        typename DGF_CON::BASE_DGF dgfConElecNew(solutionContainer.getGfsCon(), *solutionContainer.getSolutionCon());
        DGF_CON dgfConNew(dgfCon.getGridView(), dgfConElecNew, physics);

        const std::vector<bool> membraneActive = physics.getParams().isMembraneActive();
        for(typename PHYSICS::MIterator mit = physics.mBegin(); mit != physics.mEnd(); ++mit)
        {
          if(membraneActive[physics.getMembraneNumber(*mit)])
          {
            int iIndex = physics.getIntersectionIndex(*mit);

            assert(checkIndex(iIndex));
            int mIndex = membraneIndices.at(iIndex);

            typename DGF_POT::Traits::RangeType membPot;
            if(iterateMembraneFlux)
            {
              physics.getMembranePotential(*mit, dgfPotNew, membPot);
            } else {
              physics.getMembranePotential(*mit, dgfPot, membPot);
            }
            typename DGF_POT::Traits::RangeType membPot_mV = physics.convertTo_mV(membPot);

            typename DGF_CON::Traits::RangeType conCytosol, conExtra;
            if(iterateMembraneFlux)
            {
              physics.getMembraneConcentrationJump(*mit, dgfConNew, conCytosol, conExtra);
            } else {
              physics.getMembraneConcentrationJump(*mit, dgfCon, conCytosol, conExtra);
            }

            // Update channels (i.e., gating particles)
            if(not isInitialized)
            {
              // Init charge layer once at time == tEquilibrium
              if(not physics.getParams().doEquilibration() ||
                  (std::abs(time-tEquilibrium) < 1e-6 || time > tEquilibrium)) // time >= tEquilibrium
              {
                if(not physics.getParams().general.get("automaticChannelInitialization",true))
                {
                  // Initialize channels to a predefined membrane potential; override default setting
                  membPot_mV = physics.getParams().general.get("channelRestingPotential",membPot_mV);
                }

                // Init channels
                debug_info << "Initializing charge layer to membrane potential of " << membPot_mV << " [mV]" << std::endl;
                init(mIndex, membPot_mV, conCytosol, conExtra);

                partiallyInitialized = true;
              }
            } else {
              // Do one time step for charge layer
              double real_dt = dt * physics.getTimeScale();

              timeStep(mIndex, real_dt, membPot_mV, conCytosol, conExtra);
            }
          }
        }
      }
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

      acceptState();

      // Set time, dt to default values
      time = -1;
      dt = -1;
    }

    inline void updateFlux(double time_, double dt_)
    {
      time = time_;
      dt = dt_;

      updateFlux();
    }

    inline void updateFlux(const typename SolutionContainer::U& u)
    {
      assert(u == *solutionContainer.getSolution());
      updateFlux();
    }

    inline void updateFlux()
    {
      if(isActive)
      {
        assert(time > -1 && dt > -1);

        const std::vector<bool> membraneActive = physics.getParams().isMembraneActive();

        debug_verb << "membraneActive: ";
        Output::printVector(membraneActive);

        const double real_dt = dt * physics.getTimeScale();

        // We can now use the newtonSolution's pointer to the current solution vector
        // instead of the solution from the last time step. Implement this for a truly
        // fully-implicit numerical method!
        DGF_POT dgfPotNew(solutionContainer.getGfsPot(), *solutionContainer.getSolutionPot());
        typename DGF_CON::BASE_DGF dgfConElecNew(solutionContainer.getGfsCon(), *solutionContainer.getSolutionCon());
        DGF_CON dgfConNew(dgfCon.getGridView(), dgfConElecNew, physics);

        // Iterate over all membrane interfaces and calculate flux
        for(typename PHYSICS::MIterator miit = physics.mBegin(); miit != physics.mEnd(); ++miit)
        {
          int membNumber = physics.getMembraneNumber(*miit);
          // Do nothing when membrane is passive
          if(membraneActive[membNumber])
          {
            int iIndex = physics.getIntersectionIndex(*miit);

            assert(checkIndex(iIndex));
            int mIndex = membraneIndices.at(iIndex);

            // Open channels in middle membrane element only
            if(true /*mIndex == (membGV.size(0)/2 - 1)*/)
            {
              debug_verb << std::endl << "====================";
              debug_verb << "Calculating Mori flux for membrane element #" << iIndex << "[" << mIndex << "]"
                    << " (memb #" << membNumber << ") ====================" << std::endl;

              // Get concentrations at the two membrane ends, starting from the membrane element
              typename DGF_CON::Traits::RangeType conRatio;
              if(iterateMembraneFlux)
              {
                physics.getMembraneConcentrationRatio(*miit, dgfConNew, conRatio);
              } else {
                physics.getMembraneConcentrationRatio(*miit, dgfCon, conRatio);
              }

              typename DGF_CON::Traits::RangeType conCytosol;
              typename DGF_CON::Traits::RangeType conExtra;
              if(iterateMembraneFlux)
              {
                physics.getMembraneConcentrationJump(*miit, dgfConNew, conCytosol, conExtra);
              } else {
                physics.getMembraneConcentrationJump(*miit, dgfCon, conCytosol, conExtra);
              }

              // Calc ionic strength
              std::vector<double> ionicStrength(2, 0.0);
              for(int jj=0; jj<conExtra.size(); jj++)
              {
                ionicStrength[0] += 0.5 * conCytosol[jj] * physics.getValence(jj) * physics.getValence(jj);
                ionicStrength[1] += 0.5 * conExtra[jj] * physics.getValence(jj) * physics.getValence(jj);
              }

              std::vector<typename PHYSICS::FieldType> debyeLength(2, 0.0);
              physics.getDebyeLength(debyeLength, ionicStrength);

              // get membrane permittivity
              typename PHYSICS::FieldType perm = physics.getPermittivity(*miit->inside());

              typename DGF_POT::Traits::RangeType potJump;
              if(iterateMembraneFlux)
              {
                physics.getMembranePotentialJump(*miit, dgfPotNew, potJump);
              } else {
                physics.getMembranePotentialJump(*miit, dgfPot, potJump);
              }

              typename DGF_POT::Traits::RangeType potJump_mV = physics.convertTo_mV(potJump);

              // Update flux
              typename Traits::RangeType totalMoriFlux_ext(0.0);
              typename Traits::RangeType totalMoriFlux_int(0.0);

              for(int j=0; j<NUMBER_OF_SPECIES; j++)
              {
                debug_verb << "-------------------------------------" << std::endl;

                int valence = physics.getValence(j);

                // Constant C1 to be reused later
                double C1 = (con_k * temp) / con_e;

                // Calculated concentration-dependent reversal potential calculated by Nernst equation
                double reversal_potential = - std::log(conRatio[j]) * (C1 / valence);

                debug_verb << "[" << ION_NAMES[j] << "] Nernst potential: " << 1e3 * reversal_potential << " [mV]"
                    << std::endl;

                typename Traits::RangeFieldType moriFlux_ext(0.0);
                typename Traits::RangeFieldType moriFlux_int(0.0);


                /* The following block was an experiment based on Mori's concentration/potential gradient reconstruction.
                 * He actually only used this to postprocess his results, but it would be possible to correct the membrane
                 * concentrations/potential by using this reconstruction, and to use those values instead of the actually
                 * calculated unknowns in the flux calculation.
                 * In the current implementation, it just prints some debug output as to how the reconstructed values
                 * _would_ look like, they are not actually used.
                 * It would need to be cleaned up, tested and moved to MembraneFluxGF when actually being used!
                 */
                if(useMoriCorrection)
                {
                  // After all channel contributions, add Mori charge layer contribution; we will also use
                  // Mori's concentration correction to account for the influence of the concentrations
                  // on the Nernst (reversal) potential

                  debug_verb << " Mori charge layer flux contribution " << std::endl;

                  typename Traits::RangeFieldType surf_charge_ext_old =
                      getOldSurfaceChargeDensity_Ext(mIndex, j, perm);
                  typename Traits::RangeFieldType surf_charge_ext_new =
                      getNewSurfaceChargeDensity_Ext(mIndex, j, perm);
                  typename Traits::RangeFieldType surf_charge_int_old =
                      getOldSurfaceChargeDensity_Int(mIndex, j, perm);
                  typename Traits::RangeFieldType surf_charge_int_new =
                      getNewSurfaceChargeDensity_Int(mIndex, j, perm);

                  debug_verb << "  surface charge density ES (old/new): " << surf_charge_ext_old << " / "
                      << surf_charge_ext_new << std::endl;
                  debug_verb << "  surface charge density CY (old/new): "
                      << surf_charge_int_old << " / " << surf_charge_int_new << std::endl;

    //              bool dimensionless = true;
    //              if(dimensionless)
    //              {
    //                double debyeLength_int = std::sqrt( con_eps0 * perm * con_k * temp /
    //                                ( 2.0 * con_e * con_e * stdCon * ionicStrength_int) );
    //                double debyeLength_ext = std::sqrt( con_eps0 * perm * con_k * temp /
    //                                ( 2.0 * con_e * con_e * stdCon * ionicStrength_ext) );
    //
    //                // TODO This is only a rough value from Mori 2006, use real value here!
    //                double c0 = 100;
    //                double debye_length = std::sqrt( con_eps0 * perm * con_k * temp /
    //                    ( con_e * con_e * c0) );
    //                surf_charge_ext_new /= (con_e * c0 * debye_length);
    //                surf_charge_int_new /= (con_e * c0 * debye_length);
    //              }
    //              surf_charge_ext_new /= con_e;
    //              surf_charge_int_new /= con_e;

                  // Try correction from Mori 2006 to calculate membrane interface concentrations that would
                  // have established if we had resolved the Debye layer
                  //double delta_n_ext = surf_charge_ext_new * std::sqrt(valence*valence*conExtra[j]);
                  //double delta_n_int = surf_charge_int_new * std::sqrt(valence*valence*conCytosol[j]);

                  // TODO Calculate real s (Mori 2006: theta*)
                  //double s = 0.0045;

                  // Coupling factor s (see dissertation, chapter 5), equivalent to theta* form Mori 2006
                  double s_int = (perm * debyeLength[0]) / (epsElec * physics.getLengthScale() * dMemb);
                  double s_ext = (perm * debyeLength[1]) / (epsElec * physics.getLengthScale() * dMemb);

                  // Capital Gamma from Mori 2006 (although definition is ambiguous)
                  double Gamma_ext = std::sqrt(2 * ionicStrength[1]);
                  double Gamma_int = std::sqrt(2 * ionicStrength[0]);

                  // Calculate theta (effective membrane capacitance) instead of theta*
                  double theta_int = s_int * Gamma_ext * Gamma_int
                      / (Gamma_ext*Gamma_int + s_int * Gamma_ext + s_int * Gamma_int);
                  double theta_ext = s_ext * Gamma_ext * Gamma_int
                      / (Gamma_ext*Gamma_int + s_ext * Gamma_ext + s_ext * Gamma_int);

                  debug_verb << "perm (memb / elec): " << perm << " / " << epsElec << std::endl;
                  debug_verb << "d_memb: " << dMemb << std::endl;
                  debug_verb << "Debye length (ES / CY): " << debyeLength[1] << " / " << debyeLength[0] << std::endl;
                  debug_verb << "-> theta* (==s) (ES/CY): " << s_ext << " / " << s_int << std::endl;
                  debug_verb << "-> theta        (ES/CY): " << theta_ext << " / " << theta_int << std::endl;

                  double gamma_ext = getNewGamma_Ext(mIndex, j);
                  double gamma_int = getNewGamma_Int(mIndex, j);

                  // Hack in value from compared simulation
                  // plain
                  //double potJump2 = physics.convertFrom_mV(-64.9076);
                  // myelin
                  //double potJump2 = -2.71812;

                  // TODO Use Mori correction also for potential here!?
                  double potJump2 = potJump;


                  // TODO Carefully check this! Should it be sqrt(ionicStrength) or ionicStrength?
                  // It might be that the factor is indeed sqrt(2) * sqrt(100) (2 because in Mori 2006,
                  // factor 2 is missing in the definition of ionic strength, and 100 because of scaling
                  // factor c0 = 100 mM!) => Check this for again!
                  //double delta_n_ext = - gamma_ext * s * potJump * std::sqrt(ionicStrength_ext) / valence;
                  //double delta_n_int = gamma_int * s * potJump * std::sqrt(ionicStrength_int) / valence;
                  //double delta_n_ext = - gamma_ext * s * potJump * std::sqrt(ionicStrength_ext*2*100)/ valence;
                  //double delta_n_int = gamma_int * s * potJump * std::sqrt(ionicStrength_int*2*100) / valence;
                  double delta_n_ext = - gamma_ext * s_ext * potJump2 * Gamma_ext * Gamma_ext / valence;
                  double delta_n_int = gamma_int * s_int * potJump2 * Gamma_int * Gamma_int / valence;
                  double delta_n_ext2 = - gamma_ext * s_ext * potJump2 * Gamma_ext * std::sqrt(2*100) / valence;
                  double delta_n_int2 = gamma_int * s_int * potJump2 * Gamma_int * std::sqrt(2*100) / valence;

                  debug_verb << "  concentration correction (ES): "
                      << delta_n_ext << " -> " << (conExtra[j]+delta_n_ext) << std::endl;
                  debug_verb << "  concentration correction (CY): "
                      << delta_n_int << " -> " << (conCytosol[j]+delta_n_int) << std::endl;
                  debug_verb << "  concentration correction2 (ES): "
                      << delta_n_ext2 << " -> " << (conExtra[j]+delta_n_ext2) << std::endl;
                  debug_verb << "  concentration correction2 (CY): "
                      << delta_n_int2 << " -> " << (conCytosol[j]+delta_n_int2) << std::endl;
                  double new_conRatio = (conCytosol[j]+delta_n_int) / (conExtra[j]+delta_n_ext);
                  double new_nernst = - std::log(new_conRatio) * (con_k * temp) / (con_e * valence);
                  debug_verb << "  -> new Nernst potential E = " << 1e3 * new_nernst << std::endl;

                  reversal_potential = new_nernst;
                  // TODO Actually use new reversal potential calculated above in MembraneFluxGF!
                }

                /* The surface charge density returned by charge layer class has units [C / m^2],
                 * the (approximated) time derivative has units [C / (s * m^2)]. To match our membrane flux
                 * units [(TS/LS) * mol / (m^2 * s)], we need to multiply by scaling factor s = (TS/LS) / (e*N_A)!
                 */
                double scale = (physics.getTimeScale() / physics.getLengthScale()) / (con_e * con_mol);
                debug_jochen << " Calculating Mori flux for perm = " << perm  << ", real_dt = " << real_dt << std::endl;
                moriFlux_ext = scale * getSurfaceChargeDensityTimeDerivative_Ext(mIndex, j, perm, real_dt);
                moriFlux_int = scale * getSurfaceChargeDensityTimeDerivative_Int(mIndex, j, perm, real_dt);

                double cm = con_eps0 * perm / (physics.getParams().dMemb() * physics.getLengthScale());
                debug_jochen << "C_m = " << cm << ", dvdt = " << (moriFlux_ext/cm) << std::endl;

                debug_verb << " => Mori [" << ION_NAMES[j] << "] flux = "<< moriFlux_ext << " (ES), "
                     << moriFlux_int << " (CY)" << std::endl << std::endl;

                totalMoriFlux_ext[j] += moriFlux_ext;
                totalMoriFlux_int[j] += moriFlux_int;

                debug_verb << "-------------------------------------" << std::endl;
                debug_verb << "=> Total [" << ION_NAMES[j] << "] flux = " << totalMoriFlux_ext[j]
                   << " (ES),  " << totalMoriFlux_int[j] << " (CY)" << std::endl;
              }
              debug_verb << "-------------------------------------" << std::endl;
              double sumAllFluxes_ext = 0.0;
              double sumAllFluxes_int = 0.0;
              for(int j=0; j<NUMBER_OF_SPECIES; j++)
              {
                sumAllFluxes_ext += totalMoriFlux_ext[j];
                sumAllFluxes_int += totalMoriFlux_int[j];
              }
              debug_verb << "=> Total sum of all fluxes = " << sumAllFluxes_ext << " (ES),  " <<
                  sumAllFluxes_int << " (CY)" << std::endl;
              debug_verb << "==========================================================" << std::endl << std::endl;

              cachedMoriFlux_ext[mIndex] = totalMoriFlux_ext;
              cachedMoriFlux_int[mIndex] = totalMoriFlux_int;
            }
          }
        }
      }
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
      int mIndex = physics.getMembraneIndex(is);

      // Ignore memb_pot in the explicit implementation; we already have the new memb_pot as a state variable
      return getNewSurfaceChargeDensity_Ext(mIndex, j, perm);
    }

    //! \brief return interior surface charge density with units [C / m^2]
    template<typename IS>
    T getNewSurfaceChargeDensity_Int(const IS& is, const int j, const T perm, const T memb_pot) const
    {
      int mIndex = physics.getMembraneIndex(is);

      // Ignore memb_pot in the explicit implementation; we already have the new memb_pot as a state variable
      return getNewSurfaceChargeDensity_Int(mIndex, j, perm);
    }


  private:

    T getNewGamma_Ext(int m, int j) const
    {
      return gamma_ext_new[m][j];
    }

    T getOldGamma_Ext(int m, int j) const
    {
      return gamma_ext_old[m][j];
    }

    T getNewGamma_Int(int m, int j) const
    {
      return gamma_int_new[m][j];
    }

    T getOldGamma_Int(int m, int j) const
    {
      return gamma_int_old[m][j];
    }

    //! \brief return surface charge density in units [C / m^2]
    //! \param v Potential [mV]
    T getNewSurfaceChargeDensity_Ext(int m, int j, const T perm) const
    {
      // capacitance per unit area (F / m^2) (Mori 2006: C*_m)
      const T C_m = con_eps0 * perm / (dMemb * physics.getLengthScale());
      const T membPot = 1e-3 * memb_pot_new[m];
      const T sigma = (gamma_ext_new[m][j] * C_m * membPot);

      //debug_jochen << "     [getNewSurfaceChargeDensity_Ext] C_m = " << C_m << std::endl;
      //debug_jochen << "     [getNewSurfaceChargeDensity_Ext] memb_pot_new = " << membPot << std::endl;
      //debug_jochen << "     [getNewSurfaceChargeDensity_Ext] => sigma = " << sigma << std::endl;

      // Bring membrane potential to units of [V]
      return sigma;
    }

    //! \param v Potential [mV]
    T getNewSurfaceChargeDensity_Int(int m, int j, const T perm) const
    {
      // capacitance per unit area (F / m^2)
      const T C_m = con_eps0 * perm / (dMemb * physics.getLengthScale());
      // Bring membrane potential to units of [V]
      return (gamma_int_new[m][j] * C_m * 1e-3 * memb_pot_new[m]);
    }

    //! \param v Potential [mV]
    T getOldSurfaceChargeDensity_Ext(int m, int j, const T perm) const
    {
      // capacitance per unit area (F / m^2)
      const T C_m = con_eps0 * perm / (dMemb * physics.getLengthScale());
      const T membPot = 1e-3 * memb_pot_old[m];
      const T sigma = gamma_ext_old[m][j] * C_m * membPot;

      //debug_jochen << "     [getOldSurfaceChargeDensity_Ext] C_m = " << C_m << std::endl;
      //debug_jochen << "     [getOldSurfaceChargeDensity_Ext] memb_pot_old = " << membPot << std::endl;
      //debug_jochen << "     [getOldSurfaceChargeDensity_Ext] => sigma = " << sigma << std::endl;

      // Bring membrane potential to units of [V]
      return sigma;
    }

    //! \param v Potential [mV]
    T getOldSurfaceChargeDensity_Int(int m, int j, const T perm) const
    {
      // capacitance per unit area (F / m^2)
      const T C_m = con_eps0 * perm / (dMemb * physics.getLengthScale());
      // Bring membrane potential to units of [V]
      return (gamma_int_old[m][j] * C_m * 1e-3 * memb_pot_old[m]);
    }


    //! \param v Potential [mV]
    //! \param dt Real dt [s]
    T getSurfaceChargeDensityTimeDerivative_Ext(int m, int j, const T perm, T dt_) const
    {
      return (getNewSurfaceChargeDensity_Ext(m,j,perm) - getOldSurfaceChargeDensity_Ext(m,j,perm)) / dt_;
    }

    //! \param v Potential [mV]
    //! \param dt Real dt [s]
    T getSurfaceChargeDensityTimeDerivative_Int(int m, int j, const T perm, T dt_) const
    {
      return (getNewSurfaceChargeDensity_Int(m,j,perm) - getOldSurfaceChargeDensity_Int(m,j,perm)) / dt_;
    }

    // Calculate gamma_snake values and use as initial values
    void init(int m, const T v, const VCON& conCytosol, const VCON& conExtra)
    {
      // Save membrane potential for each element
      memb_pot_old[m] = v;
      memb_pot_new[m] = v;

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
    void timeStep(int m, const T dt_, const T v, const VCON& conCytosol, const VCON& conExtra)
    {
      // Save membrane potential for each element
      memb_pot_new[m] = v;

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

      // Update gamma values
      for(int j=0; j<conCytosol.size(); j++)
      {
        gamma_ext_new[m][j] = (tau * gamma_ext_old[m][j] + dt_ * gamma_ext_snake[j]) / (tau + dt_);
        gamma_int_new[m][j] = (tau * gamma_int_old[m][j] + dt_ * gamma_int_snake[j]) / (tau + dt_);

        debug_jochen << "   " << gamma_ext_old[m][j] << " -> " << gamma_ext_new[m][j] << std::endl;
        debug_jochen << "   " << gamma_int_old[m][j] << " -> " << gamma_int_new[m][j] << std::endl;
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

    void acceptState()
    {
      memb_pot_old = memb_pot_new;
      gamma_ext_old = gamma_ext_new;
      gamma_int_old = gamma_int_new;
    }

    bool checkIndex(const int i) const
    {
      if(membraneIndices.find(i) == membraneIndices.end())
      {
        DUNE_THROW(Dune::Exception, "Given element #" << i << " is not a membrane element!");
      }
      return true;
    }


    DGF_CON& dgfCon;
    DGF_POT& dgfPot;
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
    std::vector<typename Traits::RangeType> cachedMoriFlux_ext;
    std::vector<typename Traits::RangeType> cachedMoriFlux_int;

    const T dMemb;
    const T temp;
    const T stdCon;
    const T epsElec;
    T time;
    T dt;
    const T tau;

    VPOT memb_pot_old;
    VPOT memb_pot_new;
    V gamma_ext_old;
    V gamma_ext_new;
    V gamma_int_old;
    V gamma_int_new;

    const bool iterateMembraneFlux;

};

#endif /* DUNE_AX1_CHARGE_LAYER_MORI_HH */
