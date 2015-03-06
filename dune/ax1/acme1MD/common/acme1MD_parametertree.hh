/*
 * acme1MD_parametertree.hh
 *
 *  Created on: Jul 15, 2011
 *      Author: jpods
 */

#ifndef ACME1MD_PARAMETERTREE_HH_
#define ACME1MD_PARAMETERTREE_HH_

#include <dune/common/parametertree.hh>

class Acme1MDParameters : public Dune::ParameterTree
{
  public:

    void init()
    {
      general = sub("general");
      membrane = sub("membrane");
      equilibration = sub("equilibration");
      stimulation = sub("stimulation");
      solution_in = sub("solution_in");
      solution_ex = sub("solution_ex");

      // Assert tEquilibrium is a multiple of dtEquilibrium
      debug_jochen << "tEquilibrium % dtEquilibrium = " << std::fmod(tEquilibrium(),dtEquilibrium()) << std::endl;
      assert(std::abs(std::fmod(1e3*tEquilibrium(),1e3*dtEquilibrium())) < 1e-6);
    }

    std::string getConfigName()
    {
      return general.get("configName", "default");
    }

    //! \brief get the number of groups
    int nMembraneGroups()
    {
      return membrane.getSubKeys().size();
    }

    bool useLogScaling() const
    {
      return logScaling;
    }

    bool setUseLogScaling(bool logScaling_)
    {
      logScaling = logScaling_;
    }

    bool hasAnalyticalSolution() const
    {
      return analyticalSolution;
    }

    void setHasAnalyticalSolution(bool analyticalSolution_)
    {
      analyticalSolution = analyticalSolution_;
    }

    bool useMembrane() const
    {
      return general.get("useMembrane", true);
    }

    bool refineMembrane() const
    {
      return general.get("refineMembrane", true);
    }

    bool mV_output() const
    {
    	return general.get("mV_output", true);
    }

    std::string initConc() const
    {
      return general.get("initConc", "step");
    }

    double xMax() const
    {
      return general.get("xmax", 500.);
    }

    double xMin() const
    {
      return -xMax();
    }

    double dMemb() const
    {
      if(useMembrane())
        return general.get("d_memb", 50.);
      else
        return 0.0;
    }

    double getGammaCon() const
    {
      return general.get("gammaCon", 1.0);
    }

    double getGammaPot() const
    {
      return general.get("gammaPot", 1.0);
    }

    double getReduction() const
    {
      return reduction;
    }

    void setReduction(double reduction_)
    {
      reduction = reduction_;
    }

    double getAbsLimit() const
    {
      return absLimit;
    }

    void setAbsLimit(double absLimit_)
    {
      absLimit = absLimit_;
    }

    std::string getOutputPrefix() const
    {
      return general.get("outputDir", "out-operator-split") + "/";
    }

    std::string getDiagnosticsFilename() const
    {
      return general.get("diagnosticsFile", "diagnostics.dat");
    }

    double getOutputTimeInterval() const
    {
      return general.get("outputTimeInterval", 1.0);
    }

    bool useAdaptiveTimeStep() const
    {
      return general.get("adaptiveTimeStep", true);
    }

    bool doSaveState() const
    {
      return general.get("saveState", false);
    }

    bool doLoadState() const
    {
      return general.get("loadState", false);
    }

    std::string getSaveFilename() const
    {
      return general.get("saveFilename", "simulation_states/hackepeter/default_lvl6.dat");
    }

    std::string getLoadFilename() const
    {
      return general.get("loadFilename", "simulation_states/equilibrium/default_lvl6.dat");
    }

    bool isLogarithmicGrid() const
    {
      return general.get("logarithmicGrid", false);
    }

    bool doLocalRefine() const
    {
      return general.get("localRefine", false);
    }

    double tEquilibrium() const
    {
      return equilibration.get("tEquilibrium", 0.0);
    }

    double dtEquilibrium() const
    {
      return equilibration.get("dtEquilibrium", 0.1);
    }

    bool doStimulation() const
    {
      return stimulation.get("stimulation",true);
    }

    double tInj_start() const
    {
      return stimulation.get("t_inj_start",1e6);
    }

    double tInj_end() const
    {
      return stimulation.get("t_inj_end",1.5e6);
    }

    bool isMembraneActive() const
    {
      return true;
    }

    ParameterTree general;
    ParameterTree equilibration;
    ParameterTree stimulation;
    ParameterTree membrane;
    ParameterTree solution_in;
    ParameterTree solution_ex;

  private:
    bool logScaling;
    bool analyticalSolution;
    double reduction;
    double absLimit;
};


#endif /* ACME1MD_PARAMETERTREE_HH_ */
