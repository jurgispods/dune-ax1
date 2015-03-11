/*
 * acme1_parametertree.hh
 *
 *  Created on: Jul 15, 2011
 *      Author: jpods
 */

#ifndef ACME1_PARAMETERTREE_HH_
#define ACME1_PARAMETERTREE_HH_

#include <dune/common/parametertree.hh>

class Acme1Parameters : public Dune::ParameterTree
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

    ParameterTree general;
    ParameterTree membrane;
    ParameterTree equilibration;
    ParameterTree stimulation;
    ParameterTree solution_in;
    ParameterTree solution_ex;

  private:
    bool logScaling;
    bool analyticalSolution;
    double reduction;
    double absLimit;
};


#endif /* ACME1_PARAMETERTREE_HH_ */
