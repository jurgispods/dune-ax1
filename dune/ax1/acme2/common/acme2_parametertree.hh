/*
 * acme2_parametertree.hh
 *
 *  Created on: Jul 15, 2011
 *      Author: jpods
 */

#ifndef ACME2_PARAMETERTREE_HH_
#define ACME2_PARAMETERTREE_HH_

#include <dune/common/parametertree.hh>

class Acme2Parameters : public Dune::ParameterTree
{
  public:

    void init(std::string configFileName_, int argc_, char** argv_)
    {
      general = sub("general");
      membrane = sub("membrane");
      equilibration = sub("equilibration");
      stimulation = sub("stimulation");
      solution_in = sub("solution_in");
      solution_ex = sub("solution_ex");

      configFileName = configFileName_;
      argc = argc_;
      argv = argv_;

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

    double yMax() const
    {
      return general.get("ymax", 500.);
    }

    double xMin() const
    {
      return 0.0;
    }

    double yMin() const
    {
      return 0.0;
    }

    double yMemb() const
    {
      if(useMembrane())
        return general.get("y_memb", 25.);
      else
        return 0.0;
    }

    double dMemb() const
    {
      if(useMembrane())
        return general.get("d_memb", 50.);
      else
        return 0.0;
    }

    double dX() const
    {
      return general.get("dx", 100.);
    }

    double dY() const
    {
      return general.get("dy", 100.);
    }

    double dYCell() const
    {
      return general.get("dy_cell", dY());
    }

    double dYMin() const
    {
      return general.get("dy_min", 0.5);
    }

    double getGammaCon() const
    {
      return general.get("gammaCon", 1.0);
    }

    double getGammaPot() const
    {
      return general.get("gammaPot", 1.0);
    }

    //! Return reduction set in config file if found, default value form configuration header else
    double getReduction() const
    {
      return general.get("newtonReduction",reduction);
    }

    void setReduction(double reduction_)
    {
      reduction = reduction_;
    }

    //! Return abs limit set in config file if found, default value form configuration header else
    double getAbsLimit() const
    {
      return general.get("newtonAbsLimit",absLimit);
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

    bool doContinueSimulation() const
    {
      return general.get("continueSimulation", false);
    }

    std::string getSaveFilename() const
    {
      return general.get("saveFilename", "simulation_states/hackepeter/default_lvl6.dat");
    }

    std::string getLoadFilename() const
    {
      return general.get("loadFilename", "simulation_states/equilibrium/default_lvl6.dat");
    }

    bool doRefineYDirectionGeometrically() const
    {
      return general.get("refineYDirectionGeometrically", false);
    }

    bool doRefineXDirection() const
    {
      return general.get("refineXDirection", false);
    }

    bool doEquilibration() const
    {
      return equilibration.get("doEquilibration", false);
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
      return membrane.get("isActive",true);
    }

    bool doVTKOutput() const
    {
      return general.get("vtkOutput",true);
    }

    bool doFullGnuplotOutput() const
    {
      return general.get("fullGnuplotOutput",true);
    }

    int nMembranes() const
    {
      return general.get("numberMembranes", 1);
    }

    double getCellWidth() const
    {
      return general.get("cellWidth", 10.0e3);
    }

    bool doPrintMatrix() const
    {
      return general.get("printMatrix",false);
    }

    bool doPrintRhs() const
    {
      return general.get("printRhs",false);
    }

    bool doCheckpointing() const
    {
      return general.get("doCheckpointing",true);
    }

    int checkpointInterval() const
    {
      return general.get("checkpointInterval",100);
    }

    std::string getConfigFileName() const
    {
      return configFileName;
    }

    std::string getCommandLineParameters() const
    {
      std::string com = "";
      for(int i=0; i<argc; i++)
      {
        std::string arg = argv[i];
        com += arg + " ";
      }
      return com;
    }

    bool useRowNormPreconditioner() const
    {
      return general.get("useRowNormPreconditioner",true);
    }

    int maxNumberNewtonRestarts() const
    {
      return general.get("maxNumberNewtonRestarts",3);
    }

    bool doExtrapolateXData() const
    {
      return general.get("doExtrapolateXData",false);
    }

    bool closedCell() const
    {
      return general.get("closedCell",true);
    }

    bool createGridFile() const
    {
      return general.get("createGridFile",true);
    }

    std::string getSaveGridFileName() const
    {
      return general.get("saveGridFileName",std::string("acme2.dgf"));
    }

    std::string getLoadGridFileName() const
    {
      return general.get("loadGridFileName",std::string("acme2_EQUILIBRIUM.dgf"));
    }

    bool doReorderMatrix() const
    {
      return general.get("doReorderMatrix",false);
    }

    int vtkSubSamplingLevel() const
    {
      return general.get("vtkSubSamplingLevel",0);
    }


    ParameterTree general;
    ParameterTree equilibration;
    ParameterTree stimulation;
    ParameterTree membrane;
    ParameterTree solution_in;
    ParameterTree solution_ex;

  private:
    std::string configFileName;
    int argc;
    char** argv;

    bool logScaling;
    bool analyticalSolution;
    double reduction;
    double absLimit;
};


#endif /* ACME2_PARAMETERTREE_HH_ */
