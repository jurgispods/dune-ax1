/*
 * acme2_cyl_parametertree.hh
 *
 *  Created on: Jul 15, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_PARAMETERTREE_HH
#define DUNE_AX1_ACME2CYL_PARAMETERTREE_HH

#include <dune/common/parametertree.hh>

#include <dune/ax1/common/ax1_lfs_tools.hh>

class Acme2CylParameters : public Dune::ParameterTree
{
  public:

    typedef std::tuple<std::string, double, double> MembraneGroupTuple;
    typedef std::vector<MembraneGroupTuple> MembraneGroups;

    void init(std::string configFileName_, int argc_, char** argv_)
    {
      general = sub("general");
      boundary = sub("boundary");
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
      assert(tEquilibrium() >= dtEquilibrium());
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

    void setUseLogScaling(bool logScaling_)
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

    bool useMori() const
    {
      return general.get("useMoriChargeLayerContribution", false);
    }

    bool useMoriCorrection() const
    {
      return useMori() && general.get("useMoriCorrection",false);
    }

    bool useMembrane() const
    {
      return general.get("useMembrane", true);
    }

    bool refineMembrane() const
    {
      return general.get("refineMembrane", false);
    }

    int nMembraneElements() const
    {
      return general.get("nMembraneElements", 1);
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
      return general.get("ymin", 0.0);
    }

    std::vector<double> yMemb() const
    {
      // Return vector of size nMembranes() filled with zeros if no membrane is present
      std::vector<double> yMembDefault(nMembranes(), 0.0);
      if(useMembrane())
        return general.get("y_memb", yMembDefault);
      else
        return yMembDefault;
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

    double dYCellMin() const
    {
      return general.get("dy_cell_min", dYMin());
    }

    double dYMin() const
    {
      return general.get("dy_min", 0.5);
    }

    const std::vector<double>& X() const
    {
      return x;
    }

    void setX(const std::vector<double>& x_)
    {
      x = x_;
      numberElements = (x.size()-1)*(y.size()-1);
      numberNodes = x.size()*y.size();
    }

    const std::vector<double>& Y() const
    {
      return y;
    }

    void setY(const std::vector<double>& y_)
    {
      y = y_;
      numberElements = (x.size()-1)*(y.size()-1);
      numberNodes = x.size()*y.size();
    }

    int nElements() const
    {
      return numberElements;
    }

    int nNodes() const
    {
      return numberNodes;
    }

    int nElementsThisProcess() const
    {
      return numberElementsThisProcess;
    }

    void setNElementsThisProcess(int nElemsProcess)
    {
      numberElementsThisProcess = nElemsProcess;
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

    std::string getTimeStepFile() const
    {
      return general.get("timeStepFile", std::string(""));
    }

    //! Max timestep in [s]
    double getMaxTimeStep() const
    {
      return general.get("maxTimeStep", 0.05e-3);
    }

    //! Min timestep in [s]
    double getMinTimeStep() const
    {
      return general.get("minTimeStep", 0.05e-6);
    }

    //! Max timestep during AP in [s]
    double getMaxTimeStepAP() const
    {
      return general.get("maxTimeStepAP", 0.01e-3);
    }

    int getTimeStepLowerLimitNewtonIt() const
    {
      return general.get("timeStepLowerLimitNewtonIt", 3);
    }

    int getTimeStepUpperLimitNewtonIt() const
    {
      return general.get("timeStepUpperLimitNewtonIt", 5);
    }

    double getPotentialResidualScalingFactor() const
    {
      return general.get("potentialResidualScalingFactor", 1.0);
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

    bool doRefineYDirection() const
    {
      return general.get("refineYDirection", false);
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

    bool doForceEquilibration() const
    {
      return equilibration.get("forceEquilibration", false);
    }

    bool doVolumeScaling() const
    {
      return general.get("doVolumeScaling", false);
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

    std::vector<bool> isMembraneActive() const
    {
      // Return vector of size nMembranes() filled with zeros if no membrane is present
      std::vector<bool> active(nMembranes(), true);
      if(useMembrane())
      {
        active = membrane.get("isActive", active);
        if(active.size() < nMembranes())
          DUNE_THROW(Dune::Exception, "[membrane] key 'isActive' must have " << nMembranes() << " entries!");
      }
      return active;
    }

    bool doVTKOutput() const
    {
      return general.get("vtkOutput",true);
    }

    bool doFullGnuplotOutput() const
    {
      return general.get("fullGnuplotOutput",true);
    }

    bool doHDF5Output() const
    {
#if HAVE_HDF5
      return general.get("hdf5Output",false);
#else
      return false;
#endif
    }

    bool doFullHDF5Output() const
    {
#if HAVE_HDF5
      return general.get("fullHdf5Output",false);
#else
      return false;
#endif
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

    bool doPrintResidual() const
    {
      return general.get("printResidual",false);
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
      return general.get("saveGridFileName",std::string("acme2_cyl.dgf"));
    }

    std::string getLoadGridFileName() const
    {
      return general.get("loadGridFileName",std::string("acme2_cyl_EQUILIBRIUM.dgf"));
    }

    bool doReorderMatrix() const
    {
      return general.get("doReorderMatrix",false);
    }

    int getYOffsetMembrane() const
    {
      return yOffsetMembrane;
    }

    void setYOffsetMembrane(int yOffsetMembrane_)
    {
      yOffsetMembrane = yOffsetMembrane_;
    }

    int vtkSubSamplingLevel() const
    {
      return general.get("vtkSubSamplingLevel",0);
    }

    bool doLoadChannelStates() const
    {
      return general.get("loadChannelStates", false);
    }

    // ======================= BOUNDARIES ============================================
    bool isBoundaryDirichlet_Concentration(std::string boundaryLoc) const
    {
      std::string bnd = boundary.sub("concentration").get(boundaryLoc,"Dirichlet");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isTopBoundaryDirichlet_Concentration() const
    {
      std::string bnd = boundary.sub("concentration").get("top","Dirichlet");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isBottomBoundaryDirichlet_Concentration() const
    {
      std::string bnd = boundary.sub("concentration").get("bottom","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isLeftCytosolBoundaryDirichlet_Concentration() const
    {
      std::string bnd = boundary.sub("concentration").get("leftCytosol","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isRightCytosolBoundaryDirichlet_Concentration() const
    {
      std::string bnd = boundary.sub("concentration").get("rightCytosol","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isLeftExtracellularBoundaryDirichlet_Concentration() const
    {
      std::string bnd = boundary.sub("concentration").get("leftExtracellular","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isRightExtracellularBoundaryDirichlet_Concentration() const
    {
      std::string bnd = boundary.sub("concentration").get("rightExtracellular","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isLeftDebyeBoundaryDirichlet_Concentration() const
    {
      std::string bnd = boundary.sub("concentration").get("leftDebye","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isRightDebyeBoundaryDirichlet_Concentration() const
    {
      std::string bnd = boundary.sub("concentration").get("rightDebye","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isBoundaryDirichlet_Potential(std::string boundaryLoc) const
    {
      std::string bnd = boundary.sub("potential").get(boundaryLoc,"Dirichlet");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isTopBoundaryDirichlet_Potential() const
    {
      std::string bnd = boundary.sub("potential").get("top","Dirichlet");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isBottomBoundaryDirichlet_Potential() const
    {
      std::string bnd = boundary.sub("potential").get("bottom","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isLeftCytosolBoundaryDirichlet_Potential() const
    {
      std::string bnd = boundary.sub("potential").get("leftCytosol","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isRightCytosolBoundaryDirichlet_Potential() const
    {
      std::string bnd = boundary.sub("potential").get("rightCytosol","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isLeftExtracellularBoundaryDirichlet_Potential() const
    {
      std::string bnd = boundary.sub("potential").get("leftExtracellular","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isRightExtracellularBoundaryDirichlet_Potential() const
    {
      std::string bnd = boundary.sub("potential").get("rightExtracellular","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isLeftDebyeBoundaryDirichlet_Potential() const
    {
      std::string bnd = boundary.sub("potential").get("leftDebye","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isRightDebyeBoundaryDirichlet_Potential() const
    {
      std::string bnd = boundary.sub("potential").get("rightDebye","Neumann");
      if(bnd == "Dirichlet") return true;
      else return false;
    }

    bool isMembraneBoundaryNeumann_Potential() const
    {
      std::string bnd = boundary.sub("potential").get("membrane","None");
      if(bnd == "Neumann") return true;
      else return false;
    }

    bool doOutputRootNodeOnly() const
    {
      return general.get("outputRootNodeOnly", true);
    }

    double xNode() const
    {
      return general.get("node_start", 0.0);
    }

    double dNode() const
    {
      return general.get("node_width", xMax());
    }

    void setMembraneGroups(const MembraneGroups& membGroups_)
    {
      membGroups = membGroups_;
    }

    const MembraneGroups& getMembraneGroups() const
    {
      return membGroups;
    }

    ParameterTree& elementGroup(const int groupIndex)
    {
      if(groupIndex == 0)
        return solution_in;
      else if(groupIndex == 1)
        return solution_ex;
      else
      {
        std::vector<std::string> mGroups = membrane.getSubKeys();
        assert(groupIndex-2 < mGroups.size());
        return membrane.sub(mGroups[groupIndex-2]);
      }
    }

    // =================================================================================

    ParameterTree general;
    ParameterTree boundary;
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
    bool timeDependentBoundaryValues;
    double reduction;
    double absLimit;

    std::vector<double> x;
    std::vector<double> y;

    int numberElements;
    int numberNodes;
    int numberElementsThisProcess;
    int numberNodesThisProcess;
    int yOffsetMembrane;

    MembraneGroups membGroups;
};


#endif /* DUNE_AX1_ACME2CYL_PARAMETERTREE_HH */
