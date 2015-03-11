/*
 * acme2_output.hh
 *
 *  Created on: Aug 11, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2_OUTPUT_HH
#define DUNE_AX1_ACME2_OUTPUT_HH

#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/pdelab/common/clock.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/common/functionutilities.hh>
#include <dune/pdelab/function/selectcomponent.hh>

#include <dune/ax1/common/ax1_output2d.hh>
#include <dune/ax1/common/chargedensitygridfunction.hh>
#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/error_norms.hh>
#include <dune/ax1/common/ionfluxgridfunction.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/ax1_boundaryfunction_membranefunction_adapter.hh>
#include <dune/ax1/common/ax1_discretegridfunction.hh>
#include <dune/ax1/common/ax1_gridfunctionadapters.hh>
#include <dune/ax1/common/ax1_subdomainindexgridfunction.hh>
#include <dune/ax1/common/ax1_membrane_validation_gridfunctions.hh>
#include <dune/ax1/common/ax1_multidomaingridfunction.hh>
#include <dune/ax1/common/ax1_multidomaingridfunctionextension.hh>
#include <dune/ax1/common/ax1_multidomaingridfunctionrestriction.hh>
#include <dune/ax1/common/channelgridfunction.hh>
#include <dune/ax1/common/membranefluxgridfunction.hh>
#include <dune/ax1/common/membranefluxgridfunctionhelper.hh>
#include <dune/ax1/common/membranepotentialgridfunction.hh>

#include <dune/ax1/acme2/common/acme2_solutionvectors.hh>



//template<class Acme2Traits, class GV, class SubGV, class GFS_CON, class GFS_POT, class U_CON, class U_POT, class PHYSICS>
template<class Acme2Traits>
class Acme2Output
{
  public:

    // Boilerplate code to ease type usage
    typedef typename Acme2Traits::Coord Coord;
    typedef typename Acme2Traits::Real Real;
    typedef typename Acme2Traits::GridView GV;
    typedef typename Acme2Traits::SubGridView SubGV;
    typedef typename Acme2Traits::MultiGFS MultiGFS;
    typedef typename Acme2Traits::GFS_CON_SUB GFS_CON; // SubGFS!
    typedef typename Acme2Traits::GFS_POT_SUB GFS_POT; // SubGFS!
    typedef typename Acme2Traits::U U;
    //typedef typename Acme2Traits::U_CON U_CON;
    //typedef typename Acme2Traits::U_POT U_POT;
    typedef typename Acme2Traits::Physics PHYSICS;

    // Not used at the moment, but might be useful instead of using SelectComponentGridFunctionAdapter
    template<int SPECIES>
    struct GFS_ION
    {
      typedef Dune::PDELab::GridFunctionSubSpace<GFS_CON,SPECIES> Type;
    };

    static const int SUBSAMPLING_POINTS = 0;

    const std::string vtkOutputBase_Domain;
    const std::string vtkOutputBase_Membrane;

    // A bunch of typedef's for all the gridfunctions
    struct Acme2OutputTraits
    {
      // DGFs
      typedef Dune::PDELab::DiscreteGridFunction        <GFS_POT,U> DGF_POT;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_POT,U> DGF_POT_GRAD;
      typedef Ax1MultiDomainGridFunctionRestriction<SubGV,DGF_POT_GRAD,PHYSICS> DGF_POT_GRAD_ELEC;

      typedef Dune::PDELab::VectorDiscreteGridFunction  <GFS_CON,U> DGF_CON;
      typedef Dune::PDELab::VectorDiscreteGridFunctionGradient<GFS_CON,U> DGF_CON_GRAD;

      typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_CON, PHYSICS> DGF_CON_MD;
      typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_CON_GRAD, PHYSICS> DGF_CON_GRAD_MD;

      typedef ChargeDensityGridFunction<DGF_CON_MD,PHYSICS> GF_CD;

      // TODO Do not restrict and interpolate again, can this be simplified?
      typedef Dune::PDELab::SelectComponentGridFunctionAdapter<DGF_CON_MD,1> DGF_SINGLE_CON_MD;
      typedef Dune::PDELab::SelectComponentGridFunctionAdapter<DGF_CON_GRAD_MD,GV::dimension> DGF_SINGLE_CON_GRAD_MD;

      // Meta grid function for the electrolyte ion flux taking only gridfunctions defined on the elec subdomain
      typedef IonFluxGridFunction<DGF_CON, DGF_CON_GRAD, DGF_POT_GRAD_ELEC, PHYSICS> GF_ELEC_FLUX;
      typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_ELEC_FLUX,2> GF_SINGLE_ELEC_FLUX;

      typedef MembraneFluxGridFunction<DGF_CON_MD,DGF_POT,PHYSICS> BGF_MEMB_FLUX;
      typedef MembraneFluxGridFunction_DiffusionTerm<BGF_MEMB_FLUX> BGF_MEMB_FLUX_DIFF;
      typedef MembraneFluxGridFunction_DriftTerm<BGF_MEMB_FLUX> BGF_MEMB_FLUX_DRIFT;

      // Grid functions living on the membrane (only one component for the flux in normal direction)
      typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_FLUX,PHYSICS> GF_MEMB_FLUX;
      typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_FLUX_DIFF,PHYSICS> GF_MEMB_FLUX_DIFF;
      typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_FLUX_DRIFT,PHYSICS> GF_MEMB_FLUX_DRIFT;
      typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_MEMB_FLUX,1> GF_SINGLE_MEMB_FLUX;
      typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_MEMB_FLUX_DIFF,1> GF_SINGLE_MEMB_FLUX_DIFF;
      typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_MEMB_FLUX_DRIFT,1> GF_SINGLE_MEMB_FLUX_DRIFT;

      // Grid functions living on the whole domain (one component per coordinate axis)
      typedef Ax1MultiDomainGridFunction<GV,PHYSICS,GF_ELEC_FLUX,GF_MEMB_FLUX> GF_ION_FLUX;
      typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_ION_FLUX,GV::dimension> GF_SINGLE_ION_FLUX;

      typedef MembranePotentialGridFunction<DGF_POT, PHYSICS> GF_MEMB_POT_MD;
      typedef Ax1MultiDomainGridFunctionRestriction<SubGV,GF_MEMB_POT_MD,PHYSICS> GF_MEMB_POT;

      typedef ChannelGridFunction<SubGV,Real,PHYSICS> GF_CHANNEL;

      // Analytical solutions
      typedef typename PHYSICS::Traits::SOLUTION_CON ANALYTICAL_SOLUTION_CON;
      typedef typename PHYSICS::Traits::SOLUTION_POT ANALYTICAL_SOLUTION_POT;
      typedef Dune::PDELab::SelectComponentGridFunctionAdapter<ANALYTICAL_SOLUTION_CON,1> SINGLE_ANALYTICAL_SOLUTION_CON;
    };


    typedef Acme2OutputTraits Traits;
    typedef typename Traits::DGF_POT::Traits::DomainType DomainType;

    struct DiagnosticInfo
    {
      DiagnosticInfo()
      : time(0.0),
        iterations(0), dt(0.0),
        l2ErrorPot(0.0), maxErrorPot(0.0),
        l2ErrorCon(NUMBER_OF_SPECIES), maxErrorCon(NUMBER_OF_SPECIES),
        tEquilibrium(0.0),
        maxDiffCon(0),
        maxDiffPot(0)
      {}

      // diagnostic variables for one time step
      double time;
      int iterations;
      double dt;
      Real l2ErrorPot;
      Real maxErrorPot;
      std::vector<Real> l2ErrorCon;
      std::vector<Real> maxErrorCon;
      double tEquilibrium;

      std::vector<Real> maxDiffCon;
      std::vector<Real> maxDiffPot;

      std::map<std::string, Real> debugData;


      Real getl2ErrorCon(int i) const
      {
        return l2ErrorCon[i];
      }

      void setl2ErrorCon(int i, Real err)
      {
        l2ErrorCon[i] = err;
      }

      Real getMaxErrorCon(int i) const
      {
        return maxErrorCon[i];
      }

      void setMaxErrorCon(int i, Real err)
      {
        maxErrorCon[i] = err;
      }

      void setMaxDiffConPot(Real maxDiffCon_, Real maxDiffPot_)
      {
        maxDiffCon.push_back(maxDiffCon_);
        maxDiffPot.push_back(maxDiffPot_);
      }


      std::string getHeader()
      {
        std::stringstream header;
        int count = 6;
        header << "#  (1)time        (2)dt             (3)#iterations   (4)[pot] L2 error (5)[pot] max error";
        for(int i=0; i<NUMBER_OF_SPECIES; ++i)
        {
          header << " (" << (count++) << ")[" << ION_NAMES[i] << "] L2 error";
          header << " (" << (count++) << ")[" << ION_NAMES[i] << "] max error";
        }
        for(typename std::map<std::string,Real>::const_iterator it = debugData.begin();
              it != debugData.end(); ++it)
        {
          header << "   (" << (count++) << ")" << it->first;
        }
        return header.str();
      }

      void registerDebugData(std::string key, Real value = -1e100)
      {
        debugData[key] = value;
      }

      void setDebugData(std::string key, Real value)
      {
        if(debugData.count(key))

          debugData[key] = value;
        else
          DUNE_THROW(Dune::Exception, "Key '"  << key << "' is not present in DiagnosticInfo::debugData map!");
      }

      std::vector<Real> getVectorRepresentation() const
      {
        int count = 0;

        std::vector<Real> vec;
        //vec[0] = time;
        vec.push_back(dt);
        vec.push_back((Real) iterations);
        vec.push_back(l2ErrorPot);
        vec.push_back(maxErrorPot);
        for(int i=0; i<l2ErrorCon.size()*2; i=i+2)
        {
          vec.push_back(l2ErrorCon[i]);
          vec.push_back(maxErrorCon[i]);
        }
        for(typename std::map<std::string,Real>::const_iterator it = debugData.begin();
            it != debugData.end(); ++it)
        {
          vec.push_back(it->second);
        }
        return vec;
      }

      void clear()
      {
        iterations = 0;
        dt = 0.0;
        l2ErrorPot = 0.0;
        maxErrorPot = 0.0;
        l2ErrorCon.clear(); l2ErrorCon.resize(NUMBER_OF_SPECIES);
        maxErrorCon.clear(); maxErrorCon.resize(NUMBER_OF_SPECIES);

        maxDiffCon.clear(); maxDiffCon.resize(0);
        maxDiffPot.clear(); maxDiffPot.resize(0);
      }
    };

    Acme2Output(const GV& gv_, const SubGV& elecGV_, const SubGV& membGV_,
        MultiGFS& multigfs_, GFS_CON& gfsCon_, GFS_POT& gfsPot_,
        Acme2SolutionVectors<Acme2Traits>& solutionVectors_,
        //U_CON& uCon_, U_POT& uPot_,
        typename Traits::BGF_MEMB_FLUX& bgfMembraneFlux_,
        PHYSICS& physics_,
        int intorderPot_=2, int intorderCon_=2) :
      vtkOutputBase_Domain("acme2"),
      vtkOutputBase_Membrane("acme2_membrane"),
      gv(gv_),
      elecGV(elecGV_),
      membGV(membGV_),
      multigfs(multigfs_),
      gfsCon(gfsCon_),
      gfsPot(gfsPot_),
      solutionVectors(solutionVectors_),
      u(solutionVectors_.unew),
      //uCon(solutionVectors_.unewCon),
      //uPot(solutionVectors_.unewPot),
      bgfMembraneFlux(bgfMembraneFlux_), // Reference!
      physics(physics_),
      intorderPot(intorderPot_),
      intorderCon(intorderCon_),
      con(NUMBER_OF_SPECIES),
      conGrad(NUMBER_OF_SPECIES),
      membraneFlux(NUMBER_OF_SPECIES),
      membraneFlux_DiffTerm(NUMBER_OF_SPECIES),
      membraneFlux_DriftTerm(NUMBER_OF_SPECIES),
      ionFlux(NUMBER_OF_SPECIES),
      solutionCon(NUMBER_OF_SPECIES),
      totalInitParticles(NUMBER_OF_SPECIES),
      totalParticles(NUMBER_OF_SPECIES),
      sumIonFluxIntegrals(0.0),
      vtkwriterGV(gv,physics_.getParams().vtkSubSamplingLevel()),
      vtkwriterSubGV(membGV,0), // Don't do subsampling on the membrane!
      vtkBasenameDomain(physics.getParams().getOutputPrefix() + "paraview/" + vtkOutputBase_Domain),
      vtkBasenameMembrane(physics.getParams().getOutputPrefix() + "paraview/" + vtkOutputBase_Membrane),
      checkpointBasename(physics.getParams().getOutputPrefix() + "checkpoints/checkpoint"),
      fn_vtk_GV(vtkBasenameDomain),
      fn_vtk_SubGV(vtkBasenameMembrane),
      fn_checkpoints(checkpointBasename),
      filePrefix(physics.getParams().getOutputPrefix()),
      dgfConElec(gfsCon,u),
      dgfConGradElec(gfsCon,u),
      dgfCon(gv,dgfConElec,physics),
      dgfConGrad(gv,dgfConGradElec,physics),
      dgfPot(gfsPot,u),
      dgfPotGrad(gfsPot,u),
      dgfPotGradElec(elecGV,dgfPotGrad,physics),
      gfChargeDensity(dgfCon,dgfCon,physics),
      gfElecFlux(dgfConElec,dgfConGradElec,dgfPotGradElec,physics),
      gfMembranePotentialMD(dgfPot, physics),
      gfMembranePotential(membGV, gfMembranePotentialMD, physics),
      gfMembraneFlux(membGV, bgfMembraneFlux, physics),
      bgfMembraneFlux_DiffTerm(bgfMembraneFlux),
      gfMembraneFlux_DiffTerm(membGV, bgfMembraneFlux_DiffTerm, physics),
      bgfMembraneFlux_DriftTerm(bgfMembraneFlux),
      gfMembraneFlux_DriftTerm(membGV, bgfMembraneFlux_DriftTerm, physics),
      gfIonFlux(gv,gfElecFlux,gfMembraneFlux,physics),
      gfAnalyticalSolutionCon(gv,physics.getParams()),
      gfAnalyticalSolutionPot(gv,physics.getParams()),
      vtkOutput(physics.getParams().doVTKOutput()),
      doFullGnuplotOutput(physics.getParams().doFullGnuplotOutput()),
      timer(false)
    {
      timer.start();

      std::vector<std::string> outputDirs = {"paraview", "checkpoints", "config"};
      int incrementVTK = 0;
      int incrementCheckpoint = 0;

      // Check if output directory exists and create subdirs
      if(access(physics.getParams().getOutputPrefix().c_str(),0) != 0)
      {
        DUNE_THROW(Dune::IOError, "The specified output directory '" +
            physics.getParams().getOutputPrefix() + "' does not exist. Please create it first!");
      } else {

        // Create subdirs
        for(int i=0; i<outputDirs.size(); i++)
        {
          // Create single subdir and check for success
          if(access((physics.getParams().getOutputPrefix() + outputDirs[i]).c_str(),0) != 0)
          {
            int status = mkdir((physics.getParams().getOutputPrefix() + outputDirs[i]).c_str(),
              S_IRWXU | S_IRWXG | S_IROTH); // r+w+x access for user and group, r for other

            if(status != 0)
              DUNE_THROW(Dune::IOError,
                  "Could not create directory '" + (physics.getParams().getOutputPrefix() + outputDirs[i])
                  + "'");
          }
        }
      }

      // Clean up output dir only if the data shall not be appended to existing data files
      if(not (physics.getParams().doLoadState() && physics.getParams().doContinueSimulation()))
      {
        // Remove old VTK files
        std::string command = "rm " + vtkBasenameDomain + "*";
        debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
        int unusedVariable = std::system(command.c_str());

        // Remove old checkpoint files
        command = "rm " + checkpointBasename + "*";
        debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
        unusedVariable = std::system(command.c_str());

        // Remove old files from config directory
        command = "rm " + physics.getParams().getOutputPrefix() + "config/*";
        debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
        unusedVariable = std::system(command.c_str());

        // Following: Initialize (i.e., create and write headers) all gnuplot output files
        Output::gnuplotInitialize(physics.getParams().getDiagnosticsFilename(),
            physics.getParams().getOutputPrefix(), diagInfo.getHeader());

        Output::gnuplotInitialize("debug_output.dat",
            physics.getParams().getOutputPrefix(), diagInfo.getHeader());

        Output::gnuplotInitialize("operator_split_debug.dat",
            physics.getParams().getOutputPrefix());

        // gnuplot file initialization
        Output::gnuplotInitialize( "cd.dat", filePrefix );
        Output::gnuplotInitialize( "pot.dat", filePrefix);
        Output::gnuplotInitialize( "pot_grad.dat", filePrefix );
        Output::gnuplotInitialize( "memb_pot.dat", filePrefix );

        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          // gnuplot
          Output::gnuplotInitialize( physics.getIonName(j)+".dat", filePrefix);
          Output::gnuplotInitialize( physics.getIonName(j)+"_grad.dat", filePrefix);
          Output::gnuplotInitialize( "memb_flux_" + physics.getIonName(j)+".dat", filePrefix);
          Output::gnuplotInitialize( "memb_flux_diff_term_" + physics.getIonName(j)+".dat", filePrefix);
          Output::gnuplotInitialize( "memb_flux_drift_term_" + physics.getIonName(j)+".dat", filePrefix);
          Output::gnuplotInitialize( "flux_" + physics.getIonName(j)+".dat", filePrefix);
          Output::gnuplotInitialize( "con_error_" + physics.getIonName(j) + ".dat", filePrefix );
        }

        for(int k=0; k<physics.getMembrane().getChannelSet().size(); k++)
        {
          std::stringstream chStr;
          chStr << std::setw(2) << std::setfill('0') << k;
          Output::gnuplotInitialize("channel_" + physics.getMembrane().getChannelSet().getChannel(k).getName() + ".dat",
              filePrefix);
        }
      } else {
        // Continue simulation => keep old output files

        DIR *dir;
        struct dirent *ent;

        // Determine the number of written output VTK files
        dir = opendir ((physics.getParams().getOutputPrefix() + "paraview").c_str());
        if (dir != NULL) {

          /* print all the files and directories within directory */
          while ((ent = readdir(dir)) != NULL)
          {
            std::string filename(ent->d_name);
            if(filename.find(vtkOutputBase_Membrane) != std::string::npos)
            {
              incrementVTK++;
            }
          }
          closedir (dir);
        } else {
          DUNE_THROW(Dune::IOError,
            "Error reading directory '" + (physics.getParams().getOutputPrefix() + "paraview")
            + "'");
        }

        // Determine the number of written output checkpoint files
        dir = opendir ((physics.getParams().getOutputPrefix() + "checkpoints").c_str());
        if (dir != NULL) {

          /* print all the files and directories within directory */
          while ((ent = readdir(dir)) != NULL)
          {
            std::string filename(ent->d_name);
            if(filename.find("checkpoint") != std::string::npos)
            {
              incrementCheckpoint++;
            }
          }
          closedir (dir);
        } else {
          DUNE_THROW(Dune::IOError,
            "Error reading directory '" + (physics.getParams().getOutputPrefix() + "paraview")
            + "'");
        }

        debug_verb << "Found " << incrementVTK << " old VTK files, incrementing counter accordingly."
                   << std::endl;
        for(int i=0; i<incrementVTK; i++)
        {
          fn_vtk_GV.increment();
          fn_vtk_SubGV.increment();
        }
        debug_verb << "Found " << incrementCheckpoint << " old checkpoint files, incrementing counter accordingly."
                   << std::endl;
        for(int i=0; i<incrementCheckpoint; i++)
        {
          fn_checkpoints.increment();
        }
      }

      // Copy config file to config dir
      std::string command = "cp --backup=t " + physics.getParams().getConfigFileName() + " "
          + physics.getParams().getOutputPrefix() + "config";
      debug_jochen << "Copying config file using command '" << command << "'" << std::endl;
      int unusedVariable = std::system(command.c_str());

      // Create file containing the command line parameters used to start the program
      command = "echo \"" + physics.getParams().getCommandLineParameters()
          + "\" > command_line.dat && date >> command_line.dat "
          + "&& mv --backup=t command_line.dat "
          + physics.getParams().getOutputPrefix() + "config";
      debug_jochen << "Saving command line parameters using command '" << command << "'" << std::endl;
      unusedVariable = std::system(command.c_str());

      // Generate git version control info
      command = "bash git-info.bash && mv --backup=t git-info.dat "
          + physics.getParams().getOutputPrefix() + "config";
      debug_jochen << "Saving git info using command '" << command << "'" << std::endl;
      unusedVariable = std::system(command.c_str());


      // Create a DGF file from the multidomain GV
      if(physics.getParams().createGridFile() || physics.getParams().doEquilibration())
      {
        Dune::DGFWriter<GV> dgfWriter(gv);
        std::string filename = physics.getParams().getOutputPrefix()
            + "config/" + physics.getParams().getSaveGridFileName();
        dgfWriter.write(filename);
      }



      timer.stop();
    }

    void writeStep(double time_)
    {
      double start = timer.elapsed();
      timer.start();

      time = time_;
      infoStream.str("");
      infoStream << "time: " << time;
      diagInfo.time = time;

      typename Traits::GF_CD::Traits::RangeType cdIntegral;
      Dune::PDELab::integrateGridFunction(gfChargeDensity, cdIntegral, intorderCon);
      debug_verb << "^^^ TOTAL charge = " << cdIntegral << std::endl;

      Tools::integrateGridFunctionOverSubdomain(gfChargeDensity, physics, CYTOSOL, cdIntegral, intorderCon);
      debug_verb << "^^^^ [" << SUBDOMAIN_NAMES[CYTOSOL] << "] charge = " << cdIntegral << std::endl;
      Tools::integrateGridFunctionOverSubdomain(gfChargeDensity, physics, ES, cdIntegral, intorderCon);
      debug_verb << "^^^^ [" << SUBDOMAIN_NAMES[ES] << "] charge = " << cdIntegral << std::endl;
      Tools::integrateGridFunctionOverSubdomain(gfChargeDensity, physics, MEMBRANE, cdIntegral, intorderCon);
      debug_verb << "^^^^ [" << SUBDOMAIN_NAMES[MEMBRANE] << "] charge = " << cdIntegral << std::endl;
      if(std::abs(cdIntegral) > 0)
      {
        debug_warn << "WARNING: Charge density inside membrane is non-zero!" << std::endl;
      }

      double elapsed = timer.elapsed() - start;
      debug_info << "[Acme2Output] Output preparation time: " << elapsed << " s" << std::endl;

      writePotentialOutput();
      elapsed = timer.elapsed() - start;
      debug_info << "[Acme2Output] Potential output time: " << elapsed << " s" << std::endl;
      writeMembraneOutput();
      elapsed = timer.elapsed() - start;
      debug_info << "[Acme2Output] Membrane output time: " << elapsed << " s" << std::endl;
      writeConcentrationOutput();
      elapsed = timer.elapsed() - start;
      debug_info << "[Acme2Output] Concentration output time: " << elapsed << " s" << std::endl;



      // TODO **************************************************************************************
      std::vector<Real> vec = diagInfo.getVectorRepresentation();
      //debug_verb << "***" << vec[0] <<" - " << vec[1] << std::endl;
      GnuplotTools2D::gnuplotAddLine(physics.getParams().getDiagnosticsFilename(), diagInfo.time, vec, filePrefix);

      if(vtkOutput)
      {
        //std::remove(fn_vtk_GV.getName());
        vtkwriterGV.write(fn_vtk_GV.getName(),Dune::VTK::OutputType::ascii);
        vtkwriterGV.clear();
        vtkwriterSubGV.write(fn_vtk_SubGV.getName(),Dune::VTK::OutputType::ascii);
        vtkwriterSubGV.clear();
        fn_vtk_GV.increment();
        fn_vtk_SubGV.increment();
      }

      elapsed = timer.stop() - start;
      debug_info << "[Acme2Output] TOTAL output time: " << elapsed << " s" << std::endl;
    }

    DiagnosticInfo& getDiagInfo()
    {
      return diagInfo;
    }

    void initDiagInfoFile()
    {
      Output::gnuplotInitialize(physics.getParams().getDiagnosticsFilename(),
          physics.getParams().getOutputPrefix(), diagInfo.getHeader());
    }

    void printAllCoeffs()
    {
      Output::printSingleCoefficientVector(u, "u");
    }

    double getTotalOutputTime()
    {
      return timer.elapsed();
    }


    /**
     * Function to write out a complete set of solution vectors which can be used to continue a
     * simulation later or to set a defined initial state for a new simulation
     *
     * @param time
     * @param dt
     * @param filename
     */
    void saveState(double time, double dt, std::string filename = std::string(""))
    {
      timer.start();

      // No filename is given => this is a checkpoint output
      if(filename.size() == 0)
      {
        //TODO Create concurrent checkpointing filename
        filename = fn_checkpoints.getName();
        filename += ".dat";
        fn_checkpoints.increment();
      }

      Ax1SimulationData<Real> simulationData(gv.size(0), time, dt, filename);

      std::vector<Real> stdVec;

      std::vector<typename GV::template Codim<GV::dimension>::Geometry::GlobalCoordinate > nodePositions;
      // Generate vector of node coordinates
      for(typename GV::template Codim<GV::dimension>::Iterator nit
          = gv.template begin<GV::dimension>(); nit != gv.template end<GV::dimension>(); ++nit)
      {
        for(int i=0; i<nit->geometry().corners(); i++)
        {
          nodePositions.push_back(nit->geometry().corner(i));
        }
      }

      // TODO Remove this and save the DGF file properly instead
      // Node positions
      Tools::flattenVector(nodePositions, stdVec);
      simulationData.addVector("nodePositions", stdVec);

      // Fill simulation data with std::vectors containg data from solutionVectors
      solutionVectors.serialize(simulationData);

      Ax1SimulationState<Real> state;
      state.setData(simulationData);
      state.saveState();

      timer.stop();
    }

    void loadState(double& time, double& dt, std::string filename)
    {
      timer.start();

      Ax1SimulationState<Real> state;
      state.loadState(filename);

      Ax1SimulationData<Real>& simulationData = state.getData();

      // Temporary variable; all solution vector will be loaded into this single vector one after another
      std::vector<Real> stdVec;

      bool extrapolateXData = physics.getParams().doExtrapolateXData();

      //debug_jochen << "gv.size(0) = " << gv.size(0) << std::endl;
      //debug_jochen << "simulationData.getNElements() = " << simulationData.getNElements() << std::endl;
      //debug_jochen << "simulationData.getNElements() = " << simulationData.getNElements() << std::endl;
      //debug_jochen << "solutionVectors.uold.flatsize() = " << solutionVectors.uold.flatsize() << std::endl;
      //debug_jochen << "simulationData.getVector(\"uold\").data.size() = "
      //    << simulationData.getVector("uold").data.size() << std::endl;
      //debug_jochen << "solutionVectors.uoldPot.flatsize() = " << solutionVectors.uoldPot.flatsize() << std::endl;
      //debug_jochen << "simulationData.getVector(\"uoldPot\").data.size() = "
      //    << simulationData.getVector("uoldPot").data.size() << std::endl;
      //debug_jochen << "solutionVectors.uoldCon.flatsize() = " << solutionVectors.uoldCon.flatsize() << std::endl;
      //debug_jochen << "simulationData.getVector(\"uoldCon\").data.size() = "
      //    << simulationData.getVector("uoldCon").data.size() << std::endl;

      //int nElementsFactor = gv.size(0) / simulationData.getNElements();
      //debug_jochen << "nElementsFactor = " << nElementsFactor << std::endl;

      if(gv.size(0) != simulationData.getNElements()
          || solutionVectors.uold.flatsize() != simulationData.getVector("uold").data.size())
      {
        // Throw error when number of nodes does not match;
        // When extrapolateXData is desired, disable this check
        if(not extrapolateXData /*||
            (extrapolateXData &&  (gv.size(0) % simulationData.getNElements() != 0))*/)
        {
          DUNE_THROW(Dune::Exception,
            "Error [# nodes] loading saved simulation state from file '" << filename
            << "', was the saved state generated with the same grid and finite element order?");
        } else {
          debug_info << "Number of nodes for the data to be loaded does not match number of nodes in "
              << "the current grid. Will try to interpolate data now!" << std::endl;
        }
      }


      // Use a more accurate check and compare node positions!
      // In the case extrapolateXData is set, only check for y coordinate consistency!
      bool success = true;
      try {
        stdVec = simulationData.getVector("nodePositions").data;
        int count = 0;

        //debug_jochen << "Old node positions: " << std::endl;
        //Output::printVector(stdVec);

        typename GV::template Codim<GV::dimension>::Geometry::GlobalCoordinate posNew(1e100);
        typename GV::template Codim<GV::dimension>::Geometry::GlobalCoordinate posOld(1e100);

        // Iterate over _nodes_!
        for(typename GV::template Codim<GV::dimension>::Iterator nit = gv.template begin<GV::dimension>();
            nit != gv.template end<GV::dimension>(); ++nit)
        {
          if(extrapolateXData)
          {
            // NEW NODES
            // Skip nodes with identical y coordinates
            while(nit != gv.template end<GV::dimension>() && posNew[1] == nit->geometry().corner(0)[1])
            {
              //debug_jochen << "Skipping NEW node " << nit->geometry().corner(0) << std::endl;
              ++nit;
            }
            // Break search loop as soon as new nodes were completely checked
            if(nit == gv.template end<GV::dimension>()) break;

            // OLD NODES
            // Skip nodes with identical y coordinates
            while(count < stdVec.size() && posOld[1] == stdVec[count+1])
            {
              //debug_jochen << "Skipping OLD node #" << count << "-" << (count+1) << " "
              //    << stdVec[count] << " " << stdVec[count+1] << std::endl;
              count += nit->geometry().corner(0).size();
            }
            // Break search loop as soon as old node vector was completely checked
            if(count >= stdVec.size()) break;
          }

          posNew = nit->geometry().corner(0);
          posOld = stdVec[count+1];

          // j: coordinate axis (x=0, y=1)
          // Do _not_ check node x coordinates in case we extrapolate the x data!
          for(int j=(extrapolateXData ? 1 : 0); j<nit->geometry().corner(0).size(); j++)
          {
            //debug_jochen << "CHECK nit->geometry().corner(" << 0 << ")[" << j << "] "
            //  << nit->geometry().corner(0)[j] << " <-> "
            //  << stdVec[count+j] << " stdVec["
            //  << (count+j) << "]" << std::endl;

            if(nit->geometry().corner(0)[j] != stdVec[count+j])
            {
              success = false;
              break;
            }
          }

          // Cancel check as soon as one mismatch was found
          if(! success) break;

          // Increment old node vector index
          count += nit->geometry().corner(0).size();
        }
      } catch(Dune::Exception& e) {
        debug_warn << "No vector containing node positions could be found in loaded simulation data!" << std::endl;
        debug_warn << "===================================================" << std::endl;
        debug_warn << e << std::endl;
        debug_warn << "===================================================" << std::endl;
      }

      if(! success)
      {
        DUNE_THROW(Dune::Exception, "Error [node positions] loading saved simulation state from file '"
            << filename << "', was the saved state generated with the same grid and finite element order?");
      }


      // Following: Hardcoded extrapolation for uoldpot/unewpot and uoldCon/unewCon
      // Think this is ok for our needs.
      /**
       * This is the most general approach to restore data from the old grid: Restore the old grid!
       * Then wrap the loaded solutiuon vectors into grid functions and use the GF adapters from
       * ax1_gridfunctionadapters.hh to interpolate the DOFs onto the new grid. Works nicely!
       *
       * Attention: The gridfunction adapters mentioned above have a very bad performance (they
       * need to loop over the whole gridview for each function call in the worst case), so the
       * following code might be very slow for large grids. This is however nothing compared to
       * the time one would need to perform the full equilibration on the large grid.
       * TODO One might still think about more efficient implementations exploiting information about the
       * grid (Cartesian, structured, ordering of elements in the gridview)
       */
      if(extrapolateXData && (gv.size(0) != simulationData.getNElements()
          || solutionVectors.uold.flatsize() != simulationData.getVector("uold").data.size()) )
      {
        // Load grid used for equilibration from DGF file

        // TODO Get this type from somewhere
        //typedef Dune::UGGrid<2> BaseGrid;
        Dune::GridPtr<BaseGrid> gridptr(physics.getParams().getLoadGridFileName());
        BaseGrid& equilibrationGrid = *gridptr;

        debug_jochen << "Loaded old UG gridview from dgf file '"
            << physics.getParams().getLoadGridFileName() << "'." << std::endl;

        // Wrap into multidomain grid in order to be able to use all types/classes from Acme2Setup
        typename Acme2Traits::GridType md_equilibrationGrid(equilibrationGrid,false);

        //typedef typename Acme2Traits::GridType::LeafGridView EquiGV;
        GV equilibrationGV = md_equilibrationGrid.leafView();

        debug_jochen << "Loaded equilibration gridview has " << equilibrationGV.size(0) << " elements" << std::endl;

        // Define a grid function which allows tagging of elements of the just loaded equilibration gridview
        Ax1SubdomainIndexGridFunction<GV,Real,PHYSICS> subdomainIndexGF(gv, physics);
        Ax1FineGridRestrictionGridFunction<GV,Ax1SubdomainIndexGridFunction<GV,Real,PHYSICS>,
            PHYSICS> subdomainIndexGF_Equi(equilibrationGV,subdomainIndexGF,physics,2); // maxFineElements=2!

        debug_jochen << "Beginning marking equilibration md grid" << std::endl;
        // Now reconstruct the MultiDomainGrid subdoamin structure!
        typedef typename Acme2Traits::GridType::SubDomainGrid SubDomainGrid;
        SubDomainGrid& elecGrid_Equi = md_equilibrationGrid.subDomain(0);
        SubDomainGrid& membGrid_Equi = md_equilibrationGrid.subDomain(1);
        typedef typename SubDomainGrid::LeafGridView SDGV;

        SDGV elecGV_Equi = elecGrid_Equi.leafView();
        SDGV membGV_Equi = membGrid_Equi.leafView();

        md_equilibrationGrid.startSubDomainMarking();
        for (typename GV::template Codim<0>::Iterator eit = equilibrationGV.template begin<0>();
            eit != equilibrationGV.template end<0>(); ++eit)
        {

          typename Ax1FineGridRestrictionGridFunction<GV,Ax1SubdomainIndexGridFunction<GV,Real,PHYSICS>,
                  PHYSICS>::Traits::RangeType myGroupIndex(-5);

          // Evaluate restriction gridfunction at center to get subdomainIndex for the coarse entity
          subdomainIndexGF_Equi.evaluate(*eit,eit->geometry().local(eit->geometry().center()),myGroupIndex);

          //debug_jochen << "Equi Element @ " << eit->geometry().center() << " has group index "
          //              << myGroupIndex << std::endl;

         // Convert double to int
         int subdomainIndex = myGroupIndex[0];

                   // Init subgrids
          switch(subdomainIndex)
          {
            case CYTOSOL:
              md_equilibrationGrid.addToSubDomain(0,*eit);
              break;
            case ES:
              md_equilibrationGrid.addToSubDomain(0,*eit);
              break;
            case MEMBRANE:
              md_equilibrationGrid.addToSubDomain(1,*eit);
              break;
            default:
              DUNE_THROW(Dune::Exception, "Element could not be assigned to a pre-defined group!");
          }

        }
        md_equilibrationGrid.preUpdateSubDomains();
        md_equilibrationGrid.updateSubDomains();
        md_equilibrationGrid.postUpdateSubDomains();
        debug_jochen << "Equilibration md grid restored. Transferring data to new grid..." << std::endl;

        // Now we have a fully reconstructed MultiDomainGrid; we can use the types from Acme2Traits now
        // to load the concentration DOFs and interpolate them onto the new grid!

        // Reconstruct pot GFS
        typename Acme2Traits::FEM_POT femPot;
        typename Acme2Traits::GFS_POT gfsPot_Equi_TEMP(equilibrationGV,femPot);
        //typename Acme2Traits::U_POT uoldPot_Equi(gfsPot_Equi,0.0);
        //typename Acme2Traits::U_POT unewPot_Equi(gfsPot_Equi,0.0);

        // Reconstruct con GFS
        typename Acme2Traits::FEM_CON femCon;
        typename Acme2Traits::GFS_SINGLE_CON gfsCon_Equi_Single_TEMP(elecGV_Equi,femCon);
        typename Acme2Traits::GFS_CON gfsCon_Equi_TEMP(gfsCon_Equi_Single_TEMP);
        //typename Acme2Traits::U_CON uoldCon_Equi(gfsCon_Equi,0.0);
        //typename Acme2Traits::U_CON unewCon_Equi(gfsCon_Equi,0.0);

        // Reconstruct MultiGFS
        typename Acme2Traits::MultiGFS multigfs_Equi(md_equilibrationGrid, gfsCon_Equi_TEMP, gfsPot_Equi_TEMP);

        // Extract sub GFS's
        typename Acme2Traits::GFS_CON_SUB gfsConSub_Equi(multigfs_Equi);
        typename Acme2Traits::GFS_POT_SUB gfsPotSub_Equi(multigfs_Equi);

        // Load solution vectors
        typename Acme2Traits::U uold_Equi(multigfs_Equi,0.0);
        typename Acme2Traits::U unew_Equi(multigfs_Equi,0.0);

        stdVec = simulationData.getVector("uold").data;
        uold_Equi.std_copy_from(stdVec);
        stdVec = simulationData.getVector("unew").data;
        unew_Equi.std_copy_from(stdVec);
//        stdVec = simulationData.getVector("uoldPot").data;
//        uoldPot_Equi.std_copy_from(stdVec);
//        stdVec = simulationData.getVector("unewPot").data;
//        unewPot_Equi.std_copy_from(stdVec);
//        stdVec = simulationData.getVector("uoldCon").data;
//        uoldCon_Equi.std_copy_from(stdVec);
//        stdVec = simulationData.getVector("unewCon").data;
//        unewCon_Equi.std_copy_from(stdVec);

        typename Acme2Traits::DGF_POT dgfOldPot_Equi(gfsPotSub_Equi, uold_Equi);
        typename Acme2Traits::DGF_POT dgfNewPot_Equi(gfsPotSub_Equi, unew_Equi);

        //typename Acme2Traits::DGF_CON dgfOldCon_Equi(gfsConOld, uoldCon_Equi);
        //typename Acme2Traits::DGF_CON dgfNewCon_Equi(gfsConOld, unewCon_Equi);
        typename Traits::DGF_CON dgfOldCon_Equi(gfsConSub_Equi, uold_Equi);
        typename Traits::DGF_CON dgfNewCon_Equi(gfsConSub_Equi, unew_Equi);

        Ax1CoarseGridInterpolationGridFunction<GV,typename Traits::DGF_POT,PHYSICS>
          dgfOldPot_FineGrid(gv,dgfOldPot_Equi,physics);
        Ax1CoarseGridInterpolationGridFunction<GV,typename Traits::DGF_POT,PHYSICS>
          dgfNewPot_FineGrid(gv,dgfNewPot_Equi,physics);

        // Define functions which interpolate the equilibration concentrations onto the new (finer) elecGV
        Ax1CoarseGridInterpolationGridFunction<SubGV,typename Traits::DGF_CON,PHYSICS>
          dgfOldCon_FineGrid(elecGV,dgfOldCon_Equi,physics);
        Ax1CoarseGridInterpolationGridFunction<SubGV,typename Traits::DGF_CON,PHYSICS>
          dgfNewCon_FineGrid(elecGV,dgfNewCon_Equi,physics);

        // Use the GFS belonging to the new (fine) grid here for interpolation!
        // This seems to work even when writing into the 'large' solution vectors
        Dune::PDELab::interpolate(dgfOldPot_FineGrid, gfsPot, solutionVectors.uold);
        Dune::PDELab::interpolate(dgfNewPot_FineGrid, gfsPot, solutionVectors.unew);
        Dune::PDELab::interpolate(dgfOldCon_FineGrid, gfsCon, solutionVectors.uold);
        Dune::PDELab::interpolate(dgfNewCon_FineGrid, gfsCon, solutionVectors.unew);

        debug_jochen << "Data transfer completed!" << std::endl;
      }

      time = simulationData.getTime();
      dt = simulationData.getDt();

      // If extrapolateXData==true, uold/unew are not touched!
      if(not extrapolateXData)
      {
        solutionVectors.deserialize(simulationData);
      }

      //Output::printSingleCoefficientVector(solutionVectors.uold, "uold");
      //Output::printSingleCoefficientVector(solutionVectors.unew, "unew");
      //Output::printSingleCoefficientVector(solutionVectors.uoldCon, "uoldCon");
      //Output::printSingleCoefficientVector(solutionVectors.unewCon, "unewCon");
      //Output::printSingleCoefficientVector(solutionVectors.uoldPot, "uoldPot");
      //Output::printSingleCoefficientVector(solutionVectors.unewPot, "unewPot");

      // simulationData is not valid after going out of this scope anyway!
      simulationData.clear();

      timer.stop();

    }


    /**
     * Call this to initialize a new output file, which can be used
     * to store data written out by the writeGFToGnuplot() method
     *
     * @param name
     */
    void initGnuplotFile(const std::string& name, const std::string header = "")
    {
      std::stringstream filename;
      filename << name << ".dat";
      Output::gnuplotInitialize(filename.str(), filePrefix, header);
    }

    /**
     * This is a method which other classes may use to write out non-standard gridfunctions.
     * This is flexible, but maybe a little more expensive (must check if file already exists,
     * gridfunction and solution vectors are created each time, etc.)
     *
     * Make sure to call initFile() once at the beginning to initialize the output file!
     *
     * @param gf
     * @param name
     */
    template<typename GF>
    void writeGFToGnuplot(GF& gf, const std::string& name, const std::string headerSuffix = "")
    {
      //debug_jochen << "Writing out vector for " << name << std::endl;
      std::stringstream filename;
      filename << name << ".dat";

      std::vector<typename GF::Traits::RangeType> sol;
      Ax1Output2D<Traits>::getSolutionVector(physics, gf, SUBSAMPLING_POINTS, x, sol);
      std::string header = infoStream.str() + headerSuffix;
      GnuplotTools2D::gnuplotAddBlock(filename.str(), x, sol, filePrefix, header);
    }

  private:
    void writePotentialOutput()
    {
      if(vtkOutput)
      {
        vtkwriterGV.addCellData(
            new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::DGF_POT>(dgfPot,"pot"));
        vtkwriterGV.addCellData(
            new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::DGF_POT_GRAD>(dgfPotGrad,"pot_grad"));
      }

      if(doFullGnuplotOutput)
      {
        Ax1Output2D<Traits>::getSolutionVector(physics, dgfPot, SUBSAMPLING_POINTS, x, pot);
        pot = physics.convertTo_mV(pot);

        // Get analytical solution for potential
        if(physics.getParams().hasAnalyticalSolution())
        {
          gfAnalyticalSolutionPot.setTime(time);
          gfAnalyticalSolutionCon.setTime(time);
          Ax1Output2D<Traits>::getSolutionVector(physics, gfAnalyticalSolutionPot, SUBSAMPLING_POINTS,
              x, solutionPot);

          std::vector<std::vector<typename Traits::DGF_POT::Traits::RangeType> > potVec
            = {pot, physics.convertTo_mV(solutionPot)};

          GnuplotTools2D::gnuplotAddMultiColumnBlock("pot.dat", x, potVec, filePrefix, infoStream.str());

          diagInfo.l2ErrorPot = ErrorNorms::l2Norm(dgfPot, gfAnalyticalSolutionPot, gv, intorderPot);
          diagInfo.maxErrorPot = ErrorNorms::maxNorm(dgfPot, gfAnalyticalSolutionPot, gv, intorderPot);

          debug_info << "== L2 error for potential: " << diagInfo.l2ErrorPot << std::endl;
          debug_info << "== Max error for potential: " << diagInfo.maxErrorPot << std::endl;

          //debug_info << "== L2 error for concentrations vector: "
          //    << ErrorNorms::l2Norm(dgfCon, gfAnalyticalSolutionCon, gfsCon.gridview()) << std::endl;
          //debug_info << "== Max error for concentrations vector: "
          //    << ErrorNorms::maxNorm(dgfCon, gfAnalyticalSolutionCon, gfsCon.gridview()) << std::endl;

        } else {
          GnuplotTools2D::gnuplotAddBlock("pot.dat", x, pot, filePrefix, infoStream.str());
        }

        // potential gradient
        Ax1Output2D<Traits>::getSolutionVector(physics, dgfPotGrad, SUBSAMPLING_POINTS, x, potGrad);
        GnuplotTools2D::gnuplotAddBlock("pot_grad.dat", x, potGrad, filePrefix, infoStream.str());
      }
    }

    void writeMembraneOutput()
    {
      if(physics.getParams().useMembrane())
      {
        if(vtkOutput)
        {
          vtkwriterSubGV.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_MEMB_POT>
            (gfMembranePotential,"memb_pot"));
        }

        // Membrane potential
        // TODO Adapt to new output strategy
        Ax1Output2D<Traits>::getCellCenterSolutionVector(physics, gfMembranePotential, x, membranePotential);
        membranePotential = physics.convertTo_mV(membranePotential);
        GnuplotTools2D::gnuplotAddLine("memb_pot.dat", time, membranePotential, filePrefix);

        // Membrane flux
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          std::vector<size_t> componentFlux;
            //debug_jochen << ">>> Getting flux component " << (j*GV::dimension + i) << std::endl;
          //componentFlux.push_back(j*GV::dimension);     // x component
          componentFlux.push_back(j*GV::dimension + 1); // y component

          Dune::shared_ptr<typename Traits::GF_SINGLE_MEMB_FLUX> gfSingleMembraneFlux(
              new typename Traits::GF_SINGLE_MEMB_FLUX(gfMembraneFlux, componentFlux));
          Dune::shared_ptr<typename Traits::GF_SINGLE_MEMB_FLUX_DIFF> gfSingleMembraneFlux_DiffTerm(
              new typename Traits::GF_SINGLE_MEMB_FLUX_DIFF(gfMembraneFlux_DiffTerm, componentFlux));
          Dune::shared_ptr<typename Traits::GF_SINGLE_MEMB_FLUX_DRIFT> gfSingleMembraneFlux_DriftTerm(
              new typename Traits::GF_SINGLE_MEMB_FLUX_DRIFT(gfMembraneFlux_DriftTerm, componentFlux));

          if(vtkOutput)
          {
            vtkwriterSubGV.addCellData(
                new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_SINGLE_MEMB_FLUX>
                        (gfSingleMembraneFlux,"memb_flux_" + physics.getIonName(j)));
            vtkwriterSubGV.addCellData(
                new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_SINGLE_MEMB_FLUX_DIFF>
                        (gfSingleMembraneFlux_DiffTerm,"memb_flux_diff_term_" + physics.getIonName(j)));
            vtkwriterSubGV.addCellData(
                new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_SINGLE_MEMB_FLUX_DRIFT>
                        (gfSingleMembraneFlux_DriftTerm,"memb_flux_drift_term_" + physics.getIonName(j)));
          }

          Ax1Output2D<Traits>::getCellCenterSolutionVector(physics, *gfSingleMembraneFlux, x, membraneFlux[j]);
          // Print total flux
          GnuplotTools2D::gnuplotAddLine("memb_flux_" + physics.getIonName(j) + ".dat", time, membraneFlux[j],
              filePrefix);

          if(doFullGnuplotOutput)
          {
            // Print diffusion/drift terms as well

            Ax1Output2D<Traits>::getCellCenterSolutionVector(physics, *gfSingleMembraneFlux_DiffTerm, x, membraneFlux_DiffTerm[j]);
            Ax1Output2D<Traits>::getCellCenterSolutionVector(physics, *gfSingleMembraneFlux_DriftTerm, x, membraneFlux_DriftTerm[j]);

            GnuplotTools2D::gnuplotAddLine("memb_flux_diff_term_" + physics.getIonName(j) + ".dat", time,
                membraneFlux_DiffTerm[j], filePrefix);
            GnuplotTools2D::gnuplotAddLine("memb_flux_drift_term_" + physics.getIonName(j) + ".dat", time,
                          membraneFlux_DriftTerm[j], filePrefix);
          }
        }

        // Channel stuff
        int nChannels = physics.getMembrane().getChannelSet().size();
        for(int k=0; k<nChannels; k++)
        {
          // Create grid function for the k-th channel
          Dune::shared_ptr<typename Traits::GF_CHANNEL> gfChannel(
              new typename Traits::GF_CHANNEL(membGV,physics,k));
          if(vtkOutput)
          {
            vtkwriterSubGV.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_CHANNEL>
                      (gfChannel,"channel_" + physics.getMembrane().getChannelSet().getChannel(k).getName()));
          }

          std::vector<typename Traits::GF_CHANNEL::Traits::RangeType> channelInfo;
          Ax1Output2D<Traits>::getCellCenterSolutionVector(physics, *gfChannel, x, channelInfo);


          //std::vector<Real> flatChannelInfo;
          //Tools::flattenVector(channelInfo, flatChannelInfo);
          //std::stringstream chStr;
          //chStr << std::setw(2) << std::setfill('0') << k;
          GnuplotTools2D::gnuplotAddLine("channel_" + physics.getMembrane().getChannelSet().getChannel(k).getName() + ".dat",
              time, channelInfo, filePrefix);
        }
      }
    }

    void writeConcentrationOutput()
    {
      if(doFullGnuplotOutput)
      {
        Ax1Output2D<Traits>::getSolutionVector(physics, gfChargeDensity, SUBSAMPLING_POINTS, x, chargeDensity);
        GnuplotTools2D::gnuplotAddBlock("cd.dat", x, chargeDensity, filePrefix, infoStream.str());
      }

      sumIonFluxIntegrals = 0.0;
      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        // select 1 component from the concentrations vector DGF
        std::vector<size_t> component;
        component.push_back(j);

        std::vector<size_t> componentGradient;
        for(int i=0; i<GV::dimension; i++)
        {
          componentGradient.push_back(j*GV::dimension + i);
        }

        Dune::shared_ptr<typename Traits::DGF_SINGLE_CON_MD> dgfSingleCon(
            new typename Traits::DGF_SINGLE_CON_MD(dgfCon, component));
        Dune::shared_ptr<typename Traits::DGF_SINGLE_CON_GRAD_MD> dgfSingleConGrad(
            new typename Traits::DGF_SINGLE_CON_GRAD_MD(dgfConGrad, componentGradient));
        Dune::shared_ptr<typename Traits::GF_SINGLE_ION_FLUX> gfSingleIonFlux(
            new typename Traits::GF_SINGLE_ION_FLUX(gfIonFlux, componentGradient));

         // VTK output
        if(vtkOutput)
        {
          vtkwriterGV.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::DGF_SINGLE_CON_MD>
            (dgfSingleCon,physics.getIonName(j)));
          vtkwriterGV.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::DGF_SINGLE_CON_GRAD_MD>
            (dgfSingleConGrad,physics.getIonName(j) + "_grad"));
          vtkwriterGV.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_SINGLE_ION_FLUX>
            (gfSingleIonFlux,"flux_" + physics.getIonName(j)));
        }

        if(doFullGnuplotOutput)
        {
          Ax1Output2D<Traits>::getSolutionVector(physics, *dgfSingleCon, SUBSAMPLING_POINTS, x, con[j]);
          // Get analytical solution for concentrations
          if(physics.getParams().hasAnalyticalSolution())
          {
            Dune::shared_ptr<typename Traits::SINGLE_ANALYTICAL_SOLUTION_CON> gfSingleAnalyticalSolutionCon(
                new typename Traits::SINGLE_ANALYTICAL_SOLUTION_CON(gfAnalyticalSolutionCon, component));
            Ax1Output2D<Traits>::getSolutionVector(
                physics, *gfSingleAnalyticalSolutionCon, SUBSAMPLING_POINTS, x, solutionCon[j]);

            diagInfo.setl2ErrorCon(j, ErrorNorms::l2Norm(*dgfSingleCon, *gfSingleAnalyticalSolutionCon,
                gv, intorderCon));
            diagInfo.setMaxErrorCon(j, ErrorNorms::maxNorm(*dgfSingleCon, *gfSingleAnalyticalSolutionCon,
                gv, intorderCon));
            debug_info << "== L2 error for " << physics.getIonName(j) << " concentration: "
                << diagInfo.getl2ErrorCon(j) << std::endl;
            debug_info << "== Max error for " << physics.getIonName(j) << " concentration: "
                << diagInfo.getMaxErrorCon(j) << std::endl;
          }
          Ax1Output2D<Traits>::getSolutionVector(physics, *dgfSingleConGrad, SUBSAMPLING_POINTS, x, conGrad[j]);

          Ax1Output2D<Traits>::getSolutionVector(physics, *gfSingleIonFlux, SUBSAMPLING_POINTS, x_elem, ionFlux[j]);

          // gnuplot
          if(physics.getParams().hasAnalyticalSolution())
          {
            std::vector<std::vector<typename Traits::SINGLE_ANALYTICAL_SOLUTION_CON::Traits::RangeType> > conVec
                      = {con[j], solutionCon[j]};
            GnuplotTools2D::gnuplotAddBlock(physics.getIonName(j) + ".dat", x, con[j], filePrefix, infoStream.str());
          }
          GnuplotTools2D::gnuplotAddBlock(physics.getIonName(j)+"_grad.dat", x, conGrad[j], filePrefix, infoStream.str());

          // Use (possibly) different position vector here!
          GnuplotTools2D::gnuplotAddBlock("flux_" + physics.getIonName(j)+".dat", x_elem, ionFlux[j],
              filePrefix, infoStream.str());

          //totalParticles[j] = Tools::integral( x, con[j] );
          typename Traits::DGF_SINGLE_CON_MD::Traits::RangeType concIntegral;
          Dune::PDELab::integrateGridFunction(*dgfSingleCon, concIntegral, intorderCon);
          if(time == 0)
          {
            totalInitParticles[j] = concIntegral;
          } else {
            totalParticles[j] = concIntegral;
          }

          std::vector<Real> chargeChange = {totalParticles[j] / totalInitParticles[j] - 1.0};
          GnuplotTools2D::gnuplotAddLine("con_error_" + physics.getIonName(j) + ".dat",
                                time, chargeChange, filePrefix);
        }

        //typename GF_SINGLE_IONFLUX::Traits::RangeType ionFluxIntegral;
        //Dune::PDELab::integrateGridFunction(*gfSingleIonFlux, ionFluxIntegral, intorderCon);
        //debug_verb << "Integral " << ION_NAMES[j] << " ion flux: " << ionFluxIntegral << std::endl;
        //sumIonFluxIntegrals += ionFluxIntegral;
      }
      //debug_verb << "SUM integrals ion flux: " << sumIonFluxIntegrals << std::endl;
    }

    // Constructor parameters
    const GV& gv;
    const SubGV& elecGV;
    const SubGV& membGV;

    MultiGFS& multigfs;
    GFS_CON& gfsCon;
    GFS_POT& gfsPot;

    Acme2SolutionVectors<Acme2Traits>& solutionVectors;
    U&       u;
    //U_CON&   uCon;
    //U_POT&   uPot;
    typename Traits::BGF_MEMB_FLUX&   bgfMembraneFlux;
    PHYSICS& physics;
    const int intorderPot;
    const int intorderCon;

    // (gnuplot) output arrays
    std::vector<DomainType>         x;
    std::vector<DomainType>         x_elem;

    std::vector<typename Traits::DGF_POT::Traits::RangeType> pot;
    std::vector<typename Traits::DGF_POT_GRAD::Traits::RangeType> potGrad;
    std::vector<typename Traits::GF_CD::Traits::RangeType> chargeDensity;
    std::vector<typename Traits::GF_MEMB_POT::Traits::RangeType> membranePotential;
    std::vector<std::vector<typename Traits::GF_SINGLE_MEMB_FLUX::Traits::RangeType> > membraneFlux;
    std::vector<std::vector<typename Traits::GF_SINGLE_MEMB_FLUX_DIFF::Traits::RangeType> > membraneFlux_DiffTerm;
    std::vector<std::vector<typename Traits::GF_SINGLE_MEMB_FLUX_DRIFT::Traits::RangeType> > membraneFlux_DriftTerm;
    std::vector<std::vector<typename Traits::DGF_SINGLE_CON_MD::Traits::RangeType> > con;
    std::vector<std::vector<typename Traits::DGF_SINGLE_CON_GRAD_MD::Traits::RangeType> > conGrad;
    std::vector<std::vector<typename Traits::GF_SINGLE_ION_FLUX::Traits::RangeType> > ionFlux;

    std::vector<std::vector<typename Traits::SINGLE_ANALYTICAL_SOLUTION_CON::Traits::RangeType> > solutionCon;
    std::vector<typename Traits::ANALYTICAL_SOLUTION_POT::Traits::RangeType> solutionPot;

    std::vector<Real>                 totalInitParticles;
    std::vector<Real>                 totalParticles;

    Real                              sumIonFluxIntegrals;

    Dune::SubsamplingVTKWriter<GV>     vtkwriterGV;
    Dune::SubsamplingVTKWriter<SubGV>  vtkwriterSubGV;

    const std::string vtkBasenameDomain;
    const std::string vtkBasenameMembrane;
    const std::string checkpointBasename;
    Dune::PDELab::FilenameHelper      fn_vtk_GV;
    Dune::PDELab::FilenameHelper      fn_vtk_SubGV;
    Dune::PDELab::FilenameHelper      fn_checkpoints;

    const std::string filePrefix;

    typename Traits::DGF_CON         dgfConElec;
    typename Traits::DGF_CON_GRAD    dgfConGradElec;
    typename Traits::DGF_CON_MD      dgfCon;
    typename Traits::DGF_CON_GRAD_MD dgfConGrad;

    typename Traits::DGF_POT           dgfPot;
    typename Traits::DGF_POT_GRAD      dgfPotGrad;
    typename Traits::DGF_POT_GRAD_ELEC dgfPotGradElec;

    typename Traits::GF_CD           gfChargeDensity;
    typename Traits::GF_ELEC_FLUX    gfElecFlux;

    typename Traits::GF_MEMB_POT_MD  gfMembranePotentialMD;
    typename Traits::GF_MEMB_POT     gfMembranePotential;

    typename Traits::GF_MEMB_FLUX       gfMembraneFlux;
    typename Traits::BGF_MEMB_FLUX_DIFF  bgfMembraneFlux_DiffTerm;
    typename Traits::GF_MEMB_FLUX_DIFF  gfMembraneFlux_DiffTerm;
    typename Traits::BGF_MEMB_FLUX_DRIFT bgfMembraneFlux_DriftTerm;
    typename Traits::GF_MEMB_FLUX_DRIFT gfMembraneFlux_DriftTerm;

    typename Traits::GF_ION_FLUX     gfIonFlux;

    typename Traits::ANALYTICAL_SOLUTION_CON    gfAnalyticalSolutionCon;
    typename Traits::ANALYTICAL_SOLUTION_POT    gfAnalyticalSolutionPot;

    const bool vtkOutput;
    const bool doFullGnuplotOutput;

    Dune::Timer timer;

    DiagnosticInfo diagInfo;

    double time;
    std::stringstream infoStream;

};

#endif /* DUNE_AX1_ACME2_OUTPUT_HH */
