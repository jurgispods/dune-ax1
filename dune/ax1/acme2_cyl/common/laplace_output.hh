/*
 * acme2_cyl_output.hh
 *
 *  Created on: Aug 11, 2011
 *      Author: jpods
 */
#ifndef DUNE_AX1_LAPLACE_OUTPUT_HH
#define DUNE_AX1_LAPLACE_OUTPUT_HH


#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/pdelab/common/clock.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/common/functionutilities.hh>
#include <dune/pdelab/function/selectcomponent.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>

#if HAVE_HDF5
#include <dune/ax1/common/hdf5_tools.hh>
#endif
#include <dune/ax1/common/ax1_gridfunctions.hh>

#include <dune/ax1/acme2_cyl/common/laplace_simloader.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_solutionvectors.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_geometrytools.hh>



template<class Acme2CylTraits>
class LaplaceOutput
{
  public:

    // Boilerplate code to ease type usage
    typedef typename Acme2CylTraits::Coord Coord;
    typedef typename Acme2CylTraits::Real Real;
    typedef typename Acme2CylTraits::GridView GV;
    typedef typename Acme2CylTraits::SubGridView SubGV;
    typedef typename Acme2CylTraits::MultiGFS MultiGFS;
    typedef typename Acme2CylTraits::GFS_POT_SUB GFS_POT; // SubGFS!
    typedef typename Acme2CylTraits::U U;
    typedef typename Acme2CylTraits::Physics PHYSICS;

    const int SUBSAMPLING_POINTS;

    const std::string vtkOutputBase_Domain;
    const std::string vtkOutputBase_Membrane;

    // A bunch of typedef's for all the gridfunctions
    struct LaplaceOutputTraits
    {
      // DGFs
      typedef Dune::PDELab::DiscreteGridFunction        <GFS_POT,U> DGF_POT;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_POT,U> DGF_POT_GRAD;
      typedef Ax1MultiDomainGridFunctionRestriction<SubGV,DGF_POT_GRAD,PHYSICS> DGF_POT_GRAD_ELEC;

      typedef MembranePotentialGridFunction<DGF_POT, PHYSICS> GF_MEMB_POT_MD;
      typedef Ax1MultiDomainGridFunctionRestriction<SubGV,GF_MEMB_POT_MD,PHYSICS> GF_MEMB_POT;
      typedef MembranePotentialGridFunction_MembGV<SubGV,DGF_POT,PHYSICS> GF_MEMB_POT_MEMBGV;

      typedef Ax1MembraneElementIndexGridFunction<GV,PHYSICS> GF_MEMB_GROUPS;

      // Geometry helper gridfunctions
      typedef Acme2CylPlain2DAreaGridFunction<GV,Real> GF_PLAIN2D_AREA;
      typedef Acme2CylVolumeGridFunction<GV,Real> GF_VOLUME;
      typedef Acme2CylBaseSurfaceGridFunction<GV,Real> GF_BASE_SURFACE;
      typedef Acme2CylSideSurfaceGridFunction<GV,Real> GF_SIDE_SURFACE;
      typedef Acme2CylSurfaceGridFunction<GV,Real> GF_SURFACE;
      typedef Ax1PermittivityGridFunction<GV,PHYSICS> GF_PERMITTIVITY;
      typedef Ax1ProcessorGridFunction<GV> GF_PARTITION;

      // Boundary value functions
      typedef Ax1BoundaryValueFunction<DGF_POT,PHYSICS> GF_POT_BOUNDARY;
    };

    typedef LaplaceOutputTraits Traits;
    typedef typename Traits::DGF_POT::Traits::DomainType DomainType;

    struct DiagnosticInfo
    {
      DiagnosticInfo(bool hasAnalyticalSolution_)
      : time(0.0),
        iterations(0), dt(0.0),
        l2ErrorPot(0.0), maxErrorPot(0.0),
        l2ErrorCon(NUMBER_OF_SPECIES,0.0), maxErrorCon(NUMBER_OF_SPECIES,0.0),
        tEquilibrium(0.0),
        maxDiffCon(0),
        maxDiffPot(0),
        hasAnalyticalSolution(hasAnalyticalSolution_)
      {}

      enum ReductionFunction { Function_Identity = 0, Function_Min = 1, Function_Max = 2,
        Function_Mean = 3, Function_Sum = 4  };

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
      bool hasAnalyticalSolution;

      std::map<std::string, Real> debugData;
      std::map<std::string, int> debugDataReductionFunctions;


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
        int count = (hasAnalyticalSolution ? 6 : 4);
        if(hasAnalyticalSolution)
        {
          header << "#  (1)time        (2)dt             (3)#iterations   (4)[pot] L2 error (5)[pot] max error";
          for(int i=0; i<NUMBER_OF_SPECIES; ++i)
          {
            header << " (" << (count++) << ")[" << ION_NAMES[i] << "] L2 error";
            header << " (" << (count++) << ")[" << ION_NAMES[i] << "] max error";
          }
        } else {
          header << "#  (1)time        (2)dt             (3)#iterations";
        }
        for(typename std::map<std::string,Real>::const_iterator it = debugData.begin();
              it != debugData.end(); ++it)
        {
          header << "   (" << (count++) << ")" << it->first;
        }
        return header.str();
      }

      void registerDebugData(std::string key, Real value = std::numeric_limits<Real>::lowest(),
          const int reductionFunction = ReductionFunction::Function_Identity)
      {
        debugData[key] = value;
        debugDataReductionFunctions[key] = reductionFunction;
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
        if(hasAnalyticalSolution)
        {
          vec.push_back(l2ErrorPot);
          vec.push_back(maxErrorPot);
          for(int i=0; i<l2ErrorCon.size(); i++)
          {
            vec.push_back(l2ErrorCon[i]);
            vec.push_back(maxErrorCon[i]);
          }
        }
        for(typename std::map<std::string,Real>::const_iterator it = debugData.begin();
            it != debugData.end(); ++it)
        {
          vec.push_back(it->second);
        }
        return vec;
      }

      //! This returns a vector of the same size as getVectorRepresentation containing a flag
      //! for each diagnostic data entry which function to be used when reducing the values
      //! on all processors to the root node.
      std::vector<int> getReductionFunctions() const
      {
        int count = 0;

        std::vector<int> vec;
        //vec[0] = time;
        vec.push_back(ReductionFunction::Function_Identity); // dt
        vec.push_back(ReductionFunction::Function_Identity); // iterations
        if(hasAnalyticalSolution)
        {
          vec.push_back(ReductionFunction::Function_Identity); // l2ErrorPot
          vec.push_back(ReductionFunction::Function_Max); // maxErrorPot
          for(int i=0; i<l2ErrorCon.size(); i++)
          {
            vec.push_back(ReductionFunction::Function_Identity); // l2ErrorCon[i]
            vec.push_back(ReductionFunction::Function_Max); // maxErrorCon[i]
          }
        }
        for(typename std::map<std::string,int>::const_iterator it = debugDataReductionFunctions.begin();
            it != debugDataReductionFunctions.end(); ++it)
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

        for(typename std::map<std::string, Real>::iterator it = debugData.begin();
            it != debugData.end(); ++it)
        {
          it->second = std::numeric_limits<Real>::lowest();
        }
      }
    };

    LaplaceOutput(const GV& gv_, const SubGV& elecGV_, const SubGV& membGV_,
        MultiGFS& multigfs_, GFS_POT& gfsPot_,
        Acme2CylSolutionVectors<Acme2CylTraits>& solutionVectors_,
        PHYSICS& physics_,
        int intorderPot_=2) :
      SUBSAMPLING_POINTS(0),
      vtkOutputBase_Domain("acme2_cyl"),
      vtkOutputBase_Membrane("acme2_cyl_membrane"),
      gv(gv_),
      elecGV(elecGV_),
      membGV(membGV_),
      multigfs(multigfs_),
      gfsPot(gfsPot_),
      solutionVectors(solutionVectors_),
      u(solutionVectors_.unew),
      physics(physics_),
      intorderPot(intorderPot_),
      vtkwriterGV(gv,physics_.getParams().vtkSubSamplingLevel()),
      vtkwriterSubGV(membGV,0), // Don't do subsampling on the membrane!
      vtkBasenameDomain(physics.getParams().getOutputPrefix() + "paraview/" + vtkOutputBase_Domain),
      vtkBasenameMembrane(physics.getParams().getOutputPrefix() + "paraview/" + vtkOutputBase_Membrane),
      checkpointBasename(physics.getParams().getOutputPrefix() + "checkpoints/checkpoint"),
      hdf5Basename(physics.getParams().getOutputPrefix() + "hdf5/output"),
      fn_vtk_GV(vtkBasenameDomain),
      fn_vtk_SubGV(vtkBasenameMembrane),
      fn_checkpoints(checkpointBasename,-1), // checkpoints counter is always 1 lower than the others
      fn_hdf5(hdf5Basename),
      filePrefix(physics.getParams().getOutputPrefix()),
      dgfPot(gfsPot,u),
      dgfPotGrad(gfsPot,u),
      dgfPotGradElec(elecGV,dgfPotGrad,physics),
      gfMembranePotentialMD(dgfPot, physics),
      gfMembranePotential(membGV, gfMembranePotentialMD, physics),
      vtkOutput(physics.getParams().doVTKOutput()),
      doFullGnuplotOutput(physics.getParams().doFullGnuplotOutput()),
      doHDF5Output(physics.getParams().doHDF5Output()),
      hdf5FileInitialized(false),
      to_mV(physics),
      timer(false),
      diagInfo(physics.getParams().hasAnalyticalSolution()),
      nTimeSteps(-1)
    {
      timer.start();

      // All standard file I/O is handled by rank 0 only!
      // Exceptions: HDF5 output (MPI-parallel) and boundary value output (also rank size-1 for right boundary)
      if(gv.comm().rank() == 0)
      {
        std::vector<std::string> outputDirs = {"paraview", "checkpoints", "config", "hdf5"};
        int existingVTK = 0;
        int existingCheckpoints = 0;
        int existingHDF5 = 0;

        // Set up list of gnuplot files
        gnuplotFiles.push_back("cd.dat");
        gnuplotFiles.push_back("pot.dat");
        gnuplotFiles.push_back("pot_grad.dat");
        gnuplotFiles.push_back("memb_pot.dat");

        gnuplotFiles.push_back("boundary_pot_bottom.dat");
        gnuplotFiles.push_back("boundary_pot_top.dat");
        gnuplotFiles.push_back("boundary_pot_left.dat");
        gnuplotFiles.push_back("boundary_pot_right.dat");

        // 'One time' output stuff (only done once in the beginning)
        oneTimeGnuplotFiles.push_back("area_2d.dat");
        oneTimeGnuplotFiles.push_back("permittivity.dat");
        oneTimeGnuplotFiles.push_back("partition.dat");
        oneTimeGnuplotFiles.push_back("memb_groups.dat");
        oneTimeGnuplotFiles.push_back("volume.dat");
        oneTimeGnuplotFiles.push_back("base_surface.dat");
        oneTimeGnuplotFiles.push_back("side_surface.dat");
        oneTimeGnuplotFiles.push_back("surface.dat");


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

          // Remove old VTK files
          command = "rm " + hdf5Basename + "*";
          debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
          unusedVariable = std::system(command.c_str());

          // Remove old files from config directory
          command = "rm " + physics.getParams().getOutputPrefix() + "config/*";
          debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
          unusedVariable = std::system(command.c_str());

          // Following: Initialize (i.e., create and write headers) all gnuplot output files
          for(std::vector<std::string>::const_iterator it = gnuplotFiles.begin(); it != gnuplotFiles.end(); ++it )
          {
            Output::gnuplotInitialize(*it, filePrefix );
          }
        } else {
          // Continue simulation => keep old output files and count

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
                existingVTK++;
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
              // Count number of checkpoint files generated by rank 0 only, each processor creates
              // its own checkpoint file!
              if(filename.find("checkpoint") != std::string::npos && filename.find("_p0") != std::string::npos)
              {
                existingCheckpoints++;
              }
            }
            closedir (dir);
          } else {
            DUNE_THROW(Dune::IOError,
              "Error reading directory '" + (physics.getParams().getOutputPrefix() + "checkpoints")
              + "'");
          }

          // Determine the number of written output checkpoint files
          dir = opendir ((physics.getParams().getOutputPrefix() + "hdf5").c_str());
          if (dir != NULL) {

            /* print all the files and directories within directory */
            while ((ent = readdir(dir)) != NULL)
            {
              std::string filename(ent->d_name);
              if(filename.find("output") != std::string::npos)
              {
                existingHDF5++;
              }
            }
            closedir (dir);
          } else {
            DUNE_THROW(Dune::IOError,
              "Error reading directory '" + (physics.getParams().getOutputPrefix() + "hdf5")
              + "'");
          }

          std::string command;
          // Backup one-time gnuplot files, as these are not handled by SimulationLoader::truncateAndDeleteExistingFiles()
          for(std::vector<std::string>::const_iterator it = oneTimeGnuplotFiles.begin(); it != oneTimeGnuplotFiles.end(); ++it )
          {
            command += "mv " + filePrefix + *it + " " + filePrefix + *it + ".old; ";
          }
          debug_jochen << "Backing up one-time gnuplot files with command '" << command << "'" << std::endl;
          int unusedVariable = std::system(command.c_str());


          // Calculate number of written time steps as the maximum of all three variants
          nExistingOutputFiles = std::max(existingVTK-1, std::max(existingCheckpoints, existingHDF5-1));
        }

        // One-time gnuplot files are always initialized (i.e. truncated), no matter if simulation is continued or not
        for(std::vector<std::string>::const_iterator it = oneTimeGnuplotFiles.begin(); it != oneTimeGnuplotFiles.end(); ++it )
        {
         Output::gnuplotInitialize(*it, filePrefix );
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

        // Generate detailed git version control info for this module
        command = "bash git-info.bash && mv --backup=t git-info.dat "
            + physics.getParams().getOutputPrefix() + "config";
        debug_jochen << "Saving Git VC info using command '" << command << "'" << std::endl;
        unusedVariable = std::system(command.c_str());

        // Generate version control info for all Dune modules
        command = "bash dune-info.bash && mv --backup=t dune-info.dat "
            + physics.getParams().getOutputPrefix() + "config";
        debug_jochen << "Saving Dune VC info using command '" << command << "'" << std::endl;
        unusedVariable = std::system(command.c_str());

        // Create a DGF file from the multidomain GV
        if(physics.getParams().createGridFile() || physics.getParams().doEquilibration())
        {
          debug_jochen << "Saving DGF file..." << std::endl;
          Dune::DGFWriter<GV> dgfWriter(gv);

          std::string filename = physics.getParams().getOutputPrefix()
              + "config/" + physics.getParams().getSaveGridFileName();
          dgfWriter.write(filename);

#if USE_GRID == 1
          debug_jochen << "Saving special YaspGrid DGF file..." << std::endl;
          // Write out 'INTERVAL'-style DGF file in case of Yasp base grid
          int index = filename.find(".dgf");
          assert(index != std::string::npos);

          std::string basegrid_filename(filename.substr(0,index));
          basegrid_filename += "_basegrid.dgf";

          std::ofstream out_basegrid(basegrid_filename);
          if(! out_basegrid.good())
            DUNE_THROW(Dune::Exception, "DGF file '" << basegrid_filename << "' could not be created!");

          out_basegrid << "DGF" << std::endl;
          out_basegrid << std::endl;
          out_basegrid << "INTERVAL" << std::endl;
          out_basegrid << "0 0" << std::endl;
          out_basegrid << physics.getParams().xMax() << " " << physics.getParams().yMax() << std::endl;
          out_basegrid << (physics.getParams().X().size() - 1) << " "
                       << (physics.getParams().Y().size() - 1) << std::endl;
          out_basegrid << "# INTERVAL" << std::endl;
          out_basegrid << std::endl;
          out_basegrid << "GRIDPARAMETER" << std::endl;
          out_basegrid << "name Yasp_BaseGrid" << std::endl;
          out_basegrid << "overlap " << gv.overlapSize(0) << std::endl;
          out_basegrid << "# GRIDPARAMETER" << std::endl;
          out_basegrid << std::endl;
          out_basegrid << "# DGF";
#endif
        }
      }
      debug_jochen << "Writing one-time gridfunction output..." << std::endl;

      // Write out volume gf's (this is done only once in the beginning)
      typename Traits::GF_PLAIN2D_AREA area2DGf(gv);
      GnuplotTools2D<Traits>::writeGF(physics, area2DGf, SUBSAMPLING_POINTS,
          filePrefix + "area_2d.dat", time, "plain 2D element area", true);

      // Write permittivity (gnuplot debug output)
      typename Traits::GF_PERMITTIVITY permittivityGf(gv,physics);
      GnuplotTools2D<Traits>::writeGF(physics, permittivityGf, SUBSAMPLING_POINTS,
          filePrefix + "permittivity.dat", time, "permittivity", true);

      typename Traits::GF_PARTITION processorGf(gv);
      GnuplotTools2D<Traits>::writeGF(physics, processorGf, SUBSAMPLING_POINTS,
            filePrefix + "partition.dat", time, "partition", true);

      if(physics.getParams().useMembrane())
      {
        // RangeType is tuple!
        typename Traits::GF_MEMB_GROUPS membGroupGf(gv,physics);

        std::string membGroupFilename = physics.getParams().getOutputPrefix() + "config/membrane_groups.dat";
        std::string infoMembGroups = "Membrane elements (1,2) coordinates (3) membrane index (4) group index (5) group name "
             "(6) permittivity (7) processor";
        infoMembGroups += "\n# " + physics.getGroupInfo();
        GnuplotTools2D<Traits>::writeGF(physics, membGroupGf, SUBSAMPLING_POINTS,
              filePrefix + "memb_groups.dat", time, infoMembGroups, true);
      }

#if HAVE_HDF5
      if(doHDF5Output)
      {
        std::string filename(physics.getParams().getOutputPrefix() + "hdf5/geometry.h5");
        HDF5Tools<Traits>::initialize(filename, 0.0,physics);
        std::string datasetName("area_2d");
        HDF5Tools<Traits>::writeGF(physics,area2DGf,SUBSAMPLING_POINTS,
            filename,datasetName);
        datasetName = "permittivity";
        HDF5Tools<Traits>::writeGF(physics,permittivityGf,SUBSAMPLING_POINTS,
            filename,datasetName);
        datasetName = "partition";
        HDF5Tools<Traits>::writeGF(physics,processorGf,SUBSAMPLING_POINTS,
            filename,datasetName);
        if(physics.getParams().useMembrane())
        {
          // This will currently not work, as the gridfunction returns a tuple of different data types
//          datasetName = "membrane_element_indices";
//          HDF5Tools<Traits>::writeGF(physics,membElementIndexGf,SUBSAMPLING_POINTS,
//              filename,datasetName);
        }
      }
#endif

      // Helper grid functions for cylinder volume/surface areas etc.
      if(USE_CYLINDER_COORDINATES)
      {
        typename Traits::GF_VOLUME volumeGf(gv);
        GnuplotTools2D<Traits>::writeGF(physics, volumeGf, SUBSAMPLING_POINTS,
            filePrefix + "volume.dat", time, "cylinder element volume", true);

        typename Traits::GF_BASE_SURFACE baseSurfaceGf(gv);
        GnuplotTools2D<Traits>::writeGF(physics, baseSurfaceGf, SUBSAMPLING_POINTS,
            filePrefix + "base_surface.dat", time, "cylinder element base (top+bottom) areas", true);

        typename Traits::GF_SIDE_SURFACE sideSurfaceGf(gv);
        GnuplotTools2D<Traits>::writeGF(physics, sideSurfaceGf, SUBSAMPLING_POINTS,
            filePrefix + "side_surface.dat", time, "cylinder element side surface area", true);

        typename Traits::GF_SURFACE surfaceGf(gv);
        GnuplotTools2D<Traits>::writeGF(physics, surfaceGf, SUBSAMPLING_POINTS,
            filePrefix + "surface.dat", time, "total cylinder element surface", true);

#if HAVE_HDF5
        if(doHDF5Output)
        {
          std::string filename(physics.getParams().getOutputPrefix() + "hdf5/geometry.h5");
          std::string datasetName = "volume";
          HDF5Tools<Traits>::writeGF(physics,volumeGf,SUBSAMPLING_POINTS,
              filename,datasetName);
          datasetName = "base_surface";
          HDF5Tools<Traits>::writeGF(physics,baseSurfaceGf,SUBSAMPLING_POINTS,
              filename,datasetName);
          datasetName = "side_surface";
          HDF5Tools<Traits>::writeGF(physics,sideSurfaceGf,SUBSAMPLING_POINTS,
              filename,datasetName);
          datasetName = "surface";
          HDF5Tools<Traits>::writeGF(physics,surfaceGf,SUBSAMPLING_POINTS,
              filename,datasetName);
        }
#endif
      }

      debug_jochen << "...done!" << std::endl;
      timer.stop();
    }

    void writeStep(double time_)
    {
      double start = timer.elapsed();
      double previous = start;
      double elapsed = start;
      timer.start();

      nTimeSteps++;
      time = time_;
      infoStream.str("");
      infoStream << "time: " << std::setprecision(16) << time;
      diagInfo.time = time;

#if HAVE_HDF5
      if(doHDF5Output)
      {
        if(! hdf5FileInitialized)
        {
          debug_jochen << "Initializing HDF file '" << getHDF5FileName() << "'" << std::endl;
          HDF5Tools<Traits>::initialize(getHDF5FileName(),time,physics);
          hdf5FileInitialized = true;
        }
        debug_jochen << "Writing to HDF5 file '" << getHDF5FileName() << "'" << std::endl;
      }
#endif

      previous = elapsed;
      elapsed = timer.elapsed();
      debug_info << "[Acme2CylOutput] Output preparation time: " << (elapsed-previous) << " s" << std::endl;

      writePotentialOutput();
      previous = elapsed;
      elapsed = timer.elapsed();
      debug_info << "[Acme2CylOutput] Potential output time: " << (elapsed-previous) << " s" << std::endl;

      // Only do membrane output if there is a membrane!
      if(physics.getParams().useMembrane())
      {
        writeMembraneOutput();
        previous = elapsed;
        elapsed = timer.elapsed();
        debug_info << "[Acme2CylOutput] Membrane output time: " << (elapsed-previous) << " s" << std::endl;
      }

      writeBoundaryOutput();
      previous = elapsed;
      elapsed = timer.elapsed();
      debug_info << "[Acme2CylOutput] Boundary output time: " << (elapsed-previous) << " s" << std::endl;

      writeDiagnosticOutput();
      previous = elapsed;
      elapsed = timer.elapsed();
      debug_info << "[Acme2CylOutput] Diagnostic output time: " << (elapsed-previous) << " s" << std::endl;

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
      if(doHDF5Output)
      {
        fn_hdf5.increment();
        hdf5FileInitialized = false;
      }
      if(physics.getParams().doCheckpointing())
      {
        fn_checkpoints.increment();
      }

      elapsed = timer.stop() - start;
      debug_info << "[Acme2CylOutput] TOTAL output time: " << elapsed << " s" << std::endl;
    }

    DiagnosticInfo& getDiagInfo()
    {
      return diagInfo;
    }

    void initDiagInfoFile()
    {
      if(gv.comm().rank() == 0)
      {
        Output::gnuplotInitialize(physics.getParams().getDiagnosticsFilename(),
          physics.getParams().getOutputPrefix(), diagInfo.getHeader());
      }
    }


    void printAllCoeffs()
    {
      Output::printSingleCoefficientVector(u, "u");
    }

    double getTotalOutputTime()
    {
      return timer.elapsed();
    }

    int getNTimeSteps() const
    {
      return nTimeSteps;
    }

    /**
     * Function to write out a complete set of solution vectors which can be used to continue a
     * simulation later or to set a defined initial state for a new simulation
     *
     * @param time
     * @param dt
     * @param filename
     */
    void saveState(double time, double dt, std::string& filename)
    {
      timer.start();

      // Empty filename is given => this is a checkpoint output!
      if(filename.size() == 0)
      {
        std::stringstream checkpoint_filename;
        checkpoint_filename << fn_checkpoints.getName()
          << "_p" << physics.gridView().comm().rank()
          << ".dat";
        filename = checkpoint_filename.str();
      } else {
        int index = filename.find(".dat");
        assert(index != std::string::npos);
        std::stringstream save_filename;
        save_filename << filename.substr(0,index)
            << "_p" << physics.gridView().comm().rank()
            << ".dat";
        filename = save_filename.str();
      }

      Ax1SimulationData<Real> simulationData(gv.size(0), time, dt, filename);
      simulationData.setTimeStep(nTimeSteps);
      simulationData.setMpiRank(gv.comm().rank());
      simulationData.setMpiSize(gv.comm().size());

      // NEW method: save separate grid vectors for x,y coordinates
      simulationData.addVector("x", physics.getParams().X());
      simulationData.addVector("y", physics.getParams().Y());

      // Fill simulation data with std::vectors containing data from solutionVectors
      solutionVectors.serialize(simulationData);

      Ax1SimulationState<Real> state;
      state.setData(simulationData);
      state.saveState();

      timer.stop();
    }

    std::map<std::string,std::string> loadState(double& time, double& dt)
    {
      timer.start();

      // This class contains several helper functions when loading previously saved simulation data
      LaplaceSimulationLoader<Acme2CylTraits, Traits> simLoader(physics, gv, elecGV, gfsPot, solutionVectors);

      bool extrapolateXData = physics.getParams().doExtrapolateXData();

      std::map<std::string,std::string> loadFileNames = simLoader.getLoadFileNames();

      if(loadFileNames.size() > 2)
        DUNE_THROW(Dune::Exception, "Loading from more than 2 simulation states not implemented yet!");
      if(loadFileNames.size() < 1)
        DUNE_THROW(Dune::Exception, "Did not find any 'loadFilename' in config file!");

      bool doAllCoordinatesMatch = true;

      // Store all states we possibly load from in a vector and perform checks an them
      std::vector<Ax1SimulationState<Real> > states(loadFileNames.size());
      int i=0;
      for(std::map<std::string,std::string>::iterator it = loadFileNames.begin(); it != loadFileNames.end(); ++it)
      {
        debug_jochen << "Processing loadFilename " << it->second << " with group name '" << it->first << "'" << std::endl;

        // Get consistent parallel filename
        simLoader.getParallelFilename(it->second);

        // Load the simulation state and store it in a vector
        states[i].loadState(it->second);

        debug_jochen << "Ok, loaded state: nElems = " << states[i].getData().getNElements() << std::endl;

        states[i].getData().addMetaInfo("groupName", it->first);

        // Check coordinates of the underlying grid; is it compatible?
        bool doCoordinatesMatch = simLoader.checkCoordinates(states[i].getData());

        doAllCoordinatesMatch = doAllCoordinatesMatch && doCoordinatesMatch;

        assert(doAllCoordinatesMatch == doCoordinatesMatch);

        i++;
      }


      // Get simulation data from first entry in states vector
      Ax1SimulationData<Real>& simulationData = states[0].getData();

      if(extrapolateXData && !(doAllCoordinatesMatch) )
      {
        // Load from vector of simulation states!
        simLoader.extrapolateData(states);

      } else {
        // If data didn't need to be extrapolated, simply copy loaded solution vectors!
        solutionVectors.deserialize(simulationData);
      }

      time = simulationData.getTime();
      dt = simulationData.getDt();

      // Special case 'continueSimulation': Keep all data generated until the timestep we are loading from!
      if(physics.getParams().doLoadState() && physics.getParams().doContinueSimulation())
      {
        int loadTimeStep = simulationData.getTimeStep();

        // Set the timestep counter to the loaded value
        nTimeSteps = loadTimeStep;

        debug_jochen << "Continue simulation from timestep #" << loadTimeStep << " (t=" << time
            << ")" << std::endl;

        // Calculate number of written time steps as the maximum of all three variants
        nExistingOutputFiles = std::max(loadTimeStep, nExistingOutputFiles);


        std::vector<std::string> gnuplotFilesToBeTruncated = gnuplotFiles;
        gnuplotFilesToBeTruncated.push_back(physics.getParams().getDiagnosticsFilename());
        // Truncate gnuplot files / delete spurious checkpoint and HDF5 files
        simLoader.truncateAndDeleteExistingFiles(loadTimeStep, time, nExistingOutputFiles,
            gnuplotFilesToBeTruncated, vtkBasenameDomain, vtkBasenameMembrane,
            checkpointBasename, hdf5Basename, filePrefix);

        // Set filename helper counters to the correct value; do this for all processors!
        for(int i=0; i<loadTimeStep+1; i++)
        {
          fn_vtk_GV.increment();
          fn_vtk_SubGV.increment();
          fn_checkpoints.increment();
          fn_hdf5.increment();
        }
      }

      // simulationData is not valid after going out of this scope anyway!
      simulationData.clear();

      timer.stop();

      return loadFileNames;
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
    void writeGFToGnuplot(GF& gf, const std::string& name, const std::string headerSuffix = "",
        const bool forceWrite = false)
    {
      GnuplotTools2D<Traits>::writeGF(physics, gf, SUBSAMPLING_POINTS, filePrefix + name + ".dat",
          time, infoStream.str() + headerSuffix, forceWrite);

#if HAVE_HDF5
      if(doHDF5Output)
      {
        if(! hdf5FileInitialized)
        {
          // As this method is called before writeStep(), the HDF5 has to be created first
          HDF5Tools<Traits>::initialize(getHDF5FileName(), diagInfo.time, physics);
          hdf5FileInitialized = true;
        }

        HDF5Tools<Traits>::writeGF(physics,gf,SUBSAMPLING_POINTS,getHDF5FileName(),name);
      }
#endif
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
        GnuplotTools2D<Traits>::writeGF(physics, dgfPot, SUBSAMPLING_POINTS, filePrefix + "pot.dat",
            time, infoStream.str(), to_mV);

        // potential gradient
        GnuplotTools2D<Traits>::writeGF(physics, dgfPotGrad, SUBSAMPLING_POINTS, filePrefix + "pot_grad.dat",
            time, infoStream.str());
      }

#if HAVE_HDF5
      if(doHDF5Output)
      {
        std::string datasetName("pot");
        HDF5Tools<Traits>::writeGF(physics,dgfPot,SUBSAMPLING_POINTS,
            getHDF5FileName(),datasetName, to_mV);
        datasetName = "pot_grad";
        HDF5Tools<Traits>::writeGF(physics,dgfPotGrad,SUBSAMPLING_POINTS,
            getHDF5FileName(),datasetName); // Do not convert gradient to mV, as
                                            // Matlab evaluation scripts already do the conversion
      }
#endif
    }

    void writeMembraneOutput()
    {
      if(physics.getParams().useMembrane())
      {
        if(vtkOutput)
        {
          DUNE_THROW(Dune::Exception, "VTK output for membrane gridfunctions not implemented!");
          //vtkwriterSubGV.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_MEMB_POT>
          //  (gfMembranePotential,"memb_pot"));
        }

        GnuplotTools2D<Traits>::writeGF(physics, gfMembranePotentialMD, SUBSAMPLING_POINTS,
            filePrefix + "memb_pot.dat", time, std::string(""), to_mV);

#if HAVE_HDF5
        if(doHDF5Output)
        {
          std::string datasetName("memb_pot");
          HDF5Tools<Traits>::writeGF(physics,gfMembranePotentialMD,SUBSAMPLING_POINTS,
              getHDF5FileName(),datasetName,to_mV);
        }
#endif

      }
    }


    void writeBoundaryOutput()
    {
      if(doFullGnuplotOutput || physics.getParams().boundary.get("writeBoundaryOutput",false))
      {
        std::vector<std::vector<typename Traits::DGF_POT::Traits::DomainType> > pos_vec;
        std::vector<std::vector<typename Traits::DGF_POT::Traits::RangeType> > sol_vec;
        std::vector<typename Ax1Output2D<Traits>::SolutionVectorInfo> info;
        Ax1Output2D<Traits>::getBoundarySolutionVector(physics,dgfPot,pos_vec,sol_vec,info);

        std::string prefix("");
        for(int boundary = 0; boundary < 4; boundary++)
        {
          std::string filename = filePrefix + "boundary_pot_";
          switch(boundary)
          {
            case BoundaryPositions::BOUNDARY_BOTTOM :
            {
              filename += "bottom";
              break;
            }
            case BoundaryPositions::BOUNDARY_TOP :
            {
              filename += "top";
              break;
            }
            case BoundaryPositions::BOUNDARY_LEFT :
            {
              filename += "left";
              break;
            }
            case BoundaryPositions::BOUNDARY_RIGHT :
            {
              filename += "right";
              break;
            }
            default :
              DUNE_THROW(Dune::Exception, "Unexpected boundary position identifier!");
          }
          filename += ".dat";

          debug_jochen << "  Writing boundary file " << filename << std::endl;

          std::vector<typename Traits::DGF_POT::Traits::DomainType>& pos = pos_vec[boundary];
          std::vector<typename Traits::DGF_POT::Traits::RangeType>& sol = sol_vec[boundary];
          //debug_jochen << "My vector size: " << sol.size() << std::endl;

          const int rank = gv.comm().rank();
          const int size = gv.comm().size();

          //debug_jochen << "MPI rank: " << rank << std::endl;
          //debug_jochen << "MPI size: " << size << std::endl;

          if(boundary == BoundaryPositions::BOUNDARY_BOTTOM
              || boundary == BoundaryPositions::BOUNDARY_TOP)
          {
            debug_jochen << "    Collecting data..." << std::endl;
            std::vector<int> all_dimensions_x = physics.nElements_AllProcessors(0); // copy, may not be const in MPI_Gatherv
            std::vector<int> displacements = physics.nOffset_AllProcessors(0); // copy, may not be const in MPI_Gatherv
            int total_size_x = std::accumulate(all_dimensions_x.begin(),all_dimensions_x.end(),0);
            total_size_x += 1;

            // Gather partial vectors on the root node
            std::vector<typename Traits::DGF_POT::Traits::DomainType> pos_all(rank == 0 ? total_size_x : 0);
            std::vector<typename Traits::DGF_POT::Traits::RangeType> sol_all(rank == 0 ? total_size_x : 0);

            // Parallel run: Communication is necessary!
            if(size > 1)
            {
              debug_jochen << "    Communicating..." << std::endl;
              // Let k_p be the number of elements on processor p
              // Take vertices (0,...,k_0) on root processor p=0
              // Take vertices (1,...,k_p) on each processor p >= 1
              // => get nonoverlapping vertex partition for boundary values
              if(rank == 0)
              {
                pos_all[0] = pos[0];
                sol_all[0] = sol[0];
              }

              // Old version with same number of elements on all processors
//              //debug_jochen << "Communicating vectors; my vector size: " << sol.size()
//              //    << " global vector size: " << sol_all.size() << std::endl;
//              gv.comm().gather(&pos[1],&pos_all[1],pos.size()-1,0);
//              gv.comm().gather(&sol[1],&sol_all[1],pos.size()-1,0);

              // New version which allows different numbers of elements in x-direction
              int status = MPI_Gatherv(&pos[1], info[boundary].dimensions[0]-1,
                  Dune::MPITraits<typename Traits::DGF_POT::Traits::DomainType>::getType(),
                  &pos_all[1], &all_dimensions_x[0], &displacements[0],
                  Dune::MPITraits<typename Traits::DGF_POT::Traits::DomainType>::getType(),
                  0, gv.comm());

              status = MPI_Gatherv(&sol[1], info[boundary].dimensions[0]-1,
                  Dune::MPITraits<typename Traits::DGF_POT::Traits::RangeType>::getType(),
                  &sol_all[1], &all_dimensions_x[0], &displacements[0],
                  Dune::MPITraits<typename Traits::DGF_POT::Traits::RangeType>::getType(),
                  0, gv.comm());

            } else {
              pos_all = pos;
              sol_all = sol;
            }

            debug_jochen << "    Writing output..." << std::endl;
            // Do actual output
            if(gv.comm().rank() == 0)
            {
              GnuplotTools2D<Traits>::gnuplotAddBlock(filename,pos_all,sol_all,prefix,infoStream.str());
            }
          } else {
            // Complete boundary (left or right) is present on rank 0 and rank size-1.
            // => just output them directly without communication
            if(boundary == BoundaryPositions::BOUNDARY_LEFT)
            {
              if(rank == 0)
              {
                assert(sol.size() > 0);
                GnuplotTools2D<Traits>::gnuplotAddBlock(filename,pos,sol,prefix,infoStream.str(),false);
              } else {
                assert(sol.size() == 0);
              }
            }
            if(boundary == BoundaryPositions::BOUNDARY_RIGHT)
            {
              if(rank == size-1)
              {
                assert(sol.size() > 0);
                GnuplotTools2D<Traits>::gnuplotAddBlock(filename,pos,sol,prefix,infoStream.str(),false);
              } else {
                assert(sol.size() == 0);
              }
            }
          }
        }
      }
    }

    void writeDiagnosticOutput()
    {
      // Diagnostic output
      std::vector<Real> vec = diagInfo.getVectorRepresentation();
      std::vector<int> reductionFuncs = diagInfo.getReductionFunctions();

      // Reduce diagnostic data from all processors on root node
      for(int i=0; i<vec.size(); i++)
      {
        switch(reductionFuncs[i])
        {
          case DiagnosticInfo::ReductionFunction::Function_Identity :
          {
            // Do nothing
            break;
          }
          case DiagnosticInfo::ReductionFunction::Function_Min :
          {
            vec[i] = gv.comm().min(vec[i]);
            break;
          }
          case DiagnosticInfo::ReductionFunction::Function_Max :
          {
            vec[i] = gv.comm().max(vec[i]);
            break;
          }
          case DiagnosticInfo::ReductionFunction::Function_Mean :
          {
            vec[i] = gv.comm().sum(vec[i]);
            vec[i] /= gv.comm().size();
            break;
          }
          case DiagnosticInfo::ReductionFunction::Function_Sum :
          {
            vec[i] = gv.comm().sum(vec[i]);
            break;
          }
          default:
            DUNE_THROW(Dune::Exception, "Could not find a reduction function for diagnostic data"
                << " entry #" << i << "!");
        }
      }
      // Print reduced diagnostic data on root node
      if(gv.comm().rank() == 0)
      {
        //Output::printVector(vec);
        GnuplotTools2D<Traits>::gnuplotAddLine(physics.getParams().getDiagnosticsFilename(),
          diagInfo.time, vec, filePrefix);
      }
    }


    std::string getHDF5FileName()
    {
      std::stringstream filename;
      filename << fn_hdf5.getName() << ".h5";
      return filename.str();
    }

    // Constructor parameters
    const GV& gv;
    const SubGV& elecGV;
    const SubGV& membGV;

    MultiGFS& multigfs;
    GFS_POT& gfsPot;

    Acme2CylSolutionVectors<Acme2CylTraits>& solutionVectors;
    U&       u;
    PHYSICS& physics;
    const int intorderPot;

    Dune::SubsamplingVTKWriter<GV>     vtkwriterGV;
    Dune::SubsamplingVTKWriter<SubGV>  vtkwriterSubGV;

    const std::string vtkBasenameDomain;
    const std::string vtkBasenameMembrane;
    const std::string checkpointBasename;
    const std::string hdf5Basename;
    Dune::PDELab::FilenameHelper      fn_vtk_GV;
    Dune::PDELab::FilenameHelper      fn_vtk_SubGV;
    Dune::PDELab::FilenameHelper      fn_checkpoints;
    Dune::PDELab::FilenameHelper      fn_hdf5;

    const std::string filePrefix;

    typename Traits::DGF_POT           dgfPot;
    typename Traits::DGF_POT_GRAD      dgfPotGrad;
    typename Traits::DGF_POT_GRAD_ELEC dgfPotGradElec;

    typename Traits::GF_MEMB_POT_MD  gfMembranePotentialMD;
    typename Traits::GF_MEMB_POT     gfMembranePotential;

    const bool vtkOutput;
    const bool doFullGnuplotOutput;
    const bool doHDF5Output;

    bool hdf5FileInitialized;

    // Functor for converting pot to mV
    Ax1OutputFunctors::to_millivolt<PHYSICS> to_mV;

    Dune::Timer timer;
    DiagnosticInfo diagInfo;

    int nTimeSteps;
    double time;
    std::stringstream infoStream;

    int nExistingOutputFiles;
    std::vector<std::string> gnuplotFiles;
    std::vector<std::string> oneTimeGnuplotFiles;

};



#endif /* DUNE_AX1_ACME2CYL_OUTPUT_HH */
