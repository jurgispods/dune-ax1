#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif


// Preprocessor defines
// YaspGrid=1, UG=2
#define USE_GRID 1
// Sequential=0, Parallel=1
#ifndef AX1_PARALLEL
#define USE_PARALLEL (0)
#else
#define USE_PARALLEL (1)
#define USE_OVERLAP (1) // Does not converge without overlap, why?
#define OVERLAP_SIZE (1)
#endif


#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <typeinfo>
#include <unistd.h>

// Needs to be included before debugstream.hh, otherwise the corresponding operator<< is missing!
#include <dune/common/array.hh>
#include <dune/common/debugstream.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_MULTIDOMAINGRID
#include <dune/grid/multidomaingrid.hh>
#endif

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/ax1_yaspgrid_loadbalancer.hh>
#include <dune/ax1/common/ax1_gridgenerator.hh>
#include <dune/ax1/common/ax1_gridvector.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>
#include <dune/ax1/acme2_cyl/configurations/all_configurations.hh>
#include <dune/ax1/common/ax1_parallelhelper.hh>

const bool SCALE_TO_A_M2 = false;

static const int pot = 3;

enum { Mori_Na = Na+pot+1, Mori_K = K+pot+1, Mori_Cl = Cl+pot+1};

template<typename Config>
void aggregate(Acme2CylParameters& params, Config& config, int flux_choice, std::string inputDir, std::string outputDir,
    bool useMori)
{
  std::vector<std::string> filenames;

  if(flux_choice < 2)
  {
    filenames.push_back("memb_flux_na.dat");
    filenames.push_back("memb_flux_k.dat");
    filenames.push_back("memb_flux_cl.dat");
  } else {
    filenames.push_back("dummy_Na.wurst");
    filenames.push_back("dummy_K.wurst");
    filenames.push_back("dummy_Cl.wurst");
  }
  if(flux_choice != 1)
  {
    filenames.push_back("memb_pot.dat");
  } else {
    filenames.push_back("dummy_cap.wurst");
  }
  if(useMori)
  {
    filenames.push_back("memb_flux_mori_na.dat");
    filenames.push_back("memb_flux_mori_k.dat");
    filenames.push_back("memb_flux_mori_cl.dat");
  }

  std::vector<std::ifstream*> file_in;
  for(int i=0; i<filenames.size(); i++)
  {    
    debug_info << "Creating input filestream for file " << (inputDir + "/" + filenames[i]) 
      << "..." << std::endl;
    std::ifstream* ifs = new std::ifstream(inputDir + "/" + filenames[i]);
    file_in.push_back(ifs);
  }

  int masterFile = -1;
  std::vector<bool> good(filenames.size());
  for(int i=0; i<file_in.size(); i++)
  {
    good[i] = file_in[i]->good();
  }
  if(flux_choice < 2)
  {
    if(! good[Na] || ! good[K])
    {      
      DUNE_THROW(Dune::Exception, "Na or K flux file not readable");
    }
  }
  if(flux_choice != 1)
  {
    if(! good[2])
      DUNE_THROW(Dune::Exception, "memb pot file not readable");
  }
  if(flux_choice == 2)
    masterFile = pot;
  else
    masterFile = Na;

  std::string line;

  std::vector<double> xVec;
  std::vector<double> capacity;
  double yMemb = -1;

  std::ifstream memb_groups_in(inputDir + "/memb_groups.dat");
  int nMembraneElements = 0;
  if(!memb_groups_in.good())
  {
    debug_warn << "File 'memb_groups.dat' could not be found. Will try to determine membrane element coordinates manually."
        << std::endl;

    yMemb = params.yMemb()[0] + params.dMemb();
    const std::vector<double>& xNodes = params.X();
    for(int i=1; i<xNodes.size(); i++)
    {
      xVec.push_back(xNodes[i-1] + 0.5 * (xNodes[i]-xNodes[i-1]));

      // Try to get membrane permittivity; this will fail if there are more than one defined!
      typename Acme2CylParameters::KeyVector subs = params.membrane.getSubKeys();

      double permittivity = 2.0;
      int count = 0;
      for(int k=0; k<subs.size(); k++)
      {
        if(params.membrane.sub(subs[k]).hasKey("permittivity"))
        {
          count++;
          permittivity = params.membrane.sub(subs[k]).get("permittivity", 2.0);
        }
      }
      if(count > 1)
        DUNE_THROW(Dune::Exception, "Found more than one permittivity definition in membrane structure!");

      capacity.push_back(permittivity * con_eps0 / (Config::LENGTH_SCALE * params.dMemb()));
      nMembraneElements++;
    }

  } else {

    std::getline(memb_groups_in, line);
    std::getline(memb_groups_in, line);

    while(std::getline(memb_groups_in, line) && !line.empty())
    {
      nMembraneElements++;
      double x, y, permittivity;
      int membIndex, groupIndex, processor;
      std::stringstream line_str(line);
      // Extract x, y coordinates and the rest
      line_str >> x >> y >> membIndex >> groupIndex >> permittivity >> processor;

      xVec.push_back(x);

      // Units [F/m^2]
      capacity.push_back(permittivity * con_eps0 / (Config::LENGTH_SCALE * params.dMemb()));

      if(yMemb < 0.0)
        yMemb = y;
    }
  }

  // Now that we have the corresponding coordinates, go through each of the files and read out fluxes
  debug_info << "Determined membrane coordinates, #elements = " << nMembraneElements << std::endl;
  debug_info << "xVec.size(): " << xVec.size() << std::endl;
  debug_info << "--------------------------------------------------------------------------" << std::endl;
  Output::printVector(xVec);
  debug_info << "--------------------------------------------------------------------------" << std::endl;

  std::vector<double> memb_pot(nMembraneElements, 0.0);
  std::vector<double> previous_memb_pot(nMembraneElements, 0.0);
  double time = -1;
  double previous_time = -1;
  std::ofstream file_out(outputDir + "/total_flux.dat", std::ios_base::trunc);

  std::ofstream file_out_na(outputDir + "/total_flux_na.dat", std::ios_base::trunc);
  std::ofstream file_out_k(outputDir + "/total_flux_k.dat", std::ios_base::trunc);
  std::ofstream file_out_cl(outputDir + "/total_flux_cl.dat", std::ios_base::trunc);

  std::ofstream file_out_memb_total(outputDir + "/memb_flux_total.dat", std::ios_base::trunc);
  file_out_memb_total << std::setprecision(12) << std::scientific;

  std::ofstream file_out_memb_total_na(outputDir + "/memb_flux_total_na.dat", std::ios_base::trunc);
  file_out_memb_total_na << std::setprecision(12) << std::scientific;
  std::ofstream file_out_memb_total_k(outputDir + "/memb_flux_total_k.dat", std::ios_base::trunc);
  file_out_memb_total_k << std::setprecision(12) << std::scientific;
  std::ofstream file_out_memb_total_cl(outputDir + "/memb_flux_total_cl.dat", std::ios_base::trunc);
  file_out_memb_total_cl << std::setprecision(12) << std::scientific;

  std::ofstream file_out_memb_ionic(outputDir + "/memb_flux_ionic.dat", std::ios_base::trunc);
  file_out_memb_ionic << std::setprecision(12) << std::scientific;
  std::ofstream file_out_memb_cap(outputDir + "/memb_flux_cap.dat", std::ios_base::trunc);
  file_out_memb_cap << std::setprecision(12) << std::scientific;

  // Time
  std::ofstream file_out_time(outputDir + "/timesteps.dat", std::ios_base::trunc);
  file_out_time << std::setprecision(12) << std::scientific;

  bool convertAllInputFiles = true;
  std::vector<std::ofstream*> orig_file_out;
  if(convertAllInputFiles)
  {
    for(int i=0; i<filenames.size(); i++)
    {
      std::ofstream* ofs = new std::ofstream(outputDir + "/converted_" + filenames[i]);
      orig_file_out.push_back(ofs);
    }
  }

  int nLine = 0;
  while(std::getline(*file_in[masterFile], line))
  {
    nLine++;
    std::vector<std::vector<double> > single_values(filenames.size(), std::vector<double>(nMembraneElements, 0.0));

    // Handle all good files
    for(int j=0; j<file_in.size(); j++)
    {
      // Omit dummy files
      if(! good[j]) continue;

      if(j != masterFile)
      {
        std::getline(*file_in[j], line);
      }
      std::stringstream line_str(line);
      line_str >> time;

      if(j == masterFile && nLine % 100 == 0)
      {
        debug_info << "Processing line " << nLine << ", time = " << time << "..." << std::endl;
      }


      double value;
      // Read exactly nMembraneElements values
      for(int i=0; i<nMembraneElements; i++)
      {
        line_str >> value;
        if(j != pot)
        {
          single_values[j][i] += value;
        } else {
          memb_pot[i] = value;
          // Only calculate capacitive flux from 2nd timestep on
          if(previous_time > -1)
          {
            // Units [V/µs]
            double dv_dt = 1e-3 * (memb_pot[i] - previous_memb_pot[i]) / (time - previous_time);
            //debug_jochen << "memb_pot[i] " << memb_pot[i] << ", previous_memb_pot[i] " << previous_memb_pot[i]
            //  << ", time " << time << ", previous_time " << previous_time << std::endl;
            //debug_jochen << "time " << time << ", x=" << xVec[i] << ", dv_dt = " << dv_dt << std::endl;

            // Units [F/m^2] * [V/µs] / [ LS * C/mol] = [TS/LS mol/(m^2 s)]
            single_values[j][i] += (capacity[i] * dv_dt / (Config::LENGTH_SCALE * con_e * con_mol));
          }
        }
      }
    }

    previous_memb_pot = memb_pot;
    previous_time = time;

    // Intermediate solution: Multiply by scaling factor to have units [A/m^2] for comparison
    if(SCALE_TO_A_M2)
    {
      for(int k=0; k<single_values.size(); k++)
      {
        for(int i=0; i<nMembraneElements; i++)
        {
          single_values[k][i] *= con_e * con_mol * (Config::LENGTH_SCALE / Config::TIME_SCALE);
        }
      }
    }

    // All fluxes have been aggregated, write to file
    {
      // Open new scope for make ios_base_all_saver work
      Dune::ios_base_all_saver putengeschnetzeltes(file_out);
      file_out << std::setprecision(12) << std::scientific;
      file_out << "# time: " << time << std::endl;

      file_out_na << std::setprecision(12) << std::scientific;
      file_out_na << "# time: " << time << std::endl;
      file_out_k << std::setprecision(12) << std::scientific;
      file_out_k << "# time: " << time << std::endl;
      file_out_cl << std::setprecision(12) << std::scientific;
      file_out_cl << "# time: " << time << std::endl;

      if(convertAllInputFiles)
      {
        for(int l=0; l<orig_file_out.size(); l++)
        {
          *orig_file_out[l] << std::setprecision(12) << std::scientific;
          *orig_file_out[l] << "# time: " << time << std::endl;
        }
      }

      file_out_memb_total << time;

      file_out_memb_total_na << time;
      file_out_memb_total_k << time;
      file_out_memb_total_cl << time;

      file_out_memb_ionic << time;
      file_out_memb_cap << time;

      file_out_time << nLine << " " << time;

      for(int i=0; i<nMembraneElements; i++)
      {
        double sum_ionic = single_values[Na][i] + single_values[K][i] + single_values[Cl][i];
        double sum_na = single_values[Na][i];
        double sum_k = single_values[K][i];
        double sum_cl = single_values[Cl][i];

        if(useMori)
        {
          sum_ionic += single_values[Mori_Na][i] + single_values[Mori_K][i] + single_values[Mori_Cl][i];
          sum_na += single_values[Mori_Na][i];
          sum_k += single_values[Mori_K][i];
          sum_cl += single_values[Mori_Cl][i];
        }

        double sum_all = sum_ionic;

        // In the Mori case, the capacitive flux is calculated explicitly and contained in the ionic flux
        if(! useMori)
        {
          sum_all += single_values[pot][i];
        }
        file_out << xVec[i] << " " << yMemb << " " << sum_all << std::endl;

        file_out_na << xVec[i] << " " << yMemb << " " << sum_na << std::endl;
        file_out_k << xVec[i] << " " << yMemb << " " << sum_k << std::endl;
        file_out_cl << xVec[i] << " " << yMemb << " " << sum_cl << std::endl;


        for(int l=0; l<orig_file_out.size(); l++)
        {
          *orig_file_out[l] << xVec[i] << " " << yMemb << " " << single_values[l][i] << std::endl;
        }

        file_out_memb_total << " " << sum_all;

        file_out_memb_total_na << " " << sum_na;
        file_out_memb_total_k << " " << sum_k;
        file_out_memb_total_cl << " " << sum_cl;

        file_out_memb_ionic << " " << sum_ionic;
        file_out_memb_cap << " " << single_values[pot][i];
      }
    }
    file_out << std::endl << std::endl;

    file_out_na << std::endl << std::endl;
    file_out_k << std::endl << std::endl;
    file_out_cl << std::endl << std::endl;

    for(int l=0; l<orig_file_out.size(); l++)
    {
      *orig_file_out[l] << std::endl << std::endl;
    }
    file_out_memb_total << std::endl;

    file_out_memb_total_na << std::endl;
    file_out_memb_total_k << std::endl;
    file_out_memb_total_cl << std::endl;

    file_out_memb_ionic << std::endl;
    file_out_memb_cap << std::endl;

    file_out_time << std::endl;
  }
  debug_info << "DONE, processed " << nLine << " lines." << std::endl;
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is aggregateFluxes. Weniger cremig." << std::endl;
    else
	  {
      if(helper.rank()==0)
        std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
	  }
    
    if (argc<2)
	  {
      if(helper.rank()==0)
      {
        debug_info << "usage: ./aggregateFluxes <flux_choice> [config-file] [outputDir]" << std::endl;
        debug_info << "  possible values for flux_choice:" << std::endl;
        debug_info << "    0: total trans-membrane flux" << std::endl;
        debug_info << "    1: Na+K (ionic) fluxes only" << std::endl;
        debug_info << "    2: capacitive flux only" << std::endl;
      }

      return 1;
	  }
    
    int flux_choice;
    int status = sscanf(argv[1],"%d",&flux_choice);

    if(!status)
    {
      if(helper.rank()==0)
      {
        debug_info << "usage: ./aggregateFluxes <flux_choice> [config-file] [output-dir]" << std::endl;
        debug_info << "  possible values for flux_choice:" << std::endl;
        debug_info << "    0: total trans-membrane flux" << std::endl;
        debug_info << "    1: Na+K (ionic) fluxes only" << std::endl;
        debug_info << "    2: capacitive flux only" << std::endl;
      }
      return 1;
    }
    
    debug_info << "Flux choice: " << flux_choice << std::endl;

    std::string configFileName;
    if (argc>=3)
    {
      configFileName = argv[2];
    } else {
      configFileName = "acme2_cyl_par_laplace.config";
    }
    debug_info << "Using config file " << configFileName << std::endl;

    Dune::Timer timer;

    // Read config file
    Acme2CylParameters params;
    Dune::ParameterTreeParser::readINITree(configFileName, params, true);
    params.init(configFileName, argc, argv);

    std::string outputDir = "dummy_dir";
    if (argc>=4)
    {
      outputDir = argv[3];
    } else {
      outputDir = params.general.get("outputDir", "dummy_dir");

      std::size_t pos = configFileName.rfind("/");
      if (pos != std::string::npos)
      {
        outputDir = configFileName.substr(0,pos) + "/..";
      }
    }
    if(access(outputDir.c_str(),0) != 0)
    {
      DUNE_THROW(Dune::IOError, "The specified output directory '" +
          outputDir + "' does not exist. Please create it first!");
    }
    debug_info << "Using outputDir '" << outputDir << "'" << std::endl;

    std::string inputDir = params.general.get("outputDir", "dummy_dir");
    // Output dir as specified in configfile does not exist; this might be due to a relative path
    // in a file that has been copied to the "/config" subdir of the output dir, so try the father dir!
    if(access(inputDir.c_str(),0) != 0)
    {
      std::size_t pos = configFileName.rfind("/");
      inputDir = configFileName.substr(0,pos) + "/..";
    }

    if(access(inputDir.c_str(),0) != 0)
    {
      DUNE_THROW(Dune::IOError, "The specified input directory '" +
          inputDir + "' could not be found!");
    }
    debug_info << "Using inputDir '" << inputDir << "'" << std::endl;


    // Generate grid coordinates (might need those in case memb_groups.dat is missing)
    // Use level 0 here hardcoded
    Ax1GridGenerator gridGenerator(params, 0);
    GridVector<double> x;
    CylinderGridVector<double> y(params.dX());
    std::vector<typename Ax1GridGenerator::MembraneGroupTuple> membGroups;
    gridGenerator.generateTensorGrid(x, y, membGroups);
    params.setMembraneGroups(membGroups);

    // Now copy x,y gridvectors in to a std::vector which can be passed to the parameter class
    std::vector<double> x_std(x);
    std::vector<double> y_std(y);
    params.setX(x_std);
    params.setY(y_std);


    bool useMori = params.useMori() || params.general.get("forceMoriFluxCalculation",false);

    debug_info << "!!! Using Mori fluxes: " << useMori << " !!!" << std::endl;

    bool configFound = false;
    std::string configName = params.getConfigName();

    if(configName == "default")
    {
      DefaultConfiguration<double> config(params);
      aggregate(params,config,flux_choice,inputDir,outputDir,useMori);
      configFound = true;
    }
    if(configName == "ES")
    {
      ESConfiguration<double> config;
      aggregate(params,config,flux_choice,inputDir,outputDir,useMori);
      configFound = true;
    }
    if(configName == "mori")
    {
      MoriConfiguration<double> config(params);
      aggregate(params,config,flux_choice,inputDir,outputDir,useMori);
      configFound = true;
    }
    // Other configurations disabled to speed up compilation
    /*
    if(configName == "step")
    {
      StepConfiguration<double> config;
      aggregate(params,config,flux_choice,inputDir,outputDir,useMori);
      configFound = true;
    }
    if(configName == "test_scales")
    {
      TestScalesConfiguration<double> config;
      aggregate(params,config,flux_choice,inputDir,outputDir,useMori);
      configFound = true;
    }
    */

    if(! configFound)
    {
      DUNE_THROW(Dune::Exception, "No configuration named '" << configName << "' could be found!");
    }

    double elapsed = timer.stop();
    debug_info << "Time elapsed: " << elapsed << " s" << std::endl;

  } catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    throw e; // Throw exception in order to get a full stacktrace for debugging
  }

}
