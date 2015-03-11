/*
 * ax1_simulationstate.hh
 *
 *  Created on: Feb 28, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_SIMULATIONSTATE_HH
#define DUNE_AX1_SIMULATIONSTATE_HH

#include <dune/ax1/common/ax1_simulationdata.hh>

template<typename T>
class Ax1SimulationState
{
  public:

    typedef typename Ax1SimulationData<T>::SolutionVector SolutionVector;

    void setData(Ax1SimulationData<T>& simulationData_)
    {
      simulationData = simulationData_;
    }

    Ax1SimulationData<T>& getData()
    {
      return simulationData;
    }

    void saveState()
    {
      std::ofstream out(simulationData.getFileName());
      if(! out.good())
        DUNE_THROW(Dune::Exception, "Simulation state file " << simulationData.getFileName()
            << " could not be created!");

      out << std::scientific << std::setprecision(16);

      out << "#elements " << simulationData.getNElements() << std::endl;
      out << "timestep " << simulationData.getTimeStep() << std::endl;
      out << "time " << simulationData.getTime() << std::endl;
      out << "dt " << simulationData.getDt() << std::endl;
      out << "mpi_rank " << simulationData.getMpiRank() << std::endl;
      out << "mpi_size " << simulationData.getMpiSize() << std::endl;

      out << "#vectors " << simulationData.size() << std::endl;
      for(int i=0; i<simulationData.size(); i++)
      {
        SolutionVector vec = simulationData.getVector(i);
        out << vec.name << " " << vec.data.size() << std::endl;
        for(int j=0; j<vec.data.size(); j++)
        {
          out << vec.data[j] << " ";
        }
        //out.write((const char*) &vec.data[0], vec.data.size() * sizeof(T));
        out << std::endl;
      }
    }

    template<typename GV>
    void loadState(const GV& gv, std::string base_filename, int np)
    {
      if(np < 1)
        DUNE_THROW(Dune::Exception, "This method can only be used when loading from mutliple state files!");

      std::vector<T> uold;
      std::vector<T> unew;
      std::vector<T> channels;

      std::vector<bool> testMap;

      const int overlapSize = 1;
      const int nMembraneElements = 1;
      int nx = 0, ny = 0;

      int countDOFs = 0;
      int countChannelStates = 0;

      std::vector<int> ndofs_x(np,0);

       debug_jochen << "Loading from " << np << " different state files, base filename: " << base_filename
          << std::endl;

      for(int i=0; i<np; i++)
      {
        std::string filename = base_filename;
        std::stringstream suffix;
        suffix << "_p" << i;
        Tools::replaceAll(filename, std::string("_p0"), suffix.str());

        debug_jochen << "Reading in state file '" << filename << "'..." << std::endl;

        // After this call, the loaded data resides in class membrane 'simulationData'
        loadState(filename);
        debug_jochen << "sizeof(simulationData): " << sizeof(simulationData) << std::endl;

        debug_jochen << "jo" << std::endl;

        nx = simulationData.getVector("x").data.size();
        ny = simulationData.getVector("y").data.size();
        int ne = simulationData.getNElements();

        int nDOFs = simulationData.getVector("uold").data.size();
        int nChannelStates = simulationData.getVector("channelStates").data.size();

        debug_jochen << "Loaded file with " << nDOFs << " DOFs and " << nChannelStates << " channel states"
            << std::endl;

        int nVertices = nDOFs / (NUMBER_OF_SPECIES+1);

        // Assume partitioning in x-direction
        int nVerticesX = nVertices / ny;

        int nVerticesX_Interior = (i==0 || i==np-1) ? nVerticesX - overlapSize : nVerticesX - 2*overlapSize;

        int nDOFs_Interior = nVerticesX_Interior * ny * (NUMBER_OF_SPECIES+1);

        int nStatesPerElement = nChannelStates / (nVerticesX-1);

        if(i==0)
        {
          // Set size of solution vectors
          uold.resize(nx*ny*(NUMBER_OF_SPECIES+1));
          unew.resize(nx*ny*(NUMBER_OF_SPECIES+1));
          testMap.resize(nx*ny*(NUMBER_OF_SPECIES+1),false);
          channels.resize((nx-1)*nMembraneElements*nStatesPerElement);
        }

        debug_jochen << "p" << i << ": " << nVerticesX << " vertices in x-direction ("
            << nVerticesX_Interior << " interior vertices) => " <<  nDOFs_Interior
            << " interior DOFs!" << std::endl;

        // Now pick out interior DOFs from simulationData and put it into new vector

        // offset left: overlapSize+1, since also one column of interior DOFs is shared by two neighboring
        // processes!
        int offset_x = (i==0) ? 0 : overlapSize+1;
        // Offset right
        int offset_x_end = (i==np-1) ? 0 : overlapSize;

        int global_offset_x = 0;
        for(int k=0; k<i; k++)
        {
          global_offset_x += ndofs_x[k];
        }
        ndofs_x[i] = nVerticesX-offset_x-offset_x_end;
        debug_jochen << "ndofs_x[" << i << "] = " << ndofs_x[i] << std::endl;


        // PNP states
        int countDOFsThis = 0;
        int index = -1;
        for(int iy=0; iy<ny; iy++)
        {
          for(int ix = 0; ix<ndofs_x[i]; ix++)
          {
            // Determine local index within this processor's loaded solution vector
            int col = offset_x + ix;
            index = (iy*nVerticesX + col) * (NUMBER_OF_SPECIES+1);

            // Determine global index within the large new solution vector containing data from all processors
            // Running index (local within this processor block)
            int running_index = ix + iy*ndofs_x[i];

            // Don't ask me what this does. It does the right thing, that's all I know. I'm going to bed now.
            int global_index = (global_offset_x + iy*(nx - ndofs_x[i]) + running_index) * (NUMBER_OF_SPECIES+1);

            //debug_jochen << "DOFs: Trying to assign block starting at index " << global_index
            //    << " with old data index " << index << std::endl;

            for(int k = 0; k<(NUMBER_OF_SPECIES+1); k++)
            {
              //debug_jochen << "  " << countDOFs << " - " << index << std::endl;

              if(countDOFs > uold.size())
                DUNE_THROW(Dune::Exception, "Index " << (countDOFs) << " exceeds new vector uold of size "
                    << uold.size() << "!");

              if(index > simulationData.getVector("uold").data.size())
                DUNE_THROW(Dune::Exception, "Index " << (index) << " exceeds old vector uold of size "
                    << simulationData.getVector("uold").data.size() << "!");

              if(testMap[global_index])
              {
                DUNE_THROW(Dune::Exception, "Global index " << global_index << " has already been assigned!");
              }
              testMap[global_index] = true;

              uold[global_index] = simulationData.getVector("uold").data[index];
              unew[global_index] = simulationData.getVector("unew").data[index];
              countDOFs++;
              countDOFsThis++;
              index++;
              global_index++;
            }

          }
        }
        debug_jochen << "p" << i << ": Read out " << countDOFsThis << " DOFs! (Total count: "
            << countDOFs << " of " << uold.size() << ")" << std::endl;

        int overlapDOFs = (offset_x_end+offset_x)*ny*(NUMBER_OF_SPECIES+1);
        if(countDOFsThis+overlapDOFs != simulationData.getVector("uold").data.size())
        {
          DUNE_THROW(Dune::Exception, "Did not read out the expected number of DOFs: countDOFsThis = "
              << countDOFsThis << ", uold.size() = " << simulationData.getVector("uold").data.size()
              << ", overlap DOFs = " << overlapDOFs);
        }
//        if(countDOFsThis != nDOFs_Interior)
//        {
//          DUNE_THROW(Dune::Exception, "Number of read interior DOFs (" << countDOFsThis
//              << ") is not equal to the number of interior DOFs on thie processor (" << nDOFs_Interior << (")!"));
//        }
        if(index+(offset_x_end*(NUMBER_OF_SPECIES+1)) != simulationData.getVector("uold").data.size())
        {
          DUNE_THROW(Dune::Exception, "Did not read the expected number of DOFs from old vector uold: Read up to index "
              << index << ", uold.size() = " << simulationData.getVector("uold").data.size() << ", offset end: "
              << offset_x_end);
        }

        // Channel states
        int offset_x_channels = (i==0) ? 0 : overlapSize;
        // Offset right
        int offset_x_end_channels = (i==np-1) ? 0 : overlapSize;

        // TODO: This method will not work for nMembraneElements > 1, since data needs to be reordered in the
        // same way DOFs get reordered above!
        assert(nMembraneElements == 1);
        int countChannelThis = 0;
        for(int iy=0; iy<nMembraneElements; iy++)
        {
          for(int ix = 0; ix<nVerticesX_Interior-1; ix++)
          {
            int col = offset_x_channels + ix;
            index = (iy*(nVerticesX-1) + col) * nStatesPerElement;

            //debug_jochen << "Channels: Trying to assign index " << countChannelStates
            //    << " with old data index " << index << std::endl;

            for(int k = 0; k<nStatesPerElement; k++)
            {
              //debug_jochen << "  " << countChannelStates << " - " << index << std::endl;

              channels[countChannelStates] = simulationData.getVector("channelStates").data[index];
              countChannelStates++;
              countChannelThis++;
              index++;
            }
          }
        }
        debug_jochen << "p" << i << ": Read out " << countChannelThis << " channel states!" << std::endl;

        int overlapChannelStates = (offset_x_end_channels+offset_x_channels)*nMembraneElements*(nStatesPerElement);
        if(countChannelThis+overlapChannelStates
            != simulationData.getVector("channelStates").data.size())
        {
          DUNE_THROW(Dune::Exception, "Did not read out the expected number of channel states: countChannelThis = "
              << countChannelThis << ", channelStates.size() = "
              << simulationData.getVector("channelStates").data.size()
              << ", overlapChannelStates = " << overlapChannelStates);
        }
        if(index+(offset_x_end_channels*(nStatesPerElement))
            != simulationData.getVector("channelStates").data.size())
        {
          DUNE_THROW(Dune::Exception, "Did not read the expected number of DOFs from old vector uold: Read up to index "
              << index << ", channelStates.size() = " << simulationData.getVector("channelStates").data.size() << ", offset end: "
              << offset_x_end_channels);
        }
      }
      debug_jochen << "TOTAL: Read out " << countDOFs << " DOFs and "
                  << countChannelStates << " channel states!" << std::endl;

      debug_jochen << "  Wanted: " << uold.size() << " DOFs and " << channels.size()
          << " channel states" << std::endl;

      assert(countDOFs == uold.size());
      assert(countChannelStates == channels.size());

      for(int i=0; i<testMap.size(); i++)
      {
        if(!testMap[i])
        {
          DUNE_THROW(Dune::Exception, "Global index " << i << " has not been assigned!");
        }
      }

      simulationData.setFileName(base_filename);
      simulationData.setNElements((nx-1)*(ny-1));
      simulationData.setMpiRank(0);
      simulationData.setMpiSize(1);

      simulationData.updateVector("uold",uold);
      simulationData.updateVector("unew",unew);
      simulationData.updateVector("channelStates",channels);

    }

    void loadState(std::string filename)
    {
      simulationData.clear();
      std::ifstream in(filename);
      if(! in.good())
        DUNE_THROW(Dune::Exception, "Simulation state file " << filename << " could not be found!");

      // Read meta information
      std::string dummy;
      int nElements;
      int timeStep = -1;
      double time, dt;
      int mpi_rank = 0;
      int mpi_size = 1;

      std::string line("");
      std::getline(in, line);
      double value = -1;
      while(line.find("#vectors") == std::string::npos && in.good())
      {
        std::stringstream line_str(line);
        line_str >> dummy >> value;
        if(dummy == "#elements")
          nElements = (int) value;
        if(dummy == "timestep")
          timeStep = value;
        if(dummy == "time")
          time = value;
        if(dummy == "dt")
          dt = value;
        if(dummy == "mpi_rank")
          mpi_rank = (int) value;
        if(dummy == "mpi_size")
          mpi_size = (int) value;

        std::getline(in, line);
      }
      if(! in.good())
        DUNE_THROW(Dune::Exception, "Bad file format, expected line beginning with '#vectors' at some point!");
      //debug_verb << "time = " << time << ", dt = " << dt << std::endl;

      simulationData.setFileName(filename);
      simulationData.setNElements(nElements);
      simulationData.setTimeStep(timeStep);
      simulationData.setTime(time);
      simulationData.setDt(dt);
      simulationData.setMpiRank(mpi_rank);
      simulationData.setMpiSize(mpi_size);

      int nVectors;
      std::stringstream line_str(line);
      line_str >> dummy >> nVectors;

      //debug_verb << "nVectors = " << nVectors << std::endl;

      for(int i=0; i<nVectors; i++)
      {
        std::string name;
        int size;
        in >> name >> size;

        std::vector<T> vec;
        vec.resize(size);
        //in.read((char*)&vec.data[0], size * sizeof(T));

        //debug_verb << "reading " << name << "[" << size << "]" << std::endl;
        for(int j=0; j<size; j++)
        {
          in >> vec[j];
          //debug_verb << "read " << vec[j] << std::endl;
        }
        simulationData.addVector(name, vec);
      }
    }

  private:
    Ax1SimulationData<T> simulationData;

};

#endif /* DUNE_AX1_SIMULATIONSTATE_HH */
