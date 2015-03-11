/*
 * simulation_data.hh
 *
 *  Created on: Feb 28, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_SIMULATION_DATA_HH
#define DUNE_AX1_SIMULATION_DATA_HH

template<typename T>
class Ax1SimulationData
{
  public:

    struct SolutionVector
    {
      std::string name;
      std::vector<T> data;
    };

    Ax1SimulationData()
    : timeStep(0),
      time(0),
      dt(0),
      mpi_rank(0),
      mpi_size(1),
      filename("")
    {
    }

    Ax1SimulationData(int nElems_, double time_, double dt_, std::string filename_)
    : nElements(nElems_),
      timeStep(0),
      time(time_),
      dt(dt_),
      mpi_rank(0),
      mpi_size(1),
      filename(filename_)
    {
    }

    void addVector(std::string name, const std::vector<T>& vec)
    {
      SolutionVector solVector;
      solVector.name = name;
      solVector.data = vec;
      solutionVectors.push_back(solVector);
    }

    void updateVector(std::string name, const std::vector<T>& vec)
    {
      getVector(name).data = vec;
    }

    void updateVector(int i, const std::vector<T>& vec)
    {
      getVector(i).data = vec;
    }

    int size() const
    {
      return solutionVectors.size();
    }

    const SolutionVector& getVector(int i) const
    {
      assert(i < solutionVectors.size());
      return solutionVectors[i];
    }

    const SolutionVector& getVector(std::string name) const
    {
      for(int i=0; i<solutionVectors.size(); i++)
      {
        if(solutionVectors[i].name == name)
          return solutionVectors[i];
      }
      DUNE_THROW(Dune::Exception, "Requested solution vector '" << name  << "' was not found!");
    }

    SolutionVector& getVector(int i)
    {
      assert(i < solutionVectors.size());
      return solutionVectors[i];
    }

    SolutionVector& getVector(std::string name)
    {
      for(int i=0; i<solutionVectors.size(); i++)
      {
        if(solutionVectors[i].name == name)
          return solutionVectors[i];
      }
      DUNE_THROW(Dune::Exception, "Requested solution vector '" << name  << "' was not found!");
    }

    bool hasVector(std::string name) const
    {
      for(int i=0; i<solutionVectors.size(); i++)
      {
        if(solutionVectors[i].name == name)
          return true;
      }
      return false;
    }

    std::string getFileName() const
    {
      return filename;
    }

    void setFileName(std::string filename_)
    {
      filename = filename_;
    }

    int getNElements() const
    {
      return nElements;
    }

    void setNElements(int nElems_)
    {
      nElements = nElems_;
    }

    int getTimeStep() const
    {
      return timeStep;
    }

    void setTimeStep(int timeStep_)
    {
      timeStep = timeStep_;
    }

    double getTime() const
    {
      return time;
    }

    void setTime(double time_)
    {
      time = time_;
    }

    double getDt() const
    {
      return dt;
    }

    void setDt(double dt_)
    {
      dt = dt_;
    }

    int getMpiRank() const
    {
      return mpi_rank;
    }

    void setMpiRank(int mpi_rank_)
    {
      mpi_rank = mpi_rank_;
    }

    int getMpiSize() const
    {
      return mpi_size;
    }

    void setMpiSize(int mpi_size_)
    {
      mpi_size = mpi_size_;
    }

    void addMetaInfo(std::string key, std::string value)
    {
      metaInfo[key] = value;
    }

    std::string getMetaInfo(std::string key) const
    {
      if(metaInfo.count(key) > 0)
        return metaInfo.at(key);
      else
        DUNE_THROW(Dune::Exception, "Could not find entry '" << key << "' in metaInfo map!");
    }

    void clear()
    {
      nElements = 0;
      time = 0;
      dt = 0;
      mpi_rank = 0;
      mpi_size = 1;
      filename = "";
      solutionVectors.clear();
    }


  private:
    int nElements;
    int timeStep;
    double time;
    double dt;
    int mpi_rank;
    int mpi_size;
    std::string filename;
    std::vector<SolutionVector> solutionVectors;
    std::map<std::string, std::string> metaInfo;

};


#endif /* DUNE_AX1_SIMULATION_DATA_HH */
