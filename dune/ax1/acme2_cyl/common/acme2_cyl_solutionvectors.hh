/*
 * acme2_cyl_solutionvectors.hh
 *
 *  Created on: Feb 28, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2CYL_SOLUTIONVECTORS_HH
#define DUNE_AX1_ACME2CYL_SOLUTIONVECTORS_HH

//! \brief Simple structure to store references to all used solution vectors
//! in one place
template<typename Acme2CylTraits>
struct Acme2CylSolutionVectors
{
  typedef Acme2CylTraits Traits;

  Acme2CylSolutionVectors(typename Traits::U& uold_,
                         typename Traits::U& unew_)
    : uold(uold_),
      unew(unew_)
  {}

  void printAll() const
  {
    printMultiDomainVectors();
  }

  void printMultiDomainVectors() const
  {
    debug_verb << "[uold:]" << std::endl;
    Output::printSingleCoefficientVector(uold, "uold");
    debug_verb << "[unew:]" << std::endl;
    Output::printSingleCoefficientVector(unew, "unew");
  }

  void serialize(Ax1SimulationData<typename Acme2CylTraits::Real>& simulationData)
  {
    std::vector<typename Acme2CylTraits::Real> stdVec(uold.flatsize());

    //uold.std_copy_to(stdVec);
    std::size_t i = 0;
    for (typename Acme2CylTraits::U::iterator it = uold.begin(); it != uold.end(); ++it)
    {
       stdVec[i] = *it;
       i++;
    }
    simulationData.addVector("uold", stdVec);

    //unew.std_copy_to(stdVec);
    i = 0;
    for (typename Acme2CylTraits::U::iterator it = unew.begin(); it != unew.end(); ++it)
    {
       stdVec[i] = *it;
       i++;
    }
    simulationData.addVector("unew", stdVec);
  }

  void deserialize(const Ax1SimulationData<typename Acme2CylTraits::Real>& simulationData)
  {
    std::vector<typename Acme2CylTraits::Real> stdVec;

    stdVec = simulationData.getVector("uold").data;

    //uold.std_copy_from(stdVec);
    std::size_t i = 0;
    for (typename Acme2CylTraits::U::iterator it = uold.begin(); it != uold.end(); ++it)
    {
       *it = stdVec[i];
       i++;
    }

    if(i != stdVec.size())
      DUNE_THROW(Dune::Exception, "Size of loaded solution vector uold (" << stdVec.size()
          << ") does not match size of uold on the current grid (" << i << ")! Did you accidentally load the solution vector of "
          << "one processor's solution vector onto a different processor?");

    stdVec = simulationData.getVector("unew").data;
    //unew.std_copy_from(stdVec);
    i = 0;
    for (typename Acme2CylTraits::U::iterator it = unew.begin(); it != unew.end(); ++it)
    {
       *it = stdVec[i];
       i++;
    }

    if(i != stdVec.size())
      DUNE_THROW(Dune::Exception, "Size of loaded solution vector uold (" << stdVec.size()
        << ") does not match size of uold on the current grid (" << i << ")! Did you accidentally load the solution vector of "
        << "one processor's solution vector onto a different processor?");
  }

  typename Traits::U& uold;
  typename Traits::U& unew;
};

#endif /* DUNE_AX1_ACME2CYL_SOLUTIONVECTORS_HH */
