/*
 * acme2_solutionvectors.hh
 *
 *  Created on: Feb 28, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME2_SOLUTIONVECTORS_HH
#define DUNE_AX1_ACME2_SOLUTIONVECTORS_HH

//! \brief Simple structure to store references to all used solution vectors
//! in one place
template<typename Acme2Traits>
struct Acme2SolutionVectors
{
  typedef Acme2Traits Traits;

  Acme2SolutionVectors(typename Traits::U& uold_,
                         typename Traits::U& unew_ /*,
                         typename Traits::U_CON& uoldCon_,
                         typename Traits::U_CON& unewCon_,
                         typename Traits::U_POT& uoldPot_,
                         typename Traits::U_POT& unewPot_*/)
    : uold(uold_),
      unew(unew_)/*,
      uoldCon(uoldCon_),
      unewCon(unewCon_),
      uoldPot(uoldPot_),
      unewPot(unewPot_)*/
  {
  }

  void printAll() const
  {
    printMultiDomainVectors();
    //printConcentrationVectors();
    //printPotentialVectors();
  }

  void printMultiDomainVectors() const
  {
    debug_verb << "[uold:]" << std::endl;
    Output::printSingleCoefficientVector(uold, "uold");
    debug_verb << "[unew:]" << std::endl;
    Output::printSingleCoefficientVector(unew, "unew");
  }

  /*
  void printConcentrationVectors() const
  {
    debug_verb << "[uoldCon:]" << std::endl;
    Output::printMultipleComponentCoefficientVector(uoldCon, NUMBER_OF_SPECIES);
    debug_verb << "[unewCon:]" << std::endl;
    Output::printMultipleComponentCoefficientVector(unewCon, NUMBER_OF_SPECIES);
  }

  void printPotentialVectors() const
  {
    debug_verb << "[uoldPot:]" << std::endl;
    Output::printSingleCoefficientVector(uoldPot, "uoldPot");
    debug_verb << "[unewPot:]" << std::endl;
    Output::printSingleCoefficientVector(unewPot, "unewPot");
  }
  */

  void serialize(Ax1SimulationData<typename Acme2Traits::Real>& simulationData)
  {
    std::vector<typename Acme2Traits::Real> stdVec;

    uold.std_copy_to(stdVec);
    simulationData.addVector("uold", stdVec);

    unew.std_copy_to(stdVec);
    simulationData.addVector("unew", stdVec);

    /*
    uoldCon.std_copy_to(stdVec);
    simulationData.addVector("uoldCon", stdVec);

    unewCon.std_copy_to(stdVec);
    simulationData.addVector("unewCon", stdVec);

    uoldPot.std_copy_to(stdVec);
    simulationData.addVector("uoldPot", stdVec);

    unewPot.std_copy_to(stdVec);
    simulationData.addVector("unewPot", stdVec);
    */
  }

  void deserialize(const Ax1SimulationData<typename Acme2Traits::Real>& simulationData)
  {
    std::vector<typename Acme2Traits::Real> stdVec;

    stdVec = simulationData.getVector("uold").data;
    uold.std_copy_from(stdVec);

    stdVec = simulationData.getVector("unew").data;
    unew.std_copy_from(stdVec);

    /*
    stdVec = simulationData.getVector("uoldCon").data;
    uoldCon.std_copy_from(stdVec);

    stdVec = simulationData.getVector("unewCon").data;
    unewCon.std_copy_from(stdVec);

    stdVec = simulationData.getVector("uoldPot").data;
    uoldPot.std_copy_from(stdVec);

    stdVec = simulationData.getVector("unewPot").data;
    unewPot.std_copy_from(stdVec);
    */
  }




  typename Traits::U& uold;
  typename Traits::U& unew;
  /*typename Traits::U_CON& uoldCon;
  typename Traits::U_CON& unewCon;
  typename Traits::U_POT& uoldPot;
  typename Traits::U_POT& unewPot;*/
};

#endif /* DUNE_AX1_ACME2_SOLUTIONVECTORS_HH */
