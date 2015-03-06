/*
 * acme1_solutionvectors.hh
 *
 *  Created on: Feb 28, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME1_SOLUTIONVECTORS_HH
#define DUNE_AX1_ACME1_SOLUTIONVECTORS_HH

//! \brief Simple structure to store references to all used solution vectors
//! in one place
template<typename Acme1Traits>
struct Acme1SolutionVectors
{
  typedef Acme1Traits Traits;

  Acme1SolutionVectors(typename Traits::U_CON& uCon_,
                       typename Traits::U_POT& uPot_,
                       typename Traits::USUB_CON& uoldCon_Inside_,
                       typename Traits::USUB_CON& unewCon_Inside_,
                       typename Traits::USUB_POT& uPot_Inside_,
                       typename Traits::USUB_CON& uoldCon_Outside_,
                       typename Traits::USUB_CON& unewCon_Outside_,
                       typename Traits::USUB_POT& uPot_Outside_)
    : uCon(uCon_),
      uPot(uPot_),
      uoldCon_Inside(uoldCon_Inside_),
      unewCon_Inside(unewCon_Inside_),
      uPot_Inside(uPot_Inside_),
      uoldCon_Outside(uoldCon_Outside_),
      unewCon_Outside(unewCon_Outside_),
      uPot_Outside(uPot_Outside_)
  {
  }

  void printAll()
  {
    printHostGridVectors();
    printSubGridVectors();
  }

  void printHostGridVectors()
  {
    debug_verb << "[uCon:]" << std::endl;
    Output::printMultipleComponentCoefficientVector(uCon, NUMBER_OF_SPECIES);
    debug_verb << "[uPot:]" << std::endl;
    Output::printSingleCoefficientVector(uPot, "uPot");
  }

  void printSubGridVectors()
  {
    debug_verb << "[uoldCon_Inside:]" << std::endl;
    Output::printMultipleComponentCoefficientVector(uoldCon_Inside, NUMBER_OF_SPECIES);
    debug_verb << "[uoldCon_Outside:]" << std::endl;
    Output::printMultipleComponentCoefficientVector(uoldCon_Outside, NUMBER_OF_SPECIES);

    debug_verb << "[unewCon_Inside:]" << std::endl;
    Output::printMultipleComponentCoefficientVector(unewCon_Inside, NUMBER_OF_SPECIES);
    debug_verb << "[unewCon_Outside:]" << std::endl;
    Output::printMultipleComponentCoefficientVector(unewCon_Outside, NUMBER_OF_SPECIES);

    debug_verb << "[uPot_Inside:]" << std::endl;
    Output::printSingleCoefficientVector(uPot_Inside, "uPot_Inside");
    debug_verb << "[uPot_Outside:]" << std::endl;
    Output::printSingleCoefficientVector(uPot_Outside, "uPot_Outside");
  }


  // Host grid vectors
  typename Traits::U_CON& uCon;
  typename Traits::U_POT& uPot;

  // Subgrid vectors
  typename Traits::USUB_CON& uoldCon_Inside;
  typename Traits::USUB_CON& unewCon_Inside;
  typename Traits::USUB_POT& uPot_Inside;

  typename Traits::USUB_CON& uoldCon_Outside;
  typename Traits::USUB_CON& unewCon_Outside;
  typename Traits::USUB_POT& uPot_Outside;
};

#endif /* DUNE_AX1_ACME1_SOLUTIONVECTORS_HH */
