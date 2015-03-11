/*
 * acme1MD_solutionvectors.hh
 *
 *  Created on: Feb 28, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME1MD_SOLUTIONVECTORS_HH
#define DUNE_AX1_ACME1MD_SOLUTIONVECTORS_HH

//! \brief Simple structure to store references to all used solution vectors
//! in one place
template<typename Acme1MDTraits>
struct Acme1MDSolutionVectors
{
  typedef Acme1MDTraits Traits;

  Acme1MDSolutionVectors(typename Traits::U& uold_,
                         typename Traits::U& unew_,
                         typename Traits::U_CON& uoldCon_,
                         typename Traits::U_CON& unewCon_,
                         typename Traits::U_POT& uoldPot_,
                         typename Traits::U_POT& unewPot_)
    : uold(uold_),
      unew(unew_),
      uoldCon(uoldCon_),
      unewCon(unewCon_),
      uoldPot(uoldPot_),
      unewPot(unewPot_)
  {
  }

  void printAll()
  {
    printMultiDomainVectors();
    printConcentrationVectors();
    printPotentialVectors();
  }

  void printMultiDomainVectors()
  {
    debug_verb << "[uold:]" << std::endl;
    Output::printSingleCoefficientVector(uold, "uold");
    debug_verb << "[unew:]" << std::endl;
    Output::printSingleCoefficientVector(unew, "unew");
  }

  void printConcentrationVectors()
  {
    debug_verb << "[uoldCon:]" << std::endl;
    Output::printMultipleComponentCoefficientVector(uoldCon, NUMBER_OF_SPECIES);
    debug_verb << "[unewCon:]" << std::endl;
    Output::printMultipleComponentCoefficientVector(unewCon, NUMBER_OF_SPECIES);
  }

  void printPotentialVectors()
  {
    debug_verb << "[uoldPot:]" << std::endl;
    Output::printSingleCoefficientVector(uoldPot, "uoldPot");
    debug_verb << "[unewPot:]" << std::endl;
    Output::printSingleCoefficientVector(unewPot, "unewPot");
  }


  typename Traits::U& uold;
  typename Traits::U& unew;
  typename Traits::U_CON& uoldCon;
  typename Traits::U_CON& unewCon;
  typename Traits::U_POT& uoldPot;
  typename Traits::U_POT& unewPot;
};

#endif /* DUNE_AX1_ACME1MD_SOLUTIONVECTORS_HH */
