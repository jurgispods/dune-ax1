#ifndef DUNE_AX1_ACME0_BOUNDARY_HH
#define DUNE_AX1_ACME0_BOUNDARY_HH

#include<dune/common/fvector.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/constraintsparameters.hh>

#include <dune/ax1/common/constants.hh>

// boundary condtion types #########################################################################

// concentration

class BCTypeSingleCon
: public Dune::PDELab::DirichletConstraintsParameters            /*@\label{bcp:base}@*/
{
  double time;
public:
  template<typename I>
  bool isDirichlet(
                   const I & intersection,   /*@\label{bcp:name}@*/
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
    return false;  // all Neumann b.c.
  }
  
  // boundary flux decider
  template<typename I>
  bool nonZeroBoundaryFlux(
                           const I & intersection,   /*@\label{bcp:name}@*/
                           const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord,
                           const unsigned int ionSpecies
                           ) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension> xg = intersection.geometry().global(coord);
    
    switch(ionSpecies)
    {
        // sodium
      case Na:
        if      ( xg[0] < 0.0 ) return false;
        else if ( xg[0] > 0.0 ) return false;
        break;
        
        // potassium
      case K:
        if      ( xg[0] < 0.0 ) return false;
        else if ( xg[0] > 0.0 ) return false;
        break;
        
        // chloride
      case Cl:
        if      ( xg[0] < 0.0 ) return false;
        else if ( xg[0] > 0.0 ) return false;
        break;

      // default value (should never be used)
      default:
        return false;
    }
  }
  
  // set time for subsequent evaluation
  void setTime (double t) { time = t; }
};

// potential

class BCTypePot
: public Dune::PDELab::DirichletConstraintsParameters            /*@\label{bcp:base}@*/
{
  double time;
public:
  template<typename I>
  bool isDirichlet(const I & intersection, const Dune::FieldVector<typename I::ctype, I::dimension-1>& coord) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension> xg = intersection.geometry().global(coord);

    if( xg[0] > 0.0 )
      return true; // Dirichlet
    else
      return false;   // Neumann
    
  }
  
  // set time for subsequent evaluation
  void setTime (double t) { time = t; }
};

#endif /* DUNE_AX1_ACME0_BOUNDARY_HH */
