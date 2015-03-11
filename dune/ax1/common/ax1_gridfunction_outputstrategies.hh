/*
 * ax1_gridfunction_outputstrategies.hh
 *
 *  Created on: May 28, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_GRIDFUNCTION_OUTPUTSTRATEGIES_HH
#define DUNE_AX1_GRIDFUNCTION_OUTPUTSTRATEGIES_HH

namespace
{
  enum BoundaryPositions { BOUNDARY_BOTTOM = 0, BOUNDARY_LEFT = 1, BOUNDARY_RIGHT = 2, BOUNDARY_TOP = 3 };

  enum OutputPoints { OUTPUT_VERTICES_NONCONFORMING = 0, OUTPUT_VERTICES_CONFORMING = 1,
                      OUTPUT_CELL_CENTERS = 2, OUTPUT_MEMBRANE = 3 };

  // Default output strategy
  template<typename AcmeOutputTraits, typename GF>
  struct OutputStrategy
  {
#ifdef MULTIPLE_MEMBRANE_ELEMENTS
    enum { value = OutputPoints::OUTPUT_CELL_CENTERS };
    // TODO Remove this, this was only inserted for one special simulation run
    //enum { value = OUTPUT_VERTICES_NONCONFORMING };
#else
    enum { value = OutputPoints::OUTPUT_CELL_CENTERS };
    // TODO Remove this, this was only inserted for one special simulation run
    //enum { value = OUTPUT_VERTICES_NONCONFORMING };
#endif
  };

  // Specialization for the ion flux; ion flux is discontinuous over element borders and therefore
  // looks as if it was oscillatory when visualized: evaluate it on cell centers for a 'smoother'
  // output
  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_ION_FLUX>
  {
   enum { value = OutputPoints::OUTPUT_CELL_CENTERS };
  };

  // Specialization for element-wise constant geometry GFs such as volume/surface/...
  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_VOLUME>
  {
   enum { value = OutputPoints::OUTPUT_CELL_CENTERS };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_PLAIN2D_AREA>
  {
   enum { value = OutputPoints::OUTPUT_CELL_CENTERS };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_BASE_SURFACE>
  {
   enum { value = OutputPoints::OUTPUT_CELL_CENTERS };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_SIDE_SURFACE>
  {
   enum { value = OutputPoints::OUTPUT_CELL_CENTERS };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_SURFACE>
  {
   enum { value = OutputPoints::OUTPUT_CELL_CENTERS };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_PARTITION>
  {
   enum { value = OutputPoints::OUTPUT_CELL_CENTERS };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_PERMITTIVITY>
  {
   enum { value = OutputPoints::OUTPUT_CELL_CENTERS };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_CONDUCTIVITY>
  {
   enum { value = OutputPoints::OUTPUT_CELL_CENTERS };
  };

  /************************************* MEMBRANE STUFF **************************************/
  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_FLUX>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_FLUX_DIFF>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_FLUX_DRIFT>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_FLUX_LEAK>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_FLUX_VOLTAGE_GATED>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_FLUX_MORI>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_CURRENT>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_PARAM_MEMB_FLUX>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::BGF_PARAM_POT_NEUMANN>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_MEMB_POT_MD>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_MEMB_POT>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_MEMB_POT_MEMBGV>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_CHANNEL>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

  template<typename AcmeOutputTraits>
  struct OutputStrategy<AcmeOutputTraits, typename AcmeOutputTraits::GF_MEMB_GROUPS>
  {
   enum { value = OutputPoints::OUTPUT_MEMBRANE };
  };

}


#endif /* DUNE_AX1_GRIDFUNCTION_OUTPUTSTRATEGIES_HH */
