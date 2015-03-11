/*
 * ax1_gridgeneration.hh
 *
 *  Created on: Jun 4, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_GRIDGENERATION_HH
#define DUNE_AX1_GRIDGENERATION_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/geometrygrid/grid.hh>
#if HAVE_DUNE_MULTIDOMAINGRID
#include <dune/grid/multidomaingrid.hh>
#endif

#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>
#include <dune/ax1/common/ax1_gridtools.hh>
#include <dune/ax1/common/ax1_gridvector.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/ax1_tensorgrid_transformation.hh>
#include <dune/ax1/common/ax1_yaspgrid_loadbalancer.hh>
#include <dune/ax1/common/ax1_groupmapper.hh>
#include <dune/ax1/common/element_subdomain_mapper.hh>

class Ax1GridGenerator
{
  public:
    typedef std::tuple<std::string, double, double> MembraneGroupTuple;

    // Grid stuff
    enum { dim=2 };
    typedef Ax1YaspPartition<dim,Dune::FieldVector<int,dim> > YP;
    typedef Ax1TensorGridTransformation<double> CoordFunction;

#if USE_GRID==1
    typedef Dune::GeometryGrid<BaseGrid,CoordFunction> HostGrid;
#elif USE_GRID==2
    typedef BaseGrid HostGrid;
#endif



    //typedef Dune::MultiDomainGrid<HostGrid, Dune::mdgrid::FewSubDomainsTraits<HostGrid::dimension,2> > Grid;
    typedef Dune::MultiDomainGrid<HostGrid,
        Dune::mdgrid::FewSubDomainsTraits<HostGrid::dimension,2,Dune::mdgrid::CellAndVertexCodims> > Grid;

    struct sort_tuple
    {
      bool operator() (const MembraneGroupTuple& a, const MembraneGroupTuple& b)
      {
        return (std::get<1>(a) < std::get<1>(b));
      }
    };


    Ax1GridGenerator(Acme2CylParameters& params_, const int level_)
    : params(params_),
      level(level_),
      qTransition(0.7)
    {}

    void checkForDoublettes(GridVector<double>& v)
    {
      for(int i=1; i<v.size(); i++)
      {
        if(std::abs(v[i]-v[i-1]) < 1e-6)
          DUNE_THROW(Dune::Exception, "Found duplicate value in grid vector (v[" << i << "] == v[" << (i-1) << "] == " << v[i] << ")!");
      }
    }

    void onesidedGeometricTowardsMembrane(CylinderGridVector<double>& y, double dy_membrane, int n_dy_membrane,
        double max_dy_equidistant, double yMemb)
    {
      double yStart = y.back();
      debug_jochen << "[onesidedGeometricTowardsMembrane] " << yStart << " -> " << yMemb << std::endl;

      // How many intervals do we need to smoothly transition between dy_membrane and max_dy_equidistant?
      double q = qTransition; // 'smoothness' factor q
      int n_transition = std::ceil(std::abs(std::log(dy_membrane / max_dy_equidistant) / std::log(q)));

      // range of the geometric series
      double geoRefineRange = std::max(max_dy_equidistant, dy_membrane)
        * (std::pow(q, n_transition+1) - 1) / (q - 1);
      debug_jochen << "geoRefineRange = " << geoRefineRange << std::endl;

      // The distance from the membrane to switch from linear to geometric refinement
      double geoRefineDistance = n_dy_membrane*dy_membrane + geoRefineRange;
      debug_jochen << "Switch distance from linear spacing to geometric refinement: "
          << geoRefineDistance << std::endl;

      double yrange_equidistant = yMemb - geoRefineDistance - yStart;
      // number of equidistant element that fit into remaining space
      int nEquidistant = (int) std::floor(yrange_equidistant / max_dy_equidistant);
      //double dy_linear = yrange_linear / nLinear;

      // Handle case when dy was chosen too large to obtain a fine grid resolution near the membrane
      if(yMemb - geoRefineDistance < yStart)
      {
        std::stringstream errorMsg;
        errorMsg << "Please check your grid parameters, cannot reach a sufficiently fine resolution of dy=";
        errorMsg << dy_membrane << " near the membrane with yMemb=" << yMemb;
        errorMsg << " when using a maximum dy=" << max_dy_equidistant << "!";
        DUNE_THROW(Dune::Exception, errorMsg.str());
      }
      debug_jochen << "Equidistant y-spacing till " << (yMemb - geoRefineDistance)
                   << ", then geometric refinement" << std::endl;

      debug_jochen << "yMemb = " << yMemb << std::endl;
      debug_jochen << "dMemb = " << params.dMemb() << std::endl;

      // Add equidistant elements
      debug_jochen << nEquidistant << " elements of height " << max_dy_equidistant << std::endl;
      if(nEquidistant > 0)
      {
        y.equidistant_n_h(nEquidistant, max_dy_equidistant);
      }

      // Now refine geometrically towards membrane; enforce hend to by guaranteed
      debug_jochen << "geometric transition till " << (yMemb - n_dy_membrane*dy_membrane) << std::endl;
      y.geometric_h0_hend_xend(max_dy_equidistant, dy_membrane, yMemb - n_dy_membrane*dy_membrane);

      // Add equidistant elements close to membrane
      debug_jochen << n_dy_membrane << " equidistant elements of height " << dy_membrane << " near membrane" << std::endl;
      y.equidistant_n_h(n_dy_membrane, dy_membrane);

      debug_jochen << "[onesidedGeometricTowardsMembrane] --> " << y.back() << std::endl;
    }

    void equidistantMembrane(CylinderGridVector<double>& y)
    {
      debug_jochen << "[equidistantMembrane] " << y.back() << " -> " << (y.back()+params.dMemb()) << std::endl;

      // Add one membrane element layer
      debug_jochen << "membrane thickness " << params.dMemb() << std::endl;

      // Only add a coordinate if the membrane has non-zero thickness
      if(params.dMemb() > 0)
      {
        if(! params.refineMembrane())
          y.equidistant_n_h(1, params.dMemb());
        else
          y.equidistant_n_h(params.nMembraneElements(), (params.dMemb() / params.nMembraneElements()));
      }

      debug_jochen << "[equidistantMembrane] --> " << y.back() << std::endl;
    }

    void twosidedGeometricBetweenMembranes(CylinderGridVector<double>& y, double dy_membrane, int n_dy_membrane,
        double max_dy_equidistant, double yMemb)
    {
      double yStart = y.back();
      debug_jochen << "[twosidedGeometricBetweenMembranes] " << yStart << " -> " << yMemb << std::endl;

      double interMembDistance = yMemb - yStart;

      // How many intervals do we need to smoothly transition between dy_membrane and max_dy_equidistant?
      double q = qTransition; // 'smoothness' factor q
      int n_transition = std::ceil(std::abs(std::log(dy_membrane / max_dy_equidistant) / std::log(q)));

      // range of the geometric series
      double geoRefineRangeSingle = max_dy_equidistant * (std::pow(q, n_transition+1) - 1) / (q - 1);
      debug_jochen << "geoRefineRange = " << geoRefineRangeSingle << std::endl;

      if(std::ceil(interMembDistance / dy_membrane) < 2*n_dy_membrane)
      {

        debug_jochen << "interMembDistance: " << interMembDistance << std::endl;
        debug_jochen << "max #equidistant elements:" << std::ceil(interMembDistance / dy_membrane) << std::endl;
        DUNE_THROW(Dune::Exception,
            "Minimum of " << (2*n_dy_membrane) << "equidistant elements do not fit between membranes!");
      }

      // Add equidistant elements close to membrane
      debug_jochen << n_dy_membrane << " equidistant elements of height " << dy_membrane << " near membrane" << std::endl;
      y.equidistant_n_h(n_dy_membrane, dy_membrane);

      double cellCenter = yStart + 0.5*interMembDistance;
      // This is the amount of spare space we have between the equidistant interval ends
      double geoRefineRangeTotal = interMembDistance - 2*n_dy_membrane*dy_membrane;
      debug_jochen << "Available range geoRefineRangeTotal = " << geoRefineRangeTotal << std::endl;

      if(std::floor((geoRefineRangeTotal - 2.0 * geoRefineRangeSingle) / max_dy_equidistant) < 2.0)
      {
        // No room for at least two equidistant elements in between => Check if there is at least some
        // space left for a geometric refinement towards the middle from both sides
        if(geoRefineRangeTotal > 0.0)
        {
          // Refine geometrically towards center and mirror this towards end
          debug_warn << "No space left for equidistant elements. Trying to do geometric refinement between "
              << y.back() << " and " << cellCenter << ", please check if result is acceptable!" << std::endl;

          // 1st half of the interval
          y.geometric_n_h0_xend(n_transition, dy_membrane, cellCenter);
          // 2nd half of the interval
          y.geometric_n_hend_xend(n_transition, dy_membrane, yStart + interMembDistance - n_dy_membrane*dy_membrane);
        } else {
          debug_warn << "No space left for equidistant elements nor for geometric refinement. "
              "No elements added between Debye layers!" << std::endl;
        }

      } else { // If room for at least 2 max_dy_equidistant cells => add equidistant elements in between

        double equidistantRangeTotal = interMembDistance - 2*n_dy_membrane*dy_membrane - 2.0*geoRefineRangeSingle;
        int nEquidistant = (int) std::floor(equidistantRangeTotal / max_dy_equidistant);
        equidistantRangeTotal = nEquidistant * max_dy_equidistant;

        double geoRefineEnd = yStart + 0.5*(interMembDistance-2*n_dy_membrane*dy_membrane-equidistantRangeTotal);

        // Geometrically coarsen
        //y.geometric_n_h0_q(n_transition, dy_membrane, 1./q);
        y.geometric_h0_hend_xend(dy_membrane, max_dy_equidistant, geoRefineEnd, 1);

        // Equidistant elements in the middle if the cell
        y.equidistant_n_h(nEquidistant, max_dy_equidistant);

        // Geometrically refine
        //y.geometric_n_h0_q(n_transition, max_dy_equidistant, q);
        y.geometric_h0_hend_xend(max_dy_equidistant, dy_membrane, yStart+interMembDistance-n_dy_membrane*dy_membrane);
      }


      // Add equidistant elements close to membrane
      debug_jochen << n_dy_membrane << " equidistant elements of height " << dy_membrane << " near membrane" << std::endl;
      y.equidistant_n_h(n_dy_membrane, dy_membrane);

      debug_jochen << "[twosidedGeometricBetweenMembranes] --> " << y.back() << std::endl;
    }

    void onesidedGeometricFromMembrane(CylinderGridVector<double>& y, double dy_membrane, int n_dy_membrane,
        double max_dy_equidistant, double yEnd)
    {
      double yStart = y.back();
      debug_jochen << "[onesidedGeometricFromMembrane] " << yStart << " -> " << yEnd << std::endl;

      // Store the volume of the first cell
      double min_extracellular_volume = (2*con_pi*params.dX()*yStart
                                        + con_pi*params.dX()*yStart)*dy_membrane;

      // Now mirror the whole procedure on the other side of the membrane!
      // Add equidistant elements close to membrane
      debug_jochen << n_dy_membrane << " equidistant elements of height " << dy_membrane << " near membrane" << std::endl;
      y.equidistant_n_h(n_dy_membrane, dy_membrane);



      // How many intervals do we need to smoothly transition between dy_membrane and max_dy_equidistant?
      double q = qTransition; // 'smoothness' factor q
      int n_transition = std::ceil(std::abs(std::log(dy_membrane / max_dy_equidistant) / std::log(q)));

      // range of the geometric series
      double geoRefineRange = max_dy_equidistant * (std::pow(q, n_transition+1) - 1) / (q - 1);
      debug_jochen << "geoRefineRange = " << geoRefineRange << std::endl;

      double geoRefineDistance = n_dy_membrane*dy_membrane + geoRefineRange;
      // ============= NEW => test! ========================================================================

      // Update ranges for the part above the membrane
      double yrange_equidistant = yEnd - yStart - geoRefineDistance;
      // number of equidistant element that fit into remaining space
      int nEquidistant = (int) std::floor(yrange_equidistant / max_dy_equidistant);

      bool useRadialCylinderRefinement = params.general.get("useRadialCylinderRefinement", false);
      if(not useRadialCylinderRefinement)
      {
        // Geometric coarsening
        debug_jochen << "geometric transition till " << (yEnd - (nEquidistant*max_dy_equidistant)) << std::endl;
        y.geometric_h0_hend_xend(dy_membrane, max_dy_equidistant, yEnd - (nEquidistant*max_dy_equidistant), 1);

        // Add linear elements until upper boundary is reached
        debug_jochen << nEquidistant << " elements of height " << max_dy_equidistant << std::endl;
        if(nEquidistant > 0)
        {
          y.equidistant_n_xend(nEquidistant, yEnd);
        }
      } else {
        // This is the maximum allowed volume ratio between two extracellular cells
        double max_volume_ratio = params.general.get("max_volume_ratio", 1000.);
        int max_nelements_y = params.general.get("max_radial_points", 1000.);

        /*
        double y_geoRefineEnd = y.back() + geoRefineRange;
        double max_volume = (2*con_pi*params.dX()*(y_geoRefineEnd - max_dy_equidistant)
                             + con_pi*params.dX()*(y_geoRefineEnd - max_dy_equidistant))*max_dy_equidistant;

        debug_jochen << "Maximum volume when using geometric refinement: " << max_volume << std::endl;
        debug_jochen << " => volume ratio: " << (max_volume / min_extracellular_volume) << std::endl;
        */


        bool transitionIntervalGeometricRefinement =
            params.general.get("transitionIntervalGeometricRefinement", false);

        double y_volumeEquilibrationStart =
            params.general.get("radialCylinderRefinementStart",y.back());
        y_volumeEquilibrationStart = std::max(y_volumeEquilibrationStart, y.back());
        double transitionInterval = y_volumeEquilibrationStart - y.back();

        // Temporary vector to store results of radial refinement
        CylinderGridVector<double> y_temp(params.dX());
        y_temp.push_back(y_volumeEquilibrationStart);

        // Calculate necessary parameters for radial refinement beforehand
        double remaining_volume = y_temp.remaining_volume(yEnd);
        int nelements_y = std::ceil(remaining_volume / (min_extracellular_volume * max_volume_ratio));
        debug_jochen << "Remaining volume: " << remaining_volume << std::endl;
        debug_jochen << "Would need " << nelements_y << " to reach a max volume ratio of "
            << max_volume_ratio << "." << std::endl;
        nelements_y = std::min(nelements_y, max_nelements_y);
        debug_jochen << "Using " << nelements_y << " => volume ratio "
            << (remaining_volume / (nelements_y * min_extracellular_volume)) << std::endl;
        y_temp.equivolume_n_rend(nelements_y, yEnd, false, 1, false);
        double h_radial_start = y_temp[1] - y_temp[0];
        debug_jochen << "First radial refinement interval: " << h_radial_start << std::endl;

        // Is there a transition interval at all or do we start radial refinement right away?
        if(transitionInterval > 0 && std::abs(transitionInterval > 1e-6))
        {
          // Transition interval: Do refinement in such a way that the volume increases
          // geometrically
          if(transitionIntervalGeometricRefinement)
          {
            //double q_cyl = params.general.get("q_cyl", 2.0);
            // This is the distance over which the volumes are (more or less) smoothly increased
            // until the volume equilibration strategy takes over
            double y_extracellular_geoRefineEnd = y.back() + transitionInterval;

            debug_jochen << "Doing geometric volume refinement from "
                << y.back() << " to " << y_extracellular_geoRefineEnd << std::endl;

            //y.geometric_volume_h_q_rend(dy_membrane, q_cyl, y_extracellular_geoRefineEnd, false, 0);
            y.geometric_h0_hend_xend(dy_membrane, h_radial_start, y_extracellular_geoRefineEnd, 1);
          } else {
            // No geometric refinement on transition interval, add equidistant elements instead
            int nElemsTransitionInterval = params.general.get("nElementsTransitionInterval", 10);
            y.equidistant_n_xend(nElemsTransitionInterval, y_volumeEquilibrationStart);
          }
        }

        remaining_volume = y.remaining_volume(yEnd);
        nelements_y = std::ceil(remaining_volume / (min_extracellular_volume * max_volume_ratio));
        debug_jochen << "Remaining volume: " << remaining_volume << std::endl;
        debug_jochen << "Would need " << nelements_y << " to reach a max volume ratio of "
            << max_volume_ratio << "." << std::endl;
        nelements_y = std::min(nelements_y, max_nelements_y);
        debug_jochen << "Using " << nelements_y << " => volume ratio "
            << (remaining_volume / (nelements_y * min_extracellular_volume)) << std::endl;

        y.equivolume_n_rend(nelements_y, yEnd, false, 1);
      }
      debug_jochen << "[onesidedGeometricFromMembrane] --> " << y.back() << std::endl;
    }


    void generateTensorGrid(GridVector<double>& x, CylinderGridVector<double>& y, std::vector<MembraneGroupTuple>& groups)
    {
      // ############## grid generation ##################
      int numMembranes = params.nMembranes();

      if(numMembranes > 3)
        DUNE_THROW(Dune::NotImplemented, "Maximum number of membranes is 3!");

      // TODO Replace by vector, one coordinate for each of the numMembranes membranes
      std::vector<double> yMembDefault(numMembranes, 500.);
      std::vector<double> yMemb = params.general.get("y_memb", yMembDefault);

      for(int i=0; i<yMemb.size(); i++)
      {
        debug_jochen << "yMemb[" << i << "] = " << yMemb[i] << std::endl;
      }

      double dMemb = params.dMemb();

      x.start(params.xMin());
      y.start(params.yMin());

      // ========================= y direction ==============================
      int n_yCytosol = -1;

      if(params.doRefineYDirectionGeometrically() && params.useMembrane())
      {
        // Maxmimum grid size in y-direction close to the membrane
        // TODO This should actually depend on the electrolyte Debye length!
        double dy_membrane = params.dYCellMin();
        int n_dy_membrane = params.general.get("n_dy_cell_min",10);

        // Default: Start with intracellular domain
        double max_dy_equidistant = params.dYCell();

        // Case two membranes: Start with extracellular domain
        if(numMembranes == 2)
        {
          // Maximum grid size in y-direction remote of the membrane
          double max_dy_equidistant = params.dY();
        }

        // ELEC -> MEMBRANE (only do this when there is a lower electrolyte!)
        if(params.yMemb()[0] > y.back()+1e-6)
        {
          // This config file flag can be utilized to force the grid generation to start with a Debye layer
          // resolution, even if there is no membrane located at the lower boundary; this is useful for
          // Laplace simulations, where the domain part below the lowest membrane was cut out
          if(params.general.get("refineTwoSidedFromBottom",false))
          {
            twosidedGeometricBetweenMembranes(y, dy_membrane, n_dy_membrane, max_dy_equidistant, yMemb[0]);
          } else {
            onesidedGeometricTowardsMembrane(y, dy_membrane, n_dy_membrane, max_dy_equidistant, yMemb[0]);
          }
          n_yCytosol = y.size();
        } else {
          n_yCytosol = 0;
        }

        // MEMBRANE (only do this when the other end of the membrane is not already a coordinate in y!)
        if(params.yMemb()[0]+params.dMemb() > y.back()+1e-6)
        {
          equidistantMembrane(y);
        }

        // In case of two membranes: Handle intracellular domain
        if(numMembranes == 2)
        {
          // MEMBRANE -> MEMBRANE (intracellular)
          twosidedGeometricBetweenMembranes(y, dy_membrane, n_dy_membrane, params.dYCell(), yMemb[1]);

          // MEMBRANE
          equidistantMembrane(y);
        }
        if(numMembranes == 3)
        {
          // MEMBRANE -> MEMBRANE (extracellular); use same dy than for intracellular domain for now!
          twosidedGeometricBetweenMembranes(y, dy_membrane, n_dy_membrane, params.dYCell(), yMemb[1]);

          // MEMBRANE
          equidistantMembrane(y);

          // MEMBRANE -> MEMBRANE (inracellular)
          twosidedGeometricBetweenMembranes(y, dy_membrane, n_dy_membrane, params.dYCell(), yMemb[2]);

          // MEMBRANE
          equidistantMembrane(y);
        }

        // TODO Case numMembranes == 3
        dy_membrane = params.dYMin();
        n_dy_membrane = params.general.get("n_dy_min",10);
        onesidedGeometricFromMembrane(y, dy_membrane, n_dy_membrane, params.dY(), params.yMax());

      } else {
        // Case: Not using (geometric refinement AND membrane)

        int n_y = std::pow(2,level);
        // Equidistant elements, number specified by refinement level
        if(params.useMembrane())
        {
          if(numMembranes > 1)
            DUNE_THROW(Dune::NotImplemented, "More than one membrane not implemented for non-geometric refinement!");

          if(params.doRefineYDirection())
          {
            n_y = std::ceil((params.yMemb()[0] - params.yMin()) / params.dYCell());
            y.equidistant_n_xend(n_y,params.yMemb()[0]);

            // Add membrane elements
            equidistantMembrane(y);

            n_y = std::ceil((params.yMax() - y.back()) / params.dY());
            y.equidistant_n_xend(n_y,params.yMax());
          } else {

            y.equidistant_n_xend(n_y,yMemb[0]);

            n_yCytosol = y.size();

            // Add membrane elements
            equidistantMembrane(y);

            // Equidistant elements, number specified by refinement level
            y.equidistant_n_xend(n_y,params.yMax());
          }
        } else {
          // A little awkward, but: Allow geometric refinement even for the case of no membrane, as one might
          // like to generate the same extracellular grid for comparison with other simulations!
          if(params.doRefineYDirectionGeometrically())
          {
            // Maxmimum grid size in y-direction close to the membrane
            // TODO This should actually depend on the electrolyte Debye length!
            const double dy_membrane = params.dYMin();
            int n_dy_membrane = params.general.get("n_dy_min",10);

            // Maximum grid size in y-direction remote of the membrane
            double max_dy_equidistant = params.dY();

            onesidedGeometricFromMembrane(y, dy_membrane, n_dy_membrane, max_dy_equidistant, params.yMax());

          } else {
            // When not using membrane: when 'refineYDirection' flag is set, use fixed size dy from config file
            if(params.doRefineYDirection())
            {
              n_y = std::ceil((params.yMax() - params.yMin()) / params.dY());
              y.equidistant_n_xend(n_y,params.yMax());
            } else {
              y.equidistant_n_xend(n_y,params.yMax());
            }
          }
        }
      }

      if(std::abs(y.back() - params.yMax()) > 1e-6)
        DUNE_THROW(Dune::Exception, "The maximum y value of " << y.back() << " does not match the grid parameter ymax=" << params.yMax()
            << ", probably the grid parameters from the config file (especially the value of dy=" << params.dY() << ") were not compatible "
            << " with the chosen refinement method!");

      params.setYOffsetMembrane(n_yCytosol);

      // print vector
      //Output::printVector(y);
      /*
      // print cell heights
      typename GridVector<BaseGrid>::const_iterator yit_last = y.begin();
      for(typename GridVector<BaseGrid>::const_iterator yit = y.begin(); yit != y.end(); ++yit)
      {
        if(yit != y.begin())
        {
          debug_jochen << "h_i = " << (*yit - *yit_last) << " --> " << *yit << std::endl;
        }
        yit_last = yit;
      }
      */

      // ========================= x direction ==============================
      setupXPartition(groups);

      int n_x = 1;
      if(params.doRefineXDirection())
      {
        // FIXME: This case does not yet work with multiple membrane groups!

        // TODO Specify different levels for x and y refinement?
        // TODO Change this at least according to the xMax/yMax ratio
        n_x = std::pow(2,level);

        double xrange = params.xMax() - params.xMin();
        double d_x = params.dX();
        n_x = std::ceil(xrange / d_x);

        // Fill x range with n_x equidistant elements
        debug_jochen << "=> n_x = " << n_x << std::endl;
        x.equidistant_n_xend(n_x,params.xMax());

      } else {

        // Now refine all intervals and put everything in a gridvector
        intervalVectorToGridVector(groups, x);
      }

      if(std::abs(x.back() - params.xMax()) > 1e-6)
        DUNE_THROW(Dune::Exception, "The maximum x value of " << x.back() << " does not match the grid parameter xmax=" << params.xMax()
            << ", probably the grid parameters from the config file were not compatible with the chosen refinement method!");
    }

    void globalRefine(GridVector<double>& x, GridVector<double>& y, int level)
    {
      debug_info << "Globally refining with level " << level << "..." << std::endl;
      for(int l=0; l<level; l++)
      {
        debug_info << "Globally refining x vector of size " << x.size() << std::endl;
        globalRefine(x);
        debug_info << "  new size: " << x.size() << std::endl;
        debug_info << "Globally refining y vector of size " << y.size() << std::endl;
        globalRefine(y,true);
        debug_info << "  new size: " << y.size() << std::endl;
      }
    }


    void makeMultidomainGrid(GridVector<double>& x, GridVector<double>& y,
        Dune::MPIHelper& helper,std::vector<MembraneGroupTuple>& membGroups,
        Ax1ElementGroupMapper& elemGroupMapper, ElementSubdomainMapper& elemSubdomainMapper)
    {
      Dune::Timer timer;

      // YaspGrid + GeometryGrid stuff
#if USE_GRID==1
      //GridType hostGrid(coords);

      const int dim = 2;
      Dune::FieldVector<double,dim> L;
      L[0] = params.xMax();
      L[1] = params.yMax();
      Dune::array<int,dim> N;
      N[0] = x.size()-1;
      N[1] = y.size()-1;
      std::bitset<dim> periodic;
#if USE_OVERLAP == 1
      int overlap = OVERLAP_SIZE;
#else
      int overlap = 0;
#endif

      // Set up custom load balancing
      int px = helper.size();
      int py = 1; // Do not partition y-direction!


      if (px*py == 0 || px*py !=helper.size())
      {
        DUNE_THROW(Dune::Exception, "px*py = 0 or != np");
      }
      if(px > N[0])
      {
        DUNE_THROW(Dune::Exception, "Cannot partition grid onto " << px << " processes, as there are only "
            << N[0] << " element stripes in y-direction!");
      }

      Dune::FieldVector<int,2> yasppartitions;
      yasppartitions[0] = px;
      yasppartitions[1] = py;
      YP* yp = new YP(yasppartitions);
      if( helper.rank() == 0 )
      {
        debug_info << "Partitioning of YASP: " << yasppartitions << std::endl;
      }
      Dune::shared_ptr<BaseGrid> yaspGridPtr(new BaseGrid (helper.getCommunicator(),L,N,periodic,overlap,yp));
      baseGrids.push_back(yaspGridPtr);

      BaseGrid& yaspGrid = *yaspGridPtr;

      //yaspGrid.globalRefine(level);

      Dune::shared_ptr<CoordFunction> coordFunctionPtr(new CoordFunction(params.xMax(),params.yMax(),x,y));
      coordFunctions.push_back(coordFunctionPtr);

      Dune::shared_ptr<HostGrid> hostGridPtr(new HostGrid(yaspGrid,*coordFunctionPtr));
      hostGrids.push_back(hostGridPtr);

      HostGrid& hostGrid = *hostGridPtr;

#elif USE_GRID==2 // UG case
      //BaseGrid::setDefaultHeapSize(2000);
      typedef BaseGrid HostGrid;

      //get GridFactory instance
      Dune::GridFactory<BaseGrid> factory;

      //get geometry type
      const Dune::GeometryType gt(Dune::GeometryType::cube,2);

      Ax1GridTools::prepareGrid(factory, x, y, gt);

      //make grid
      Dune::shared_ptr<BaseGrid> pBaseGrid(factory.createGrid());
      HostGrid& hostGrid = *pBaseGrid;

#else
#warning "Selected grid could not be found!"
#endif
      debug_info << "Grid creation finished, time elapsed: " << timer.elapsed() << "s" << std::endl;
      // ############## END grid generation ###############

      // ################### Test grid decomposition ###############################
      debug_jochen << "Leaf GV size before grid.loadBalance(): " << hostGrid.leafGridView().size(0) << std::endl;
      hostGrid.loadBalance();
      debug_jochen << "Leaf GV size after grid.loadBalance(): " << hostGrid.leafGridView().size(0) << std::endl;

      typename HostGrid::Partition<Dune::PartitionIteratorType::Interior_Partition>::LeafGridView interiorGV
       = hostGrid.leafGridView<Dune::PartitionIteratorType::Interior_Partition>();
      debug_jochen << "Leaf interior GV size after grid.loadBalance(): " << interiorGV.size(0) << std::endl;

      typedef HostGrid::LeafGridView BaseGV;
      BaseGV basegv = hostGrid.leafGridView();

      debug_jochen << "Processor #" << helper.rank() << " has " << basegv.size(0)
          << "elements." << std::endl;
      debug_jochen << "Processor #" << helper.rank() << " has " << basegv.overlapSize(0)
          << " overlap element layers" << std::endl;
      debug_jochen << "Processor #" << helper.rank() << " has " << basegv.ghostSize(0)
          << " ghost element layers" << std::endl;

      typedef BaseGV::Codim<0>::EntityPointer BElementPointer;
      typedef Dune::SingleCodimSingleGeomTypeMapper<BaseGV, 0> BElementMapper;

      // Custom mapper to map element indices to subdomain indices
      BElementMapper belementMapper(basegv);

  //      typedef BaseGV::Codim<0>::Partition<
  //          Dune::PartitionIteratorType::Interior_Partition>::Iterator IElementIterator;
  //      // print info about interior elements
  //      for(IElementIterator ep = basegv.begin<0,Dune::PartitionIteratorType::Interior_Partition>();
  //          ep != basegv.end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++ep)
  //      {
  //        debug_jochen << "Processor #" << helper.rank() << " has element #"
  //            << belementMapper.map(*ep) << " @"
  //            << ep->geometry().center() << " partition: "
  //            << ep->partitionType() << std::endl;
  //      }
  //
  //      typedef BaseGV::Codim<0>::Partition<
  //          Dune::PartitionIteratorType::Overlap_Partition>::Iterator OElementIterator;
  //      // print info about overlap elements
  //      for(OElementIterator ep = basegv.begin<0,Dune::PartitionIteratorType::Overlap_Partition>();
  //          ep != basegv.end<0,Dune::PartitionIteratorType::Overlap_Partition>(); ++ep)
  //      {
  //        debug_jochen << "Processor #" << helper.rank() << " has element #"
  //            << belementMapper.map(*ep) << " @"
  //            << ep->geometry().center() << " partition: "
  //            << ep->partitionType() << std::endl;
  //      }
  //
  //      typedef BaseGV::Codim<0>::Partition<
  //          Dune::PartitionIteratorType::Ghost_Partition>::Iterator GElementIterator;
  //      // print info about ghost elements
  //      for(GElementIterator ep = basegv.begin<0,Dune::PartitionIteratorType::Ghost_Partition>();
  //          ep != basegv.end<0,Dune::PartitionIteratorType::Ghost_Partition>(); ++ep)
  //      {
  //        debug_jochen << "Processor #" << helper.rank() << " has element #"
  //            << belementMapper.map(*ep) << " @"
  //            << ep->geometry().center() << " partition: "
  //            << ep->partitionType() << std::endl;
  //      }

      // ============== Multidomaingrid stuff ================================================================


      //Grid grid(hostGrid,false);
      Dune::shared_ptr<Grid> gridptr(new Grid(hostGrid,false));
      grids.push_back(gridptr);

      Grid& grid = *gridptr;

      typedef Grid::LeafGridView GV;
      GV gv = grid.leafGridView();

      typedef GV::Codim<0>::Iterator ElementIterator;
      typedef GV::Codim<0>::EntityPointer ElementPointer;
      typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

      // Custom mapper to map element indices to subdomain indices
      ElementMapper elementMapper(gv);

      // Grid parameters
      std::vector<double> yMemb = params.yMemb();
      double dMemb = params.dMemb();
      int numMembranes = params.nMembranes();
      double cellWidth = params.getCellWidth();

      if(params.useMembrane() && params.dMemb() == 0 && params.nMembraneElements() > 0)
        DUNE_THROW(Dune::Exception, "When membrane thickness is 0, the number of membrane element layers (nMembraneElements)"
                                 << "must also be 0!");

      double xNode = params.xNode();
      double nodeWidth = params.dNode();

      double xmin = params.xMax();
      double xmax = params.xMin();

      int nElementsThisProcess = 0;
      // Tag elements by type
      for(ElementIterator ep = gv.begin<0>(); ep != gv.end<0>(); ++ep)
      {
        // Only count interior elements
        if(ep->partitionType() == Dune::PartitionType::InteriorEntity)
        {
          nElementsThisProcess++;
          // Determine min and max x coordinate belonging exclusively to this processor
          xmin = std::min(ep->geometry().corner(0)[0], xmin);
          xmax = std::max(ep->geometry().corner(ep->geometry().corners()-1)[0],xmax);
        }

        int elemIndex = elementMapper.map(*ep);

        // default: intracellular
        int subdomainIndex = CYTOSOL;

        typename ElementIterator::Entity::Geometry::GlobalCoordinate center = ep->geometry().center();
        typename ElementIterator::Entity::Geometry::GlobalCoordinate diagonal
          = ep->geometry().corner(ep->geometry().corners()-1);
        diagonal -= ep->geometry().corner(0);

        double h_y = std::abs(diagonal[1]);

        //assert(std::abs(h_y - d) < 1e-8);

        // Assume there is a single layer of elements for each membrane in the following
        if(params.useMembrane())
        {
          // The check for a membrane element is generic: Check if it is within the bounds
          // of any membrane
          for(int i=0; i<yMemb.size(); i++)
          {
            // We are on the membrane{
            if(center[1] > yMemb[i] and center[1] < yMemb[i]+dMemb)
            {
              subdomainIndex = MEMBRANE;
              break;
            }
          }

          // For the other checks, this is not as easy: Do it on a case-by-case basis for now
          if(numMembranes == 1)
          {
            // There is only one membrane => the rest is outside of the cell
            if(center[1] > yMemb[0]+dMemb)
            {
              subdomainIndex = ES;
            }
          }
          // In case of a second membrane layer
          if(numMembranes == 2)
          {
            // Tag this cell as extracellular if it is outside the two membrane layers
            if(center[1] < yMemb[0] or center[1] > yMemb[1]+dMemb)
            {
              subdomainIndex = ES;
            }
          }
          // 3 membranes: Extracellular space is between membranes 1 and 2 and above 3
          if(numMembranes == 3)
          {
            // Tag this cell as extracellular if it is between the two membrane layers
            if(center[1] > yMemb[0]+dMemb and center[1] < yMemb[1])
            {
              subdomainIndex = ES;
            }
            // Above the last membrane => extracellular
            if(center[1] > yMemb[2]+dMemb)
            {
              subdomainIndex = ES;
            }
          }
        } else {
          // If no membrane is present, everything is extracellular space
          subdomainIndex = ES;
        }
        //debug_jochen << "Tagging element @ " << center[1] << " as " << subdomainIndex << std::endl;

        // Store subdomain of this element [M->G]
        elemSubdomainMapper.setSubdomain(elemIndex, subdomainIndex);

        if(subdomainIndex != MEMBRANE)
        {
          // not membrane: subdomainIndex == groupIndex
          elemGroupMapper.setGroup(elemIndex, subdomainIndex);
        } else {

          // Config file parameters 'start', 'width', 'stride' in every group define intervals where group is defined;
          // partition was set up in Ax1GridGenerator and written into vector membrGroups
          bool foundGroup = false;
          for(std::vector<typename Ax1GridGenerator::MembraneGroupTuple>::const_iterator it = membGroups.begin(); it != membGroups.end(); ++it)
          {
            double start = std::get<1>(*it);
            double end = std::get<2>(*it);
            if(center[0] > start and center[0] < end)
            {
              std::string groupName = std::get<0>(*it);
              const Dune::ParameterTree::KeyVector& groupNames = params.membrane.getSubKeys();
              int groupIndex = -1;
              for(int i=0; i<groupNames.size(); i++)
              {
                if(groupName == groupNames[i])
                  groupIndex = i;
              }
              if(groupIndex == -1)
                DUNE_THROW(Dune::Exception, "Group index for membrane group '" << groupName << "' could not be found!");

              // Determine element group index as group index in config file +2 (solution_in, solution_ex are always there)
              elemGroupMapper.setGroup(elemIndex, groupIndex+2);

              //debug_jochen << "Assigned element @" << center << " to group index " << elemGroupMapper.map(elemIndex) << std::endl;

              foundGroup = true;
              break;
            }
          }

          if(! foundGroup)
            DUNE_THROW(Dune::Exception, "Could not assign element @" << center << " to a membrane group!");
        }
      }
      // Store number of interior elements belonging to this process
      params.setNElementsThisProcess(nElementsThisProcess);

      // Use std::cout in order to print on every processor
      std::cout << "Processor " << helper.rank()  << " has x interval [" << xmin  << " " << xmax << "]" << std::endl;

      // TODO Check if processor boundary is sufficiently far away from membrane group transition; only check xmax on all processors
      double nearestGroupTransitionToXMax = 0.0;
      std::string groupNameProcessorBoundary("");

      for(std::vector<typename Ax1GridGenerator::MembraneGroupTuple>::const_iterator it = membGroups.begin(); it != membGroups.end(); ++it)
      {
        if(std::get<1>(*it) < xmax && std::get<2>(*it) > xmax)
          groupNameProcessorBoundary = std::get<0>(*it);

        double nextGroupTransition = std::get<2>(*it);

        if(std::abs(nextGroupTransition-xmax) < std::abs(nearestGroupTransitionToXMax-xmax))
          nearestGroupTransitionToXMax = nextGroupTransition;
      }

      // Hardcoded check for node of ranvier: Processor boundary should not be here!
      bool doProcessorBoundaryCheck = params.general.get("doProcessorBoundaryCheck", true);
      // Only check when requested in config file and if multiple membrane groups are available
      doProcessorBoundaryCheck = doProcessorBoundaryCheck && params.membrane.get("useMultipleGroups", false);
      if(groupNameProcessorBoundary == "node_of_ranvier" && doProcessorBoundaryCheck)
        DUNE_THROW(Dune::Exception, "Processor boundary " << xmax << " is located at membrane group 'node_of_ranvier', please consider"
            << " adjusting your grid parameters or number of processors!");

      double safetyDistance = params.dX();
      if(params.membrane.sub(groupNameProcessorBoundary).get("smoothTransition", false))
      {
        double node_width = params.membrane.sub("node_of_ranvier").get("width",1.e3);
        safetyDistance = params.membrane.sub(groupNameProcessorBoundary).get("dx_transition", node_width);
        safetyDistance *= params.membrane.sub(groupNameProcessorBoundary).get("n_transition", 10);
      }

      // For all inner ('real') processor boundaries: Check if it sufficiently far from membrane group transition!
      if(std::abs(nearestGroupTransitionToXMax-xmax) < safetyDistance && xmax != params.xMax() && doProcessorBoundaryCheck)
        DUNE_THROW(Dune::Exception, "Processor boundary " << xmax << " is located very close to the next membrane group transition ("
            << nearestGroupTransitionToXMax << "), please consider adjusting your grid parameters or number of processors!");

      typedef Grid::SubDomainGrid SubDomainGrid;
      SubDomainGrid& elecGrid = grid.subDomain(0);
      SubDomainGrid& membGrid = grid.subDomain(1);
      typedef Grid::ctype ctype;

      typedef SubDomainGrid::LeafGridView SDGV;
      SDGV elecGV = elecGrid.leafGridView();
      SDGV membGV = membGrid.leafGridView();

      grid.startSubDomainMarking();
      // Only add interior elements, multidomaingrid handles overlap automatically!
      typedef GV::Codim<0>::Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator ElementIterator_Interior;
      for (ElementIterator_Interior eit = gv.begin<0,Dune::PartitionIteratorType::Interior_Partition>();
          eit != gv.end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++eit)
      {
        int elemIndex = elementMapper.map(*eit);
        int subdomainIndex = elemSubdomainMapper.map(elemIndex);

        // Init subgrids
        switch(subdomainIndex)
        {
          case CYTOSOL:
            grid.addToSubDomain(0,*eit);
            break;
          case ES:
            grid.addToSubDomain(0,*eit);
            break;
          case MEMBRANE:
            grid.addToSubDomain(1,*eit);
            break;
          default:
            DUNE_THROW(Dune::Exception, "Element could not be assigned to a pre-defined subdomain!");
        }

      }
      grid.preUpdateSubDomains();
      grid.updateSubDomains();
      grid.postUpdateSubDomains();
      // ===============================================================================================
    }

    std::vector<std::shared_ptr<Grid> > getGrids()
    {
      return grids;
    }

  private:

    // globally refine once
    void globalRefine(GridVector<double>& x, bool excludeMembrane = false)
    {
      debug_info << "Globally refining..." << std::endl;
      std::vector<double> yMemb = params.yMemb();

      const int oldSize = x.size();
      for(int i=0; i<oldSize-1; i++)
      {
        bool isMembraneInterface = false;

        if(excludeMembrane && params.useMembrane())
        {
          for(int j=0; j<yMemb.size(); j++)
          {
            if(std::abs(x[i]-yMemb[j]) < 1e-6)
              isMembraneInterface = true;
          }
        }

        if(!isMembraneInterface)
        {
          x.push_back(0.5*(x[i] + x[i+1]));
        }
      }

      std::sort(x.begin(),x.end());
    }

    void intervalVectorToGridVector(std::vector<MembraneGroupTuple>& groups, GridVector<double>& x)
    {
      double default_dx = params.dX();

      std::vector<MembraneGroupTuple>::const_iterator it_last = groups.begin();
      std::vector<MembraneGroupTuple>::const_iterator it_next = groups.begin();
      for(std::vector<MembraneGroupTuple>::const_iterator it = groups.begin(); it != groups.end(); ++it)
      {
        double previousSize = x.size();

        std::string groupName = std::get<0>(*it);
        double start = std::get<1>(*it);
        double end = std::get<2>(*it);
        double myWidth = end - start;

        double remaining_width = myWidth;
        double dx = params.membrane.sub(groupName).get("dx", default_dx);

        // Flag to determine if a different mesh size should be used near transitions to other groups
        bool smoothTransition = params.membrane.sub(groupName).get("smoothTransition", false);
        bool doGeometricRefinement = params.membrane.sub(groupName).get("geometricRefinement", false);

        double geoRefineRange_left = 0.0;
        double geoRefineRange_right = 0.0;

        debug_jochen << "Processing interval [" << std::get<1>(*it) << " " << end << "] (" << groupName << ")" << std::endl;

        it_last = (it != groups.begin() ? it-1 : it);
        it_next = (it != groups.end()-1 ? it+1 : it);

        bool groupTransitionAtBeginning = (it != groups.begin() && std::get<0>(*it) != std::get<0>(*it_last));
        bool groupTransitionAtEnd = (it != groups.end() && std::get<0>(*it) != std::get<0>(*it_next));

        bool smoothTransitionAtBeginning = false;
        bool smoothTransitionAtEnd = false;

        double dxGroupTransition_left = 0.0;
        double dxGroupTransition_right = 0.0;
        int nTransition_left = 0;
        int nTransition_right = 0;
        // Check for membrane group transition at the beginning of this interval
        if(groupTransitionAtBeginning && smoothTransition)
        {
          double lastWidth = std::get<2>(*it_last) - std::get<1>(*it_last);
          debug_jochen << "Group transition @ beginning, position " << std::get<1>(*it) << std::endl;
          // Transition mesh width: Use minimum of interval lengths if not specified in config file
          dxGroupTransition_left = params.membrane.sub(groupName).get("dx_transition", std::min(myWidth, lastWidth));

          nTransition_left = params.membrane.sub(groupName).get("n_transition", 10);
          nTransition_left = std::min(nTransition_left, (int) std::floor(remaining_width/dxGroupTransition_left));

          debug_jochen << "dx_transition_left = " << dxGroupTransition_left << ", nTransition_left = " << nTransition_left << std::endl;

          if(nTransition_left > 0)
          {
            smoothTransitionAtBeginning = true;

            // Refine beginning of this interval
            double range = nTransition_left*dxGroupTransition_left;
            x.equidistant_n_h(nTransition_left, dxGroupTransition_left);

            // Decrease the remaining width of this interval
            remaining_width -= range;
          }
        }

        // Check for membrane group transition at the end of this interval
        if(groupTransitionAtEnd && smoothTransition)
        {
          double nextWidth = std::get<2>(*it_next) - std::get<1>(*it_next);
          debug_jochen << "Group transition @ end, position " << std::get<2>(*it) << std::endl;
          // Transition mesh width: Use minimum of interval lengths if not specified in config file
          dxGroupTransition_right = params.membrane.sub(groupName).get("dx_transition", std::min(myWidth, nextWidth));
          nTransition_right = params.membrane.sub(groupName).get("n_transition", 10);
          nTransition_right = std::min(nTransition_right, (int) std::floor(remaining_width/dxGroupTransition_right));

          debug_jochen << "dx_transition_right = " << dxGroupTransition_right << ", nTransition_right = " << nTransition_right << std::endl;

          if(nTransition_right > 0)
          {
            smoothTransitionAtEnd = true;

            // Refine beginning of this interval
            double range = nTransition_right*dxGroupTransition_right;

            // Postpone adding elements to vector until the remaining space has been processed

            // Decrease the remaining width of this interval
            remaining_width -= range;
          }
        }

        // Handle Geometric transition
        if(doGeometricRefinement && remaining_width < myWidth)
        {
          // How many intervals do we need to smoothly transition between dy_membrane and max_dy_linear?
          double q = qTransition; // 'smoothness' factor q

          if(smoothTransitionAtBeginning)
          {
            int n_transition_geometric = std::ceil(std::abs(std::log(dxGroupTransition_left / dx) / std::log(q)));
            debug_jochen << "n_transition_geometric left = " << n_transition_geometric << std::endl;
            geoRefineRange_left = dx * (std::pow(q, n_transition_geometric+1) - 1) / (q - 1);
            debug_jochen << "geoRefineRange_left = " << geoRefineRange_left << std::endl;
          }
          if(smoothTransitionAtEnd)
          {
            int n_transition_geometric = std::ceil(std::abs(std::log(dxGroupTransition_right / dx) / std::log(q)));
            debug_jochen << "n_transition_geometric right = " << n_transition_geometric << std::endl;
            geoRefineRange_right = dx * (std::pow(q, n_transition_geometric+1) - 1) / (q - 1);
            debug_jochen << "geoRefineRange_right = " << geoRefineRange_right << std::endl;
          }

          debug_jochen << "remaining width before geometric refinement: " << remaining_width << std::endl;

          // This is the remaining space to be filled with equidistant elements if we did the naive transition, but
          // this range will probably not fit an integer number of intervals of size dx!
          double remaining_width_temp = remaining_width - (geoRefineRange_left+geoRefineRange_right);

          // So determine the number of intervals of size dx we can actually fit:
          int n_x = 0;
          if(remaining_width_temp > dx)
          {
            n_x = (int) std::floor(remaining_width_temp / dx);
            remaining_width_temp -= (n_x * dx);
          }

          debug_jochen << "remainder = " << remaining_width_temp << std::endl;

          // And then adjust geometric transition range accordingly (remaining_width_temp has to be distributed onto
          // geoRefineRange_left and geoRefineRange_right):
          if(smoothTransitionAtBeginning && smoothTransitionAtEnd)
          {
            geoRefineRange_left += 0.5*remaining_width_temp;
            geoRefineRange_right += 0.5*remaining_width_temp;
          } else {
            if(smoothTransitionAtBeginning)
              geoRefineRange_left += remaining_width_temp;
            else
              geoRefineRange_right += remaining_width_temp;
          }

          debug_jochen << "  Adjusting geoRefineRange_left  to " << geoRefineRange_left << std::endl;
          debug_jochen << "  Adjusting geoRefineRange_right to " << geoRefineRange_right << std::endl;

          if(geoRefineRange_left + geoRefineRange_right > remaining_width)
              DUNE_THROW(Dune::Exception, "Geometric refinement was specified for membrane group '" << groupName
                << "', but the remaining width (" << remaining_width << ") is smaller than the needed range of the"
                << " geometric series (left: " << geoRefineRange_left << " / right: " << geoRefineRange_left
                << ") to smoothly transition between grid sizes dx (" << dx << ") and dxTransition (left: " << dxGroupTransition_left
                << "/ right: " << dxGroupTransition_right << "for this interval. Consider adjusting your grid parameters!");

          // Substract geoRefineRange twice, as it is used for both left and right interval transitions!
          remaining_width -= (geoRefineRange_left + geoRefineRange_right);

          if(geoRefineRange_left > 0.0)
          {
            debug_jochen << "Geometric refinement from " << x.back() << " to " << (x.back() + geoRefineRange_left) << std::endl;

            // Geometrically increase grid size towards dx, guarantee dx at end
            x.geometric_h0_hend_xend(dxGroupTransition_left, dx, x.back()+geoRefineRange_left, 1);
          }

          // Postpone transition to right until after equidistant intervals have been added for the remaining width
        }

        // Process middle of the interval (complete if no group transitions): refine if possible
        // TODO Does this conflict with geometricRefinement?
        int n_x = 1;
        double dx_equidistant = remaining_width;
        if(remaining_width >= dx)
        {
          // Determine the number of points we can use to refine this interval (round; TODO would ceil be better here?)
          n_x = (int) std::round(remaining_width / dx);
          dx_equidistant = remaining_width / n_x;
        }
        double equidistant_range = n_x * dx_equidistant;

        // Refine interval if it is larger than dx
        if(equidistant_range > 0)
        {
          debug_jochen << "Adding " << n_x << " subintervals of size "<< dx_equidistant << " from " << x.back() << " till "
              << (x.back() + equidistant_range) << std::endl;
          x.equidistant_n_xend(n_x, x.back() + equidistant_range);
        }

        // At the end, add group transition refinement to the right end of this interval
        if(groupTransitionAtEnd && nTransition_right > 0)
        {
          // Geometric transition towards right transition
          if(geoRefineRange_right > 0.0)
          {
            debug_jochen << "Geometric refinement from "
                << x.back() << " to " << (end - (nTransition_right*dxGroupTransition_right)) << std::endl;

            if(std::abs(geoRefineRange_right - (end - (nTransition_right*dxGroupTransition_right) - x.back())) > 1e-6)
              DUNE_THROW(Dune::Exception, "Something went wrong when doing geometric refinement to the right for"
                  << " this membrane group, geometricRefineRange_right=" << geoRefineRange_right
                  << ", actually used range is " << (end - (nTransition_right*dxGroupTransition_right) - x.back()) << "!");

            // Geometrically increase grid size towards dx, guarantee dx at beginning
            x.geometric_h0_hend_xend(dx, dxGroupTransition_right, end - (nTransition_right*dxGroupTransition_right), 0);
          }

          x.equidistant_n_xend(nTransition_right, end);
        }

        int size = x.size() - previousSize;
        debug_jochen << "--------------------------" << std::endl;
        debug_jochen << "Current interval, added " << size << " intervals:" << std::endl;
        Output::printVectorInterval(x.end()-size-1, x.end());
        debug_jochen << "--------------------------" << std::endl;
      }
    }


    //! Create nonoverlapping partition of x vector for membrane groups and check for consistency
    void setupXPartition(std::vector<MembraneGroupTuple>& groups)
    {
      sort_tuple sort_groups;

      // Read all desired membrane groups from config file
      const std::vector<std::string>& membGroups = params.membrane.getSubKeys();

      debug_jochen << "Setting up membrane group partition..." << std::endl;
      if(params.membrane.get("useMultipleGroups", false))
      {
        for(int i=0; i<membGroups.size(); ++i)
        {
          double start = params.membrane.sub(membGroups[i]).get("start",-1.);
          double width = params.membrane.sub(membGroups[i]).get("width",-1.);
          double stride = params.membrane.sub(membGroups[i]).get("stride",-1.);

          debug_jochen << "Processing membrane group '" << membGroups[i] << "' with start=" << start << ", width=" << width
              << ", stride=" << stride << "..." << std::endl;

          if(start > -1 && width > -1)
          {
            MembraneGroupTuple mGroup(membGroups[i], start, std::min(start+width, params.xMax()));
            groups.push_back(mGroup);
            debug_jochen << "Added group interval [" << start << " " << std::min(start+width, params.xMax()) << "]" << std::endl;

            if(stride > -1)
            {
              if(width > stride)
                DUNE_THROW(Dune::Exception, "Inconsistent parameters for membrane group '" << membGroups[i]
                   << "', width (" << width << ") > stride (" << stride << ")");

              double x_next = start + stride;
              while(x_next < params.xMax())
              {
                MembraneGroupTuple mGroup(membGroups[i], x_next, std::min(x_next+width,params.xMax()));
                groups.push_back(mGroup);

                debug_jochen << "Added group interval [" << x_next << " " << std::min(x_next+width,params.xMax()) << "]" << std::endl;

                x_next += stride;
              }

            }
            debug_jochen << "...done!" << std::endl;
          } else {
            debug_jochen << "...no position information available!" << std::endl;
          }

        }
        std::sort(groups.begin(), groups.end(), sort_groups);
      }

      // Insert default group at empty intervals
      std::string defaultGroup = params.membrane.get("defaultGroup", "axon");

      double last_end = 0.0;
      for(std::vector<MembraneGroupTuple>::iterator it = groups.begin(); it != groups.end(); ++it)
      {
        // We found a non-assigned interval
        if(std::get<1>(*it) > last_end)
        {
          MembraneGroupTuple dGroup(defaultGroup, last_end, std::get<1>(*it));
          debug_jochen << "Adding default group for interval [" << last_end << " " << std::get<1>(*it) << "]" << std::endl;
          it = groups.insert(it, dGroup);
          it++;
        }
        last_end = std::get<2>(*it);
      }
      // One more time for the remaining interval at the end of the vector
      if(last_end < params.xMax())
      {
        MembraneGroupTuple dGroup(defaultGroup, last_end, params.xMax());
        debug_jochen << "Adding default group for interval [" << last_end << " " << params.xMax() << "]" << std::endl;
        groups.push_back(dGroup);
      }


      // Check for consistency (is this a complete, nonoverlapping partition of the x axis?)
      checkPartition(groups);
    }

    void checkPartition(const std::vector<MembraneGroupTuple>& groups)
    {
      // Print partition
      debug_info << "Membrane group partition:" << std::endl;
      debug_info << "--------------------------------------------------" << std::endl;

      double last_end = 0.0;
      for(std::vector<MembraneGroupTuple>::const_iterator it = groups.begin(); it != groups.end(); ++it)
      {
        debug_info << "Group '" << std::get<0>(*it) << "', start=" << std::get<1>(*it) << ", end=" << std::get<2>(*it) << std::endl;

        double this_start = std::get<1>(*it);
        if(this_start != last_end)
          DUNE_THROW(Dune::Exception, "Inconsistent partition! Start of current interval (" << this_start
              << ") does not equal end of last interval (" << last_end << ")!");

        last_end = std::get<2>(*it);
      }
      debug_info << "--------------------------------------------------" << std::endl;
    }


    Acme2CylParameters& params;
    const int level;
    const double qTransition;

    std::vector<std::shared_ptr<CoordFunction> > coordFunctions;

    std::vector<std::shared_ptr<BaseGrid> > baseGrids;
    std::vector<std::shared_ptr<HostGrid> > hostGrids;
    std::vector<std::shared_ptr<Grid> > grids;
};

#endif /* DUNE_AX1_GRIDGENERATION_HH */
