#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

// Preprocessor defines
// YaspGrid=1, UG=2
#define USE_GRID 2
// Sequential=0, Parallel=1
#define USE_PARALLEL 0

// Use LocalBasisCache in local operator
#define USE_CACHE 1
//#define DUNE_ISTL_WITH_CHECKING

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <typeinfo>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/dgfparser/dgfwriter.hh>
#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG 
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
#if HAVE_DUNE_MULTIDOMAINGRID
#include <dune/grid/multidomaingrid/multidomaingrid.hh>
#endif

// Include constants.hh here in order to make the DOF-ordering-hack in dune-multidomain work
// TODO Use PDELab native ordering once it was ported to the new branch
#include <dune/ax1/common/constants.hh>

#if HAVE_DUNE_MULTIDOMAIN
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinggridfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinglocalfunctionspace.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#endif

#include <dune/ax1/common/ax1_gridtools.hh>
#include <dune/ax1/common/ax1_gridvector.hh>
#include <dune/ax1/common/ax1_subgridview.hh>
#include <dune/ax1/common/ax1_subgrid_hostelement_mapper.hh>
#include <dune/ax1/common/gnuplot_tools.hh>

#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_factory.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_geometrywrapper.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_physics.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_setup.hh>
#include <dune/ax1/channels/channel_builder.hh>

#include <dune/ax1/common/ax1_lfs_tools.hh>



//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is acme2_cyl. Extrem cremig." << std::endl;
    else
    {
      if(helper.rank()==0)
        std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
    }
    
    if (argc<4)
    {
      if(helper.rank()==0)
        std::cout << "usage: ./acme2_cyl <level> <dtstart> <tend> [config-file]" << std::endl;
      return 1;
    }
    
    int level;
    sscanf(argv[1],"%d",&level);
    
    double dtstart;
    sscanf(argv[2],"%lg",&dtstart);
    
    double tend;
    sscanf(argv[3],"%lg",&tend);
    
    std::string configFileName;
    if (argc>=5)
    {
      configFileName = argv[4];
    } else {
      if(USE_CYLINDER_COORDINATES)
        configFileName = "acme2_cyl.config";
      else
        configFileName = "acme2_2d.config";
    }
    debug_verb << "Using config file " << configFileName << std::endl;

    Dune::Timer timer;

    // sequential version
    if (1 && helper.size()==1)
    {
      // Read config file
      Acme2CylParameters params;
      Dune::ParameterTreeParser::readINITree(configFileName, params, true);
      params.init(configFileName, argc, argv);
      
      //assert(params.useMembrane() == true);

      // ############## grid generation ##################
      double yMemb = params.yMemb();
      double dMemb = params.dMemb();
      double cellWidth = params.getCellWidth();

#if USE_GRID==1
      typedef BaseGrid GridType;
      //GridType baseGrid(coords);

      const int dim = 2;
      Dune::FieldVector<double,2*1> L;
      L[0] = params.xMax();
      L[1] = params.yMax();
      Dune::FieldVector<int,dim> N (1);
      Dune::FieldVector<bool,dim> periodic (false);
      int overlap = 0;
      GridType baseGrid (L,N,periodic,overlap);
      baseGrid.globalRefine(level);
#elif USE_GRID==2 // UG case
      typedef BaseGrid GridType;

      //BaseGrid::setDefaultHeapSize(2000);

      //get GridFactory instance
      Dune::GridFactory<BaseGrid> factory;

      int numMembranes = params.nMembranes();

      GridVector<BaseGrid> x;
      CylinderGridVector<BaseGrid> y(params.dX());
      x.start(params.xMin());
      y.start(params.yMin());

      // ========================= y direction ==============================
      // Implicitly starts at x = 0.0
      //y.start(0.0);
      double yMembEnd = params.yMemb() + params.dMemb();

      if(params.doRefineYDirectionGeometrically() && params.useMembrane())
      {
        // Maxmimum grid size in y-direction close to the membrane
        // TODO This should actually depend on the electrolyte Debye length!
        const double dy_membrane = params.dYMin();

        // Maximum grid size in y-direction remote of the membrane
        double max_dy_equidistant = params.dY();

        // Maximum grid size in y-direction within the cell
        if(params.nMembranes() == 1)
        {
          max_dy_equidistant = params.dYCell();
        }

        // How many intervals do we need to smoothly transition between dy_membrane and max_dy_linear?
        double q = 0.7; // 'smoothness' factor q
        int n_transition = std::ceil(std::abs(std::log(dy_membrane / max_dy_equidistant) / std::log(q)));

        // range of the geometric series
        double geoRefineRange = max_dy_equidistant * (std::pow(q, n_transition+1) - 1) / (q - 1);
        debug_jochen << "geoRefineRange = " << geoRefineRange << std::endl;

        // The distance from the membrane to switch from linear to geometric refinement
        double geoRefineDistance = 10.0*dy_membrane + geoRefineRange;
        debug_jochen << "Switch distance from linear spacing to geometric refinement: "
            << geoRefineDistance << std::endl;

        double yrange_equidistant = params.yMemb() - geoRefineDistance;
        // number of equidistant element that fit into remaining space
        int nEquidistant = (int) std::floor(yrange_equidistant / max_dy_equidistant);
        //double dy_linear = yrange_linear / nLinear;

        // Handle case when dy was chosen too large to obtain a fine grid resolution near the membrane
        if(params.yMemb() - geoRefineDistance < params.yMin()
            || params.yMemb() + geoRefineDistance > params.yMax())
        {
          std::stringstream errorMsg;
          errorMsg << "Please check your grid parameters, cannot reach a sufficiently fine resolution of dy=";
          errorMsg << dy_membrane << " near the membrane with yMemb=" << params.yMemb();
          errorMsg << " when using a maximum dy=" << max_dy_equidistant << "!";
          DUNE_THROW(Dune::Exception, errorMsg.str());
        }
        debug_jochen << "Equidistant y-spacing till " << (params.yMemb() - geoRefineDistance)
                     << ", then geometric refinement" << std::endl;

        debug_jochen << "params.yMemb() = " << params.yMemb() << std::endl;
        debug_jochen << "params.dMemb() = " << params.dMemb() << std::endl;
        debug_jochen << "yMembEnd = " << yMembEnd << std::endl;

        // Add equidistant elements for electrolyte
        debug_jochen << nEquidistant << " elements of height " << max_dy_equidistant << std::endl;
        if(nEquidistant > 0)
        {
          y.equidistant_n_h(nEquidistant, max_dy_equidistant);
        }

        // Now refine geometrically towards membrane; enforce hend to by guaranteed
        debug_jochen << "geometric transition till " << (params.yMemb() - 10.0*dy_membrane) << std::endl;
        y.geometric_h0_hend_xend(max_dy_equidistant, dy_membrane, params.yMemb() - 10.0*dy_membrane);

        // Add equidistant elements close to membrane
        debug_jochen << "10 equidistant elements of height " << dy_membrane << " near membrane" << std::endl;
        y.equidistant_n_h(10, dy_membrane);

        // Add one membrane element layer
        debug_jochen << "membrane thickness " << params.dMemb() << std::endl;
        y.equidistant_n_h(1, params.dMemb());

        // In case of two membranes: Handle intracellular domain
        if(numMembranes == 2)
        {
          // Maximum grid size in y-direction within the cell
          max_dy_equidistant = params.dYCell();

          if(std::ceil(params.getCellWidth() / dy_membrane) < 20)
          {

            debug_jochen << params.getCellWidth() << std::endl;
            debug_jochen << std::ceil(params.getCellWidth() / dy_membrane) << std::endl;
            DUNE_THROW(Dune::Exception,
                "Minimum of 20 equidistant elements do not fit into intracellular domain!");
          }

          // Add equidistant elements close to membrane
          debug_jochen << "10 equidistant elements of height " << dy_membrane << " near membrane" << std::endl;
          y.equidistant_n_h(10, dy_membrane);

          double cellCenter = yMembEnd + 0.5*params.getCellWidth();
          double geoRefineRangeIntracellular = params.getCellWidth() - 20*dy_membrane;

          // Anders: Raum noch für mind. 2 max_dy Zellen => zwischendurch äquidistant
          // Sonst: Geometrisch bis zur Hälfte verfeinern und wieder spiegeln
          if(std::floor((geoRefineRangeIntracellular - 2.0 * geoRefineRange) / max_dy_equidistant) < 2.0)
          {
//            std::stringstream errorMsg;
//            errorMsg << "geoRefineRangeIntracellular (" << geoRefineRangeIntracellular
//                << ") < 2.0 * geoRefineRange (" << (2.0 * geoRefineRange) << ")";
//
//            DUNE_THROW(Dune::NotImplemented, errorMsg.str());

            y.geometric_n_h0_xend(n_transition, dy_membrane, cellCenter);
            y.geometric_n_hend_xend(n_transition, dy_membrane, yMembEnd + params.getCellWidth() - 10*dy_membrane);
          } else {
            double equidistantRangeIntracellular = params.getCellWidth() - 20*dy_membrane - 2.0*geoRefineRange;
            int nEquidistantIntracellular = (int) std::floor(equidistantRangeIntracellular / max_dy_equidistant);
            equidistantRangeIntracellular = nEquidistantIntracellular * max_dy_equidistant;

            double geoRefineEndIntracellular = yMembEnd + 0.5*(cellWidth-20.0*dy_membrane-equidistantRangeIntracellular);

            // Geometrically coarsen
            //y.geometric_n_h0_q(n_transition, dy_membrane, 1./q);
            y.geometric_h0_hend_xend(dy_membrane, max_dy_equidistant, geoRefineEndIntracellular, 1);

            // Equidistant elements in the middle if the cell
            y.equidistant_n_h(nEquidistantIntracellular, max_dy_equidistant);

            // Geometrically refine
            //y.geometric_n_h0_q(n_transition, max_dy_equidistant, q);
            y.geometric_h0_hend_xend(max_dy_equidistant, dy_membrane, yMembEnd+cellWidth-10*dy_membrane);
          }


          // Add equidistant elements close to membrane
          debug_jochen << "10 equidistant elements of height " << dy_membrane << " near membrane" << std::endl;
          y.equidistant_n_h(10, dy_membrane);

          // Add one membrane element layer
          debug_jochen << "membrane thickness " << params.dMemb() << std::endl;
          y.equidistant_n_h(1, params.dMemb());

          yMembEnd += params.getCellWidth() + params.dMemb();
        }

        // Store the volume of the first extracellular cell
        double min_extracellular_volume = (2*con_pi*params.dX()*y[y.size()-1]
                                           + con_pi*params.dX()*y[y.size()-1])*dy_membrane;

        // Now mirror the whole procedure on the other side of the membrane!
        // Add equidistant elements close to membrane
        debug_jochen << "10 equidistant elements of height " << dy_membrane << " near membrane" << std::endl;
        y.equidistant_n_h(10, dy_membrane);

        // Maximum grid size in y-direction for the extracellular space
        max_dy_equidistant = params.dY();

        // ============= NEW => test! ========================================================================
        // Recalculate geoRefineRange, max_dy_equidistent might have changed for this side of the membrane!
        n_transition = std::ceil(std::abs(std::log(dy_membrane / max_dy_equidistant) / std::log(q)));

        // range of the geometric series
        geoRefineRange = max_dy_equidistant * (std::pow(q, n_transition+1) - 1) / (q - 1);
        debug_jochen << "geoRefineRange = " << geoRefineRange << std::endl;
        geoRefineDistance = 10.0*dy_membrane + geoRefineRange;
        // ============= NEW => test! ========================================================================

        // Update ranges for the part above the membrane
        yrange_equidistant = params.yMax() - yMembEnd - geoRefineDistance;
        // number of equidistant element that fit into remaining space
        nEquidistant = (int) std::floor(yrange_equidistant / max_dy_equidistant);
        //dy_linear = yrange_linear / nLinear;

        bool useRadialCylinderRefinement = params.general.get("useRadialCylinderRefinement", false);
        if(not useRadialCylinderRefinement)
        {
          // Geometric coarsening
          debug_jochen << "geometric transition till " << (params.yMax() - (nEquidistant*max_dy_equidistant)) << std::endl;
          y.geometric_h0_hend_xend(dy_membrane, max_dy_equidistant, params.yMax() - (nEquidistant*max_dy_equidistant), 1);

          // Add linear elements until upper boundary is reached
          debug_jochen << nEquidistant << " elements of height " << max_dy_equidistant << std::endl;
          if(nEquidistant > 0)
          {
            y.equidistant_n_xend(nEquidistant, params.yMax());
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
          CylinderGridVector<BaseGrid> y_temp(params.dX());
          y_temp.push_back(y_volumeEquilibrationStart);

          // Calculate necessary parameters for radial refinement beforehand
          double remaining_volume = y_temp.remaining_volume(params.yMax());
          int nelements_y = std::ceil(remaining_volume / (min_extracellular_volume * max_volume_ratio));
          debug_jochen << "Remaining volume: " << remaining_volume << std::endl;
          debug_jochen << "Would need " << nelements_y << " to reach a max volume ratio of "
              << max_volume_ratio << "." << std::endl;
          nelements_y = std::min(nelements_y, max_nelements_y);
          debug_jochen << "Using " << nelements_y << " => volume ratio "
              << (remaining_volume / (nelements_y * min_extracellular_volume)) << std::endl;
          y_temp.equivolume_n_rend(nelements_y, params.yMax(), false, 1, false);
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

          remaining_volume = y.remaining_volume(params.yMax());
          nelements_y = std::ceil(remaining_volume / (min_extracellular_volume * max_volume_ratio));
          debug_jochen << "Remaining volume: " << remaining_volume << std::endl;
          debug_jochen << "Would need " << nelements_y << " to reach a max volume ratio of "
              << max_volume_ratio << "." << std::endl;
          nelements_y = std::min(nelements_y, max_nelements_y);
          debug_jochen << "Using " << nelements_y << " => volume ratio "
              << (remaining_volume / (nelements_y * min_extracellular_volume)) << std::endl;

          y.equivolume_n_rend(nelements_y, params.yMax(), false, 1);
        }
      } else {
        int n_y = std::pow(2,level);
        // Equidistant elements, number specified by refinement level
        if(params.useMembrane())
        {
          y.equidistant_n_xend(n_y,params.yMemb());
          // One layer of membrane elements of thickness dMemb
          y.equidistant_n_h(1, params.dMemb());
          // Equidistant elements, number specified by refinement level
          y.equidistant_n_xend(n_y,params.yMax());
        } else {
          // When not using membrane: when 'refineY...' flag is set, use fixed size dy from config file
          if(params.doRefineYDirectionGeometrically())
          {
            n_y = std::ceil((params.yMax() - params.yMin()) / params.dY());
            y.equidistant_n_xend(n_y,params.yMax());
          } else {
            y.equidistant_n_xend(n_y,params.yMax());
          }
        }
      }

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
      int n_x = 1;
      if(params.doRefineXDirection())
      {
        // TODO Specify different levels for x and y refinement?
        // TODO Change this at least according to the xMax/yMax ratio
        n_x = std::pow(2,level);

        double xrange = params.xMax() - params.xMin();
        double d_x = xrange / n_x;

        if(params.useMembrane())
        {
          double yrange_half = params.yMemb() - params.yMin();
          double rough_d_y = yrange_half / ((y.size() / 2.0) - 1);

          debug_jochen << "original d_x = " << d_x << std::endl;
          d_x = std::min(d_x, rough_d_y);
        }

        debug_jochen << "new d_x = " << d_x << std::endl;
        n_x = std::ceil(xrange / d_x);

      } else {
        double dx = params.dX();
        n_x = (int) std::ceil(params.xMax() / dx);
      }

      // Fill x range with n_x equidistant elements
      debug_jochen << "=> n_x = " << n_x << std::endl;
      x.equidistant_n_xend(n_x,params.xMax());


      // HACK: Insert some extra nodes
      if(params.general.get("doLocalRefineHack", false))
      {
        double stim_x = params.stimulation.get("position_x", 0.0);
        double stim_y = params.stimulation.get("position_y", 0.0);
        double d_stim = params.general.get("d_memb", 1.0);
        double dx_fine = params.dYMin(), dy_fine = params.dYMin();
        for(int i=1; i<(int) std::ceil(2 * d_stim/dx_fine); i++)
        {
          x.push_back(stim_x - d_stim + i*dx_fine);
        }
        for(int i=1; i<(int) std::ceil(2 * d_stim/dy_fine); i++)
        {
          y.push_back(stim_y - d_stim + i*dy_fine);
        }
        std::sort(x.begin(), x.end());
        std::sort(y.begin(), y.end());
      }

      // ========================= Grid generation ==============================
      // Now copy x,y gridvectors in to a std::vector which can be passed to the parameter class
      std::vector<double> x_std(x);
      std::vector<double> y_std(y);
      params.setX(x_std);
      params.setY(y_std);

      debug_jochen << "===== x: ===================" << std::endl;
      Output::printVector(x_std);
      debug_jochen << "============================" << std::endl;
      debug_jochen << "===== y: ===================" << std::endl;
      Output::printVector(y_std);
      debug_jochen << "============================" << std::endl;

      debug_jochen << "x.size() = " << x.size() << std::endl;
      debug_jochen << "y.size() = " << y.size() << std::endl;
      debug_jochen << "=> Will give grid of " << ((x.size()-1) * (y.size()-1)) << " cells!" << std::endl;

      if(params.doEquilibration() && x.size() > 2 && !params.doForceEquilibration())
        DUNE_THROW(Dune::Exception, std::string("The new data transfer only works with grids having exactly 1 element ")
           + std::string("in x-direction! Please consider adjusting your grid parameters."));

      //get geometry type
      const Dune::GeometryType gt(Dune::GeometryType::cube,2);

      Ax1GridTools::prepareGrid(factory, x, y, gt);

      //make grid
      Dune::shared_ptr<BaseGrid> pBaseGrid(factory.createGrid());
      BaseGrid& baseGrid = *pBaseGrid;

#else
#warning "Selected grid could not be found!"
#endif
      debug_info << "Grid creation finished, time elapsed: " << timer.elapsed() << "s" << std::endl;
      // ############## END grid generation ###############
      
      
      // ============== Multidomaingrid stuff ================================================================
      typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,2> > Grid;
      Grid grid(baseGrid,false);

      typedef Grid::LeafGridView GV;
      GV gv = grid.leafView();

      typedef GV::Codim<0>::Iterator ElementIterator;
      typedef GV::Codim<0>::EntityPointer ElementPointer;
      typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;

      // Custom mapper to map element indices to group indices
      ElementSubdomainMapper elemSubdomainMapper;
      ElementMapper elementMapper(gv);

      // Tag elements by type
      for(ElementIterator ep = gv.begin<0>(); ep != gv.end<0>(); ++ep)
      {
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
          // We are on the membrane{
          if(center[1] > yMemb and center[1] < yMemb+dMemb)
          {
            subdomainIndex = MEMBRANE;
          }

          if(numMembranes == 1)
          {
            // There is only one membrane => the rest is outside of the cell
            if(center[1] > yMemb+dMemb)
            {
              subdomainIndex = ES;
            }
          }
          // In case of a second membrane layer
          if(numMembranes == 2)
          {
            // Tag second membrane
            if(center[1] > yMemb+dMemb+cellWidth and center[1] < yMemb+dMemb+cellWidth+dMemb)
            {
              subdomainIndex = MEMBRANE;
            }
            // Tag this cell as extracellular if it is outside the two membrane layers
            if(center[1] < yMemb or center[1] > yMemb+dMemb+cellWidth+dMemb)
            {
              subdomainIndex = ES;
            }
          }
        } else {
          // If no membrane is present, everything is extracellular space
          subdomainIndex = ES;
        }
        //debug_jochen << "Tagging element @ " << center[1] << " as " << subdomainIndex << std::endl;

        // Store group of this element [M->G]
        elemSubdomainMapper.setGroup(elemIndex, subdomainIndex);
      }

      typedef Grid::SubDomainGrid SubDomainGrid;
      SubDomainGrid& elecGrid = grid.subDomain(0);
      SubDomainGrid& membGrid = grid.subDomain(1);
      typedef Grid::ctype ctype;

      typedef SubDomainGrid::LeafGridView SDGV;
      SDGV elecGV = elecGrid.leafView();
      SDGV membGV = membGrid.leafView();

      grid.startSubDomainMarking();
      for (GV::Codim<0>::Iterator eit = gv.begin<0>(); eit != gv.end<0>(); ++eit)
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
              DUNE_THROW(Dune::Exception, "Element could not be assigned to a pre-defined group!");
          }

        }
      grid.preUpdateSubDomains();
      grid.updateSubDomains();
      grid.postUpdateSubDomains();
      // ===============================================================================================

      // BEGIN TESTS ===================================================================================
      const int N = 50;

      // (0) "Naked" grid (UG without multidomaingrid)
      typedef typename BaseGrid::LeafGridView BaseGV;
      BaseGV basegv = baseGrid.leafView();
      double volume0 = 0.0;
      Dune::Timer timer_ug;
      for(int i=0; i<N; i++)
      {
        volume0 = 0.0;
        for (BaseGV::Codim<0>::Iterator eit = basegv.begin<0>(); eit != basegv.end<0>(); ++eit)
        {
          typedef typename BaseGV::Codim<0>::Iterator::Entity::Geometry GeoUG;
          GeoUG geo = eit->geometry();

          volume0 += geo.volume();

          //geo.integrationElement(geo.local(geo.center()));
          //geo.corner(0);
          //geo.corner(geo.corners()-1);

        }
      }
      timer_ug.stop();
      debug_info << "UG volume: " << volume0 << ", calculation time: " << timer_ug.elapsed() << "s." << std::endl;

      // (1) "Normal" 2D geometry
      double volume1 = 0.0;
      Dune::Timer timer_2d;
      for(int i=0; i<N; i++)
      {
        volume1 = 0.0;
        for (GV::Codim<0>::Iterator eit = gv.begin<0>(); eit != gv.end<0>(); ++eit)
        {
          typedef typename GV::Codim<0>::Iterator::Entity::Geometry Geo2D;
          Geo2D geo = eit->geometry();

          volume1 += geo.volume();

          //geo.integrationElement(geo.local(geo.center()));
          //geo.corner(0);
          //geo.corner(geo.corners()-1);

        }
      }
      timer_2d.stop();
      debug_info << "multidomaingrid volume: " << volume1 << ", calculation time: " << timer_2d.elapsed() << "s." << std::endl;

      // (2) Wrapped cylinder geometry
      double volume2 = 0.0;
      Dune::Timer timer_cyl;
      for(int i=0; i<N; i++)
      {
        volume2 = 0.0;
        for (GV::Codim<0>::Iterator eit = gv.begin<0>(); eit != gv.end<0>(); ++eit)
        {
          typedef Acme2CylGeometry<typename GV::Codim<0>::Iterator::Entity::Geometry> GeoCyl;
          //typedef Acme2CylGeometryWrapper<typename GV::Codim<0>::Iterator::Entity::Geometry> GeoCyl;
          GeoCyl geo(eit->geometry()); // wrap it


          volume2 += geo.volume();
          //geo.integrationElement(geo.local(geo.center()));
          //geo.corner(0);
          //geo.corner(geo.corners()-1);

        }
      }
      timer_cyl.stop();
      debug_info << "Cyl volume: " << volume2 << ", calculation time: " << timer_cyl.elapsed() << "s." << std::endl;

      double perc1 = 100 * (timer_2d.elapsed() - timer_ug.elapsed()) / timer_ug.elapsed();
      double perc2 = 100 * (timer_cyl.elapsed() - timer_ug.elapsed()) / timer_ug.elapsed();

      debug_info << "Diff multidomaingrid/UG time: " << perc1 << " %" << std::endl;
      debug_info << "Diff Cyl/UG time            : " << perc2 << " %" << std::endl;

      // END TESTS =====================================================================================

      double elapsed = timer.stop();
      debug_info << "Time elapsed: " << elapsed << " s" << std::endl;
    }
  } catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 0;
  }

}
