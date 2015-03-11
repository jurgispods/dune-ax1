/*
 * acme2_cyl_simulationloader.hh
 *
 *  Created on: Jul 10, 2013
 *      Author: jpods
 */

#ifndef DUNE_ACME2_CYL_SIMULATIONLOADER_HH
#define DUNE_ACME2_CYL_SIMULATIONLOADER_HH

#include <dune/ax1/acme2_cyl/common/acme2_cyl_solutionvectors.hh>

// Helper class for Acme2CylOutput managing loading a saved simulation state
template<typename Acme2CylTraits, typename Acme2CylOutputTraits>
class SimulationLoader
{
  public:
    typedef Acme2CylTraits Traits;
    typedef Acme2CylOutputTraits OutputTraits;

    typedef typename Traits::Real Real;
    typedef typename Traits::Physics PHYSICS;
    typedef typename Traits::GridView GV;
    typedef typename Traits::SubGridView SubGV;
    typedef typename Traits::MultiGFS MultiGFS;
    typedef typename Traits::GFS_CON_SUB GFS_CON; // SubGFS!
    typedef typename Traits::GFS_POT_SUB GFS_POT; // SubGFS!

    SimulationLoader(const PHYSICS& physics_,
        const GV& gv_,
        const SubGV& elecGV_,
        const GFS_POT& gfsPot_,
        const GFS_CON& gfsCon_,
        Acme2CylSolutionVectors<Traits>& solutionVectors_)
    : physics(physics_),
      gv(gv_),
      elecGV(elecGV_),
      gfsPot(gfsPot_),
      gfsCon(gfsCon_),
      solutionVectors(solutionVectors_)
    {}

    std::map<std::string,std::string> getLoadFileNames() const
    {
      std::map<std::string,std::string> loadFileNames;

      const typename Acme2CylParameters::MembraneGroups& membGroups = physics.getParams().getMembraneGroups();
      for(int i=0; i<membGroups.size(); i++)
      {
        std::string groupName = std::get<0>(membGroups[i]);

        std::string loadFilename = "dummy.dat";

        bool extrapolateXData = physics.getParams().doExtrapolateXData();

        // Only activate loading multiple files when extrapolating data; otherwise we always load from one single state
        bool loadMultipleStates = extrapolateXData && physics.getParams().membrane.get("loadMultipleStates", false);

        // When not loading multiple states: Use the loadFilename from 'general' section in config file!
        if(! loadMultipleStates && physics.getParams().general.hasKey("loadFilename"))
        {
          loadFilename = physics.getParams().general.get("loadFilename","dummy.dat");
        // When loading from multiple states, each membrane group must define its own loadFilename!
        } else if (loadMultipleStates && physics.getParams().membrane.sub(groupName).hasKey("loadFilename"))
        {
          loadFilename = physics.getParams().membrane.sub(groupName).get("loadFilename", "dummy.dat");
        } else
        {
          DUNE_THROW(Dune::Exception, "Did not find a 'loadFilename' entry for membrane group '" << groupName << "'!");
        }
        // Only add the filename if it does not yet exist in vector
		if(loadFileNames.count(groupName) == 0)
		{
		  loadFileNames[groupName] = loadFilename;
		}
      }
      return loadFileNames;
    }

    void getParallelFilename(std::string& filename) const
    {
      // Check if a state file for the current rank exists; if not, load from rank-0 file
      std::size_t index = filename.find(".dat");
      assert(index != std::string::npos);
      std::stringstream rank_filename;
      rank_filename << filename.substr(0, index)
          << "_p" << physics.gridView().comm().rank()
          << ".dat";

      std::ifstream in(rank_filename.str());
      if(! in.good())
      {
        rank_filename.str("");
        rank_filename << filename.substr(0, index)
          << "_p0.dat";
      }
      filename = rank_filename.str();
    }

    bool checkCoordinates(Ax1SimulationData<Real>& simulationData)
    {
      debug_jochen << "[SimLoader::checkCoordinates()]" << std::endl;
      bool extrapolateXData = physics.getParams().doExtrapolateXData();

      if(gv.size(0) != simulationData.getNElements()
          || solutionVectors.uold.flatsize() != simulationData.getVector("uold").data.size())
      {
        // Throw error when number of nodes does not match;
        // When extrapolateXData is desired, disable this check
        if(not extrapolateXData /*||
            (extrapolateXData &&  (gv.size(0) % simulationData.getNElements() != 0))*/)
        {
          DUNE_THROW(Dune::Exception,
            "Error [# nodes] loading saved simulation state from file '" << simulationData.getFileName()
            << "', was the saved state generated with the same grid and finite element order?");
        } else {
          debug_info << "Number of nodes for the data to be loaded does not match number of nodes in "
              << "the current grid. Will try to interpolate data now!" << std::endl;
        }
      }

      // Use a more accurate check and compare node positions!
      // In the case extrapolateXData is set, only check for y coordinate consistency!
      bool isGridCompatible = true;
      bool xCoordinatesMatch = true;
      bool yCoordinatesMatch = true;
      std::vector<Real> x, y;

      // NEW method (load separate gridvectors for x, y coordinates)
      if(simulationData.hasVector("x") && simulationData.hasVector("y"))
      {
        x = simulationData.getVector("x").data;
        y = simulationData.getVector("y").data;

        debug_jochen << "Loaded separate x, y coordinate vectors!" << std::endl;
        debug_jochen << "x:" << std::endl;
        Output::printVector(x);
        debug_jochen << "y:" << std::endl;
        Output::printVector(y);

        if(y != physics.getParams().Y())
        {
          yCoordinatesMatch = false;
        }
        // Grid must at least have the same extents
        if(y[0] != physics.getParams().Y()[0] || y.back() != physics.getParams().Y().back())
        {
          isGridCompatible = false;
        }

        if(x != physics.getParams().X())
        {
          xCoordinatesMatch = false;
          if(not extrapolateXData)
          {
            // Grid must at least have the same extents
            if(x[0] != physics.getParams().X()[0] || x.back() != physics.getParams().X().back())
            {
              isGridCompatible = false;
            }
          }
        }

        if(! isGridCompatible)
        {
          DUNE_THROW(Dune::Exception, "Error loading saved simulation state from file '" << simulationData.getFileName()
              << "', was the saved state generated with a compatible grid and same finite element order?"
              << " [Do coordinates match? x : " << xCoordinatesMatch << ", y: " << yCoordinatesMatch
              << "]");
        } else {
          debug_info << "Grids are compatible! "
            << " [Do coordinates match? x : " << xCoordinatesMatch << ", y: " << yCoordinatesMatch
            << "]" << std::endl;
        }
      } else {
        DUNE_THROW(Dune::Exception,
            "No vector containing x, y coordinate vectors could be found in loaded simulation data!");
      }

      return (xCoordinatesMatch && yCoordinatesMatch);
    }

    void transferData(std::vector<Ax1SimulationState<Real> >& states, const MultiGFS& multigfs_Equi,
        bool isRefinement)
    {
      // Temporary variable; all solution vector will be loaded into this single vector one after another
      std::vector<Real> stdVec;

      std::vector<Real>& x = states[0].getData().getVector("x").data;
      std::vector<Real>& y = states[0].getData().getVector("y").data;

      bool doYCoordinatesMatch = (y == physics.getParams().Y());

      // Extract sub GFS's
      typename Acme2CylTraits::GFS_CON_SUB gfsConSub_Equi(multigfs_Equi);
      typename Acme2CylTraits::GFS_POT_SUB gfsPotSub_Equi(multigfs_Equi);

      // Load solution vectors (one per loaded simulation state)
      std::vector<typename Acme2CylTraits::U> uold_Equi(states.size(), typename Acme2CylTraits::U(multigfs_Equi,0.0));
      std::vector<typename Acme2CylTraits::U> unew_Equi(states.size(), typename Acme2CylTraits::U(multigfs_Equi,0.0));

      // One potential/concentration gridfunction for each loaded solution vector
      std::vector<Dune::shared_ptr<typename OutputTraits::DGF_POT> > dgfOldPot_Equi;
      std::vector<Dune::shared_ptr<typename OutputTraits::DGF_POT> > dgfNewPot_Equi;

      std::vector<Dune::shared_ptr<typename OutputTraits::DGF_CON> > dgfOldCon_Equi;
      std::vector<Dune::shared_ptr<typename OutputTraits::DGF_CON> > dgfNewCon_Equi;

      // This vector contains the membrane groups associated to the loaded state files, i.e. 'equiMembraneGroups[j]' is
      // the membrane group belonging to loaded state 'states[j]'
      std::vector<int> equiMembraneGroups;

      for(int j=0; j<states.size(); j++)
      {
        stdVec = states[j].getData().getVector("uold").data;
        //uold_Equi.std_copy_from(stdVec);
        std::size_t i = 0;
        for (typename Acme2CylTraits::U::iterator it = uold_Equi[j].begin(); it != uold_Equi[j].end(); ++it)
        {
           *it = stdVec[i];
           i++;
        }
        stdVec = states[j].getData().getVector("unew").data;
        //unew_Equi.std_copy_from(stdVec);
        i = 0;
        for (typename Acme2CylTraits::U::iterator it = unew_Equi[j].begin(); it != unew_Equi[j].end(); ++it)
        {
           *it = stdVec[i];
           i++;
           //debug_jochen << "unew_Equi[" << i << "] = " << *it << std::endl;
        }

        // Create one set of gridfunctions for each membrane group here and let Ax1CoarseToFineGridTransferGridFunction
        // interpolate the values from two gridfunctions when the x-coordinate is 'between' two membrane groups!
        Dune::shared_ptr<typename OutputTraits::DGF_POT> dgfOldPot_Equi_ptr(
            new typename OutputTraits::DGF_POT(gfsPotSub_Equi, uold_Equi[j]));
        dgfOldPot_Equi.push_back(dgfOldPot_Equi_ptr);

        Dune::shared_ptr<typename OutputTraits::DGF_POT> dgfNewPot_Equi_ptr(
            new typename OutputTraits::DGF_POT(gfsPotSub_Equi, unew_Equi[j]));
        dgfNewPot_Equi.push_back(dgfNewPot_Equi_ptr);

        Dune::shared_ptr<typename OutputTraits::DGF_CON> dgfOldCon_Equi_ptr(
            new typename OutputTraits::DGF_CON(gfsConSub_Equi, uold_Equi[j]));
        dgfOldCon_Equi.push_back(dgfOldCon_Equi_ptr);

        Dune::shared_ptr<typename OutputTraits::DGF_CON> dgfNewCon_Equi_ptr(
            new typename OutputTraits::DGF_CON(gfsConSub_Equi, unew_Equi[j]));
        dgfNewCon_Equi.push_back(dgfNewCon_Equi_ptr);

        equiMembraneGroups.push_back(physics.getGroupIndexByName(states[j].getData().getMetaInfo("groupName")));

        debug_jochen << "Associating loaded state file #" << j << " ('" << states[j].getData().getFileName() << "') with group name '"
            << states[j].getData().getMetaInfo("groupName") << "' with group index " << equiMembraneGroups[j] << std::endl;
      }

      // NEW shit
      if(x.size() == 2 && (doYCoordinatesMatch || isRefinement))
      {
        int nNodesXCoarse = x.size();
        int nNodesYCoarse = y.size();
        int nNodesCoarse = nNodesXCoarse * nNodesYCoarse;

        bool doInterpolate = physics.getParams().membrane.get("interpolateMultipleStates", false);

        // Define functions which interpolate the equilibration values onto the new (finer) elecGV
        Ax1CoarseToFineGridTransferGridFunction<GV,typename OutputTraits::DGF_POT,PHYSICS>
          dgfOldPot_FineGrid(gv,dgfOldPot_Equi,physics,nNodesXCoarse-1,nNodesYCoarse-1,
          equiMembraneGroups,doInterpolate);
        //dgfOldPot_FineGrid.setVerbose(true);
        Ax1CoarseToFineGridTransferGridFunction<GV,typename OutputTraits::DGF_POT,PHYSICS>
          dgfNewPot_FineGrid(gv,dgfNewPot_Equi,physics,nNodesXCoarse-1,nNodesYCoarse-1,
          equiMembraneGroups,doInterpolate);
        //dgfNewPot_FineGrid.setVerbose(true);
        Ax1CoarseToFineGridTransferGridFunction<SubGV,typename OutputTraits::DGF_CON,PHYSICS>
          dgfOldCon_FineGrid(elecGV,dgfOldCon_Equi,physics,nNodesXCoarse-1,nNodesYCoarse-1,
          equiMembraneGroups,doInterpolate);
        Ax1CoarseToFineGridTransferGridFunction<SubGV,typename OutputTraits::DGF_CON,PHYSICS>
          dgfNewCon_FineGrid(elecGV,dgfNewCon_Equi,physics,nNodesXCoarse-1,nNodesYCoarse-1,
          equiMembraneGroups,doInterpolate);

        debug_jochen << "Doing interpolation: " << doInterpolate << std::endl;

        // Use the GFS belonging to the new (fine) grid here for interpolation!
        // This seems to work even when writing into the 'large' solution vectors
        debug_jochen << "Interpolating pot into uold (new method)..." << std::endl;
        Dune::PDELab::interpolate(dgfOldPot_FineGrid, gfsPot, solutionVectors.uold);
        debug_jochen << "Interpolating pot into unew (new method)..." << std::endl;
        Dune::PDELab::interpolate(dgfNewPot_FineGrid, gfsPot, solutionVectors.unew);
        debug_jochen << "Interpolating con into uold (new method)..." << std::endl;
        Dune::PDELab::interpolate(dgfOldCon_FineGrid, gfsCon, solutionVectors.uold);
        debug_jochen << "Interpolating con into unew (new method)..." << std::endl;
        Dune::PDELab::interpolate(dgfNewCon_FineGrid, gfsCon, solutionVectors.unew);

      } else {      // OLD shit
        // Fall back to this method in case there are no x, y coordinate vectors
        // OR when the loaded grid does have more than 1 element in x-direction!

        if(dgfOldPot_Equi.size() > 1)
          DUNE_THROW(Dune::Exception,
              "Multiple gridfunction interpolation not implemented when using old extrapolation method!");

        // The following (commented-out) stuff does not work when using a multidomain grid, as the
        // subdomain grids (on which the elements of dgfOldCon_Equi and dgfNewCon_Equi are living)
        // do not support level index sets, yielding an error when calling Dune::PDELab::interpolate.
//        // Make a function from the gridfunction living on the old (coarse) grid
//        typedef Dune::PDELab::GridFunctionToFunctionAdapter<typename OutputTraits::DGF_POT> F_POT;
//        typedef Dune::PDELab::GridFunctionToFunctionAdapter<typename OutputTraits::DGF_CON> F_CON;
//        F_POT functionOldPot(*dgfOldPot_Equi[0]);
//        F_POT functionNewPot(*dgfNewPot_Equi[0]);
//        F_CON functionOldCon(*dgfOldCon_Equi[0]); // TODO Make this a GF living in the multidomain grid
//        F_CON functionNewCon(*dgfNewCon_Equi[0]); // TODO Make this a GF living in the multidomain grid
//
//        // Now make again a gridfunction from this global-valued function, this time living on the new (fine) grid
//        Dune::PDELab::FunctionToGridFunctionAdapter<GV,F_POT> dgfOldPot_FineGrid(gv, functionOldPot);
//        Dune::PDELab::FunctionToGridFunctionAdapter<GV,F_POT> dgfNewPot_FineGrid(gv, functionNewPot);
//        Dune::PDELab::FunctionToGridFunctionAdapter<SubGV,F_CON> dgfOldCon_FineGrid(elecGV, functionOldCon);
//        Dune::PDELab::FunctionToGridFunctionAdapter<SubGV,F_CON> dgfNewCon_FineGrid(elecGV, functionNewCon);
//
//        // Use the GFS belonging to the new (fine) grid here for interpolation!
//        // This seems to work even when writing into the 'large' solution vectors
//        Dune::PDELab::interpolate(dgfOldPot_FineGrid, gfsPot, solutionVectors.uold);
//        Dune::PDELab::interpolate(dgfNewPot_FineGrid, gfsPot, solutionVectors.unew);
//        Dune::PDELab::interpolate(dgfOldCon_FineGrid, gfsCon, solutionVectors.uold);
//        Dune::PDELab::interpolate(dgfNewCon_FineGrid, gfsCon, solutionVectors.unew);

        // Define functions which interpolate the equilibration values onto the new (finer) gridview
        Ax1CoarseGridInterpolationGridFunction<GV,typename OutputTraits::DGF_POT,PHYSICS>
          dgfOldPot_FineGrid(gv,*dgfOldPot_Equi[0],physics);
        Ax1CoarseGridInterpolationGridFunction<GV,typename OutputTraits::DGF_POT,PHYSICS>
          dgfNewPot_FineGrid(gv,*dgfNewPot_Equi[0],physics);
        Ax1CoarseGridInterpolationGridFunction<SubGV,typename OutputTraits::DGF_CON,PHYSICS>
          dgfOldCon_FineGrid(elecGV,*dgfOldCon_Equi[0],physics);
        Ax1CoarseGridInterpolationGridFunction<SubGV,typename OutputTraits::DGF_CON,PHYSICS>
          dgfNewCon_FineGrid(elecGV,*dgfNewCon_Equi[0],physics);

        // Use the GFS belonging to the new (fine) grid here for interpolation!
        // This seems to work even when writing into the 'large' solution vectors
        debug_jochen << "Interpolating pot into uold (old method)..." << std::endl;
        Dune::PDELab::interpolate(dgfOldPot_FineGrid, gfsPot, solutionVectors.uold);
        debug_jochen << "Interpolating pot into unew (old method)..." << std::endl;
        Dune::PDELab::interpolate(dgfNewPot_FineGrid, gfsPot, solutionVectors.unew);
        debug_jochen << "Interpolating con into uold (old method)..." << std::endl;
        Dune::PDELab::interpolate(dgfOldCon_FineGrid, gfsCon, solutionVectors.uold);
        debug_jochen << "Interpolating con into unew (old method)..." << std::endl;
        Dune::PDELab::interpolate(dgfNewCon_FineGrid, gfsCon, solutionVectors.unew);
      }
      debug_jochen << "Data transfer completed!" << std::endl;
    }

    // Following: Hardcoded extrapolation for uold/unew. Think this is ok for our needs.
    /**
     * This is the most general approach to restore data from the old grid: Restore the old grid!
     * Then wrap the loaded solution vectors into grid functions and use the GF adapters from
     * ax1_gridfunctionadapters.hh to interpolate the DOFs onto the new grid. Works nicely!
     *
     * Attention: The gridfunction adapters assume the same ordering of elements in the loaded
     * gridview as in the one onto which the interpolation/extrapolation shall be done! In
     * addition, they are only able to handle a grid with exactly one element in x-direction.
     * The benefit is a greatly increased extrapolation procedure in comparison to the old
     * method, where the gridfunction adapters had to loop over the whole grid each time to find
     * a matching entity.
     */
    void extrapolateData(std::vector<Ax1SimulationState<Real> >& states)
    {
      std::vector<Real>& x = states[0].getData().getVector("x").data;
      std::vector<Real>& y = states[0].getData().getVector("y").data;

      bool doYCoordinatesMatch = (y == physics.getParams().Y());

      // Check coordinate compatibility of simulation states
      for(int i=1; i<states.size(); i++)
      {
        std::vector<Real>& x_other = states[i].getData().getVector("x").data;

        if(x_other.size() != x.size())
          DUNE_THROW(Dune::Exception, "You are trying to load from two states which have different numbers of x elements, this will not work!");
        else if(x.size() > 2 && x_other != x)
          DUNE_THROW(Dune::Exception, "Two states have non-matching x-coordinates, this will not work!");
      }

      if(x.size() != 2 || !doYCoordinatesMatch)
        debug_warn << "Equilibration grid has either more than 1 element in x-direction or y coordinates to not match, this might result "
          << "in a slow extrapolation!" << std::endl;

      if(physics.getParams().doLoadChannelStates())
        DUNE_THROW(Dune::NotImplemented,
            "Cannot load channel states when interpolating data from a different grid!");

#if USE_GRID == 1
      // Create YaspGrid manually
      std::ifstream in_yasp(physics.getParams().getLoadGridFileName());
      if(! in_yasp.good())
        DUNE_THROW(Dune::Exception, "Yasp DGF dile '"
            << physics.getParams().getLoadGridFileName() << "' could not be read!");

      std::string dummy("");
      while(dummy.find("INTERVAL") == std::string::npos && !in_yasp.eof())
        std::getline(in_yasp, dummy);
      debug_jochen << "dummy: " << dummy << std::endl;
      assert(dummy.find("INTERVAL") != std::string::npos);

      const int dim = 2;
      Dune::FieldVector<double,dim> L;
      std::getline(in_yasp, dummy); // ignore the "0 0" line in DGF file for lower left corner
      in_yasp >> L[0] >> L[1];
      debug_jochen << L << std::endl;

      // This check is absolutely necessary, otherwise the geometry transformation used by GeometryGrid
      // will fail to map vetex coordinates and throw a 'std::out_of_range error' for map::at!
      if(L[0] != x.back() || L[1] != y.back())
        DUNE_THROW(Dune::Exception, "Loaded grid from DGF file '" << physics.getParams().getLoadGridFileName()
          << "' does not match the state file ('" << physics.getParams().getLoadFilename()
          << "') coordinates! Grid file: xmax = " << L[0] << ", ymax = " << L[1]
          << ", state file: xmax = " << x.back() << ", ymax = " << y.back());


      Dune::array<int,dim> N;
      in_yasp >> N[0] >> N[1];
      debug_jochen << N << std::endl;
      std::bitset<dim> periodic;

      while(dummy.find("GRIDPARAMETER") == std::string::npos)
        std::getline(in_yasp, dummy);
      debug_jochen << "dummy: " << dummy << std::endl;
      assert(dummy.find("GRIDPARAMETER") != std::string::npos);

      std::string yasp_name;
      int overlap;
      in_yasp >> dummy >> yasp_name;
      in_yasp >> dummy >> overlap;

      // Set up custom load balancing
      int px = gv.comm().size();
      int py = 1; // Do not partition y-direction!

      // If equilibration grid consists of only one element in x-direction (the usual case),
      // do not partition on processors, but keep full grid on each processor!
      if(N[0] < gv.comm().size())
        px = N[0];

      // We don't need overlap when the loaded grid was generated by a single processor
      if(N[0] == 1)
        overlap = 0;

      if (px*py > 1 && px*py < gv.comm().size())
      {
        DUNE_THROW(Dune::Exception, "Loaded equilibration grid can not be partitioned with np = "
            << gv.comm().size() << " processors, it only has " << N[0] << " stripes in y-direction!");
      }

      typedef Ax1YaspPartition<2,Dune::FieldVector<int,2> > YP;
      Dune::FieldVector<int,2> yasppartitions;
      yasppartitions[0] = px;
      yasppartitions[1] = py;
      YP* yp = new YP(yasppartitions);
      if(gv.comm().rank() == 0)
      {
        debug_info << "Partitioning of YASP: " << yasppartitions << std::endl;
      }

      Dune::GridPtr<Dune::YaspGrid<2> > yaspPtr;
      if(N[0] == 1)
      {
        // Create sequential grid!
        yaspPtr = new Dune::YaspGrid<2>(L,N,periodic,overlap,yp);
      } else {
        yaspPtr = new Dune::YaspGrid<2>(gv.comm(),L,N,periodic,overlap,yp);
      }

      // Wrap loaded YaspGrid by GeometryGrid
      typedef Ax1TensorGridTransformation<double> CoordFunction;
      CoordFunction coordFunction(x.back(),y.back(),x,y);
      typedef Dune::GeometryGrid<BaseGrid,CoordFunction> HostGrid;
      HostGrid equilibrationGrid(*yaspPtr,coordFunction);

#else
      // Load grid used for equilibration from DGF file
      Dune::GridPtr<BaseGrid> gridptr(physics.getParams().getLoadGridFileName());

      debug_jochen << "Loaded old base grid from dgf file '"
          << physics.getParams().getLoadGridFileName() << "'." << std::endl;


      typedef BaseGrid HostGrid;
      HostGrid& equilibrationGrid = *gridptr;
#endif
      // Wrap into multidomain grid in order to be able to use all types/classes from Acme2CylSetup
      typename Acme2CylTraits::GridType md_equilibrationGrid(equilibrationGrid,false);

      //typedef typename Acme2CylTraits::GridType::LeafGridView EquiGV;
      GV equilibrationGV = md_equilibrationGrid.leafGridView();

      debug_jochen << "Loaded equilibration gridview has " << equilibrationGV.size(0) << " elements" << std::endl;

      if(equilibrationGV.size(0) != (x.size()-1)*(y.size()-1) )
        DUNE_THROW(Dune::Exception, "Loaded grid from DGF file '" << physics.getParams().getLoadGridFileName()
          << "' does not match the state file ('" << physics.getParams().getLoadGridFileName()
          << "') number of elements! Grid file: " << equilibrationGV.size(0)
          << ", state file: " << ((x.size()-1)*(y.size()-1)));

//        for (typename PHYSICS::ElementIterator eit
//            = equilibrationGV.template begin<0,Dune::PartitionIteratorType::Interior_Partition>();
//            eit != equilibrationGV.template end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++eit)
//        {
//          debug_jochen << "Equi Element @ " << eit->geometry().center()
//              << ", partitionType: " << eit->partitionType()
//              << std::endl;
//        }

      // Define a grid function which allows tagging of elements of the just loaded equilibration gridview
      Ax1SubdomainIndexGridFunction<GV,Real,PHYSICS> subdomainIndexGF(gv, physics);
      Ax1FineGridRestrictionGridFunction<GV,Ax1SubdomainIndexGridFunction<GV,Real,PHYSICS>,
          PHYSICS> subdomainIndexGF_Equi(equilibrationGV,subdomainIndexGF,physics,2); // maxFineElements=2!

      // Special case (and dirty hack) when loading an equilibration grid with 1 element in x-direction:
      // Do not care about x positions, just interpolate from fine elements with matching
      // y positions!
      Ax1FineToCoarseGridTransferGridFunction<GV,Ax1SubdomainIndexGridFunction<GV,Real,PHYSICS>,
          PHYSICS> subdomainIndexGF_Equi_oneXElem(equilibrationGV,subdomainIndexGF,physics,
              1, y.size()-1); // nElementsXCoarse = 1, nElementsYCoarse = y.size()-1; is this general enough?

      bool isRefinement = subdomainIndexGF_Equi_oneXElem.isCompatible();
      debug_jochen << "Current grid is a refinement of the loaded grid. Using new extrapolation method!"
          << std::endl;

      debug_jochen << "Beginning marking equilibration md grid" << std::endl;
      // Now reconstruct the MultiDomainGrid subdoamin structure!
      typedef typename Acme2CylTraits::GridType::SubDomainGrid SubDomainGrid;
      SubDomainGrid& elecGrid_Equi = md_equilibrationGrid.subDomain(0);
      SubDomainGrid& membGrid_Equi = md_equilibrationGrid.subDomain(1);
      typedef typename SubDomainGrid::LeafGridView SDGV;

      SDGV elecGV_Equi = elecGrid_Equi.leafGridView();
      SDGV membGV_Equi = membGrid_Equi.leafGridView();

      md_equilibrationGrid.startSubDomainMarking();
      typename Ax1FineGridRestrictionGridFunction<GV,Ax1SubdomainIndexGridFunction<GV,Real,PHYSICS>,
          PHYSICS>::Traits::DomainType local_center(0.5);
      for (typename PHYSICS::ElementIterator eit
          = equilibrationGV.template begin<0,Dune::PartitionIteratorType::Interior_Partition>();
          eit != equilibrationGV.template end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++eit)
      {

        typename Ax1FineGridRestrictionGridFunction<GV,Ax1SubdomainIndexGridFunction<GV,Real,PHYSICS>,
                PHYSICS>::Traits::RangeType mySubdomainIndex(-5);


        // New approach: Equilibration grid has exactly 1 element in x-direction
        if(x.size() == 2 && (doYCoordinatesMatch || isRefinement))
        {
          // Evaluate restriction gridfunction at center to get subdomainIndex for the coarse entity
          subdomainIndexGF_Equi_oneXElem.evaluate(*eit,local_center,mySubdomainIndex);
        } else {
          // Fall back to old (slow) extrapolation strategy
          subdomainIndexGF_Equi.evaluate(*eit,local_center,mySubdomainIndex);
        }

        // Convert double to int
        int subdomainIndex = mySubdomainIndex[0];

        //debug_jochen << "Marking element @" << eit->geometry().center() << " as " << subdomainIndex << std::endl;

        // Init subdomains
        switch(subdomainIndex)
        {
          case CYTOSOL:
            md_equilibrationGrid.addToSubDomain(0,*eit);
            break;
          case ES:
            md_equilibrationGrid.addToSubDomain(0,*eit);
            break;
          case MEMBRANE:
            md_equilibrationGrid.addToSubDomain(1,*eit);
            break;
          default:
            DUNE_THROW(Dune::Exception, "Element could not be assigned to a pre-defined subdomain!");
        }
      }
      md_equilibrationGrid.preUpdateSubDomains();
      md_equilibrationGrid.updateSubDomains();
      md_equilibrationGrid.postUpdateSubDomains();
      debug_jochen << "Equilibration md grid restored. Transferring data to new grid..." << std::endl;

      // Now we have a fully reconstructed MultiDomainGrid; we can use the types from Acme2CylTraits now
      // to load the concentration DOFs and interpolate them onto the new grid!

      // Reconstruct pot GFS
      typename Acme2CylTraits::FEM_POT femPot(equilibrationGV);
      typename Acme2CylTraits::GFS_POT gfsPot_Equi_TEMP(equilibrationGV,femPot);

      // Reconstruct con GFS
      typename Acme2CylTraits::FEM_CON femCon(equilibrationGV);
      typename Acme2CylTraits::GFS_SINGLE_CON gfsCon_Equi_Single_TEMP(elecGV_Equi,femCon);
      typename Acme2CylTraits::GFS_CON gfsCon_Equi_TEMP(gfsCon_Equi_Single_TEMP);

      // Reconstruct MultiGFS
#ifdef MULTIPLE_MEMBRANE_ELEMENTS

      typename Acme2CylTraits::OrderingTag orderingTag;
      typename Acme2CylTraits::MultiGFS multigfs_Equi(md_equilibrationGrid, orderingTag,
          gfsCon_Equi_TEMP, gfsPot_Equi_TEMP);

      // This is necessary in order to have the same DOF ordering in the equilibration grid...
      // assuming of course the equilibration grid used the same permutation strategy!
      if(physics.getParams().doReorderMatrix())
      {
        multigfs_Equi.update();

        std::vector<size_t> permutation(gfsCon_Equi_TEMP.size() + gfsPot_Equi_TEMP.size());
        std::vector<size_t> inv_permutation(gfsCon_Equi_TEMP.size() + gfsPot_Equi_TEMP.size());

        physics.setupDOFPermutation(multigfs_Equi, permutation, inv_permutation);
        //multigfs_Equi.orderingTag().updatePermutation(inv_permutation);

        // Following: new method using Steffen's PermutedOrderingTag
        std::vector<std::size_t>& tagPerm = multigfs_Equi.orderingTag().permutation();
        if(tagPerm.size() != inv_permutation.size())
        {
          debug_jochen << "[MME permutation] Permuted ordering vector has to be resized from size " << tagPerm.size()
              << " to size " << inv_permutation.size() << std::endl;
          tagPerm.resize(inv_permutation.size());
        }
        for(int i=0; i<tagPerm.size(); ++i)
        {
          tagPerm[i] = inv_permutation[i];
        }
        debug_jochen << "[MME permutation] Updated ordering permutation!" << std::endl;
      }

      // Check calculated permutation
//        for(int i=0; i<permutation.size(); i++)
//        {
//          debug_jochen << "perm[" << i << "] = " << permutation[i] << " -- " << "inv_perm[" << i << "] = " << inv_permutation[i] << std::endl;
//          assert(inv_permutation[permutation[i]] == i);
//        }
#else
#ifdef USE_MORI_OPERATOR_SPLIT
      // Use LexicographicOrderingTag
      typename Acme2CylTraits::MultiGFS multigfs_Equi(md_equilibrationGrid,
          gfsCon_Equi_TEMP, gfsPot_Equi_TEMP);
#else
      // Use InterleavedOrdering
      typename Acme2CylTraits::MultiGFS multigfs_Equi(md_equilibrationGrid, {NUMBER_OF_SPECIES,1},
          gfsCon_Equi_TEMP, gfsPot_Equi_TEMP);
      // Fallback call when using LexicographicOrderingTag
      //typename Acme2CylTraits::MultiGFS multigfs_Equi(md_equilibrationGrid,
      //      gfsCon_Equi_TEMP, gfsPot_Equi_TEMP);
#endif
#endif

      // Do actual data transfer
      transferData(states, multigfs_Equi,isRefinement);

      //Output::printSingleCoefficientVector(uold_Equi, "COARSE");
      //debug_jochen << " ====> " << std::endl;
      //Output::printSingleCoefficientVector(solutionVectors.uold, "FINE");
    }

    void truncateAndDeleteExistingFiles(int loadTimeStep, double time, int nExistingOutputFiles,
        const std::vector<std::string>& gnuplotFiles,
        const std::string& vtkBasenameDomain, const std::string& vtkBasenameMembrane,
        const std::string& checkpointBasename, const std::string& hdf5Basename,
        const std::string& filePrefix)
    {
      debug_jochen << "# old files: " << nExistingOutputFiles << std::endl;

      if(gv.comm().rank() == 0)
      {
        std::stringstream suffix("");

        // Clean up output directories: Remove all files that have been written between the checkpoint
        // we are loading from was written and the end of the actual end of that previous simulation
        for(int i=loadTimeStep+1; i<=nExistingOutputFiles; i++)
        {
          // Do it the C++ way!
          suffix.str("");
          suffix << std::setw(5) << std::setfill('0') << i << "*";
          //sprintf(suffix,"%05d",i);

          // Remove old VTK files
          std::string command = "rm " + vtkBasenameDomain + "-" + suffix.str();
          debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
          int unusedVariable = std::system(command.c_str());

          command = "rm " + vtkBasenameMembrane + "-" + suffix.str();
          debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
          unusedVariable = std::system(command.c_str());
        }
        for(int i=loadTimeStep+1; i<=nExistingOutputFiles; i++)
        {
          suffix.str("");
          suffix << std::setw(5) << std::setfill('0') << i << "*";

          // Remove spurious checkpoint files
          std::string command = "rm " + checkpointBasename + "-" + suffix.str();
          debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
          int unusedVariable = std::system(command.c_str());
        }
        for(int i=loadTimeStep+1; i<=nExistingOutputFiles; i++)
        {
          suffix.str("");
          suffix << std::setw(5) << std::setfill('0') << i << "*";

          // Remove spurious HDF5 files
          std::string command = "rm " + hdf5Basename + "-" + suffix.str();
          debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
          int unusedVariable = std::system(command.c_str());
        }

         // Special case gnuplot output: Since output for all timesteps is written into one single file,
         // we have to parse through all gnuplot files and search for the line where the first timestep
         // after the last checkpoint output has been written; strip everything following from that line
         // from the file

        // Take this as the tolerance when checking for equality of timestep values
         Real tol = 0.1 * (physics.getParams().getMinTimeStep() / physics.getTimeScale());
         for(int i=0; i<gnuplotFiles.size(); i++)
         {
           debug_jochen << "Processing file '" << gnuplotFiles[i] << "'..." << std::endl;
           int truncateFrom = -1;

           std::ifstream file_in(filePrefix + gnuplotFiles[i]);
           std::string prefix("# time: ");
           // Search for line with the matching time
           double searchTime = time;
           int line = Tools::findLine(file_in,searchTime,prefix,true,0,tol);
           debug_jochen << "Found value " << searchTime << " in line " << line << std::endl;

           if(line > -1)
           {
             // Find the beginning of next time output in this file_in
             double foundValue = -1;
             line = Tools::findLine(file_in,foundValue,prefix,false,line);
             debug_jochen << "Found value " << foundValue << " in line " << line << std::endl;

             if(line > -1)
               assert(foundValue > time);

             // Truncate file from line 'line' on!
             truncateFrom = line;

           } else {
             debug_verb << "Could not find a matching time step line in gnuplot file_in '"
                << gnuplotFiles[i] << "', assuming it does not have comment lines and truncating!"
                << std::endl;

             double foundValue = -1;
             // Assume that this file is free of '# time: ...' comment lines!
             assert(Tools::findLine(file_in,foundValue,prefix,false) == -1);

             // Truncate file from line belonging to timestep #loadTimeStep+1' on!
             truncateFrom = loadTimeStep+2;
           }

           debug_jochen << "Truncating file '" << gnuplotFiles[i] << "' from line " << truncateFrom << std::endl;

           // Create indermediate file and copy line 1 to line-1 into it
           file_in.close();
           file_in.open((filePrefix + gnuplotFiles[i]).c_str());

           std::string strLine("");
           std::ofstream file_out(filePrefix + gnuplotFiles[i] + ".new");

           int j=0;
           while(file_in.good())
           {
             std::getline(file_in,strLine);
             j++;
             // Note: if truncatFrom==-1, the entire file is read
             if(j == truncateFrom)
               break;
             file_out << strLine << std::endl;
           }

           file_in.close();
           file_out.close();

           // Truncate gnuplot files
           std::string command = "mv " + filePrefix + gnuplotFiles[i]
               + " " + filePrefix + gnuplotFiles[i] + ".old; "
               + "mv " + filePrefix + gnuplotFiles[i] + ".new " + filePrefix + gnuplotFiles[i];
           debug_jochen << "Truncating old file with command '" << command << "'" << std::endl;
           int unusedVariable = std::system(command.c_str());
         }
      }

    }


  private:
    const PHYSICS& physics;
    const GV& gv;
    const SubGV& elecGV;
    const GFS_POT& gfsPot;
    const GFS_CON& gfsCon;
    Acme2CylSolutionVectors<Traits>& solutionVectors;
};
#endif /* DUNE_ACME2_CYL_SIMULATIONLOADER_HH */
