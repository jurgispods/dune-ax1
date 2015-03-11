/*
 * hdf5_tools.hh
 *
 *  Created on: Feb 7, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_HDF5_TOOLS_HH
#define DUNE_AX1_HDF5_TOOLS_HH

#include <hdf5.h>

#include <dune/ax1/common/ax1_output2d.hh>
#include <dune/ax1/common/ax1_outputfunctors.hh>

template<typename OutputTraits>
class HDF5Tools
{
public:
	typedef unsigned int UINT;
	typedef double REAL;

	template<typename PHYSICS>
  static void initialize(const std::string& filename, const double time,
      const PHYSICS& physics)
  {
	  std::vector<std::string> groups(3);
	  groups[0] = "DOMAIN";
	  groups[1] = "MEMBRANE";
	  groups[2] = "ELEC";
	  HDF5Tools<OutputTraits>::initialize(filename, time, physics, groups);
  }


	template<typename PHYSICS>
	static void initialize(const std::string& filename, const double time,
	    const PHYSICS& physics,
	    const std::vector<std::string>& groups)
	{
	  herr_t status;             /* Generic return value */
#if USE_PARALLEL == 1
	  MPI_Comm communicator = physics.gridView().comm();

	  //Info variable needed for the HDF5
	  MPI_Info mpiInfo = MPI_INFO_NULL;

    // Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate( H5P_FILE_ACCESS );

    // Set up file access property list with parallel I/O access
    status = H5Pset_fapl_mpio(plist_id, communicator, mpiInfo);  //collective MPI!!! needs to be called on each processor of the communicator
#else
    hid_t plist_id = H5P_DEFAULT;
#endif
    assert(plist_id>-1);


		/* Create a new file using default properties. */
		hid_t file_id = H5Fcreate(
		  filename.c_str(),    // name of the file to be created
		  H5F_ACC_TRUNC,          // if you are trying to create a file that exists already, the existing file will be truncated, i.e., all data stored on the original file will be erased
		  H5P_DEFAULT,            // default file creation property list
		  plist_id);            // default file access property list

		assert( file_id > -1 );
    H5Pclose(plist_id);

		// Create groups (folders)
		for(int i=0; i<groups.size(); i++)
		{
		  hid_t group_id = H5Gcreate(file_id, groups[i].c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		  assert(group_id > -1);
		  //debug_jochen << "Created group " << groups[i] << std::endl;

		  htri_t exists = H5Lexists(file_id, groups[i].c_str(), H5P_DEFAULT);
      assert(exists > -1);

      status = H5Gclose(group_id);
      assert( status > -1 );
		}

    // Close the file.
    status = H5Fclose( file_id );
    assert( status > -1 );

		// Add meta information
		std::vector<int> dims(1);
		dims[0] = 1;
		std::string dataset_name("/time");
		std::vector<double> time_vec(1);
		time_vec[0] = time;
#if USE_PARALLEL == 1
		std::vector<int> global_dims = dims;
		std::vector<int> offsets(1);
		offsets[0] = 0;
		HDF5Tools<OutputTraits>::writeVector_MPI(filename, dims, global_dims, offsets, time_vec,
		    dataset_name, physics.gridView().comm());
#else
		HDF5Tools<OutputTraits>::writeVector(filename, dims, time_vec, dataset_name);
#endif
	}

	template<typename PHYSICS, typename GF>
  static void writeGF(const PHYSICS& physics, const GF& gf,
      const int SUBSAMPLING_POINTS, const std::string filename, std::string dataset)
	{
	  Ax1OutputFunctors::IdentityFunctor f;
	  HDF5Tools<OutputTraits>::writeGF(physics,gf,SUBSAMPLING_POINTS,filename,dataset,f);
	}

	template<typename PHYSICS, typename GF, typename FUNC>
	static void writeGF(const PHYSICS& physics, const GF& gf,
	    const int SUBSAMPLING_POINTS, const std::string filename, std::string dataset, FUNC& f)
	{
	  Dune::Timer timer;
	  std::vector<typename PHYSICS::Element::Geometry::GlobalCoordinate> pos;
	  std::vector<typename GF::Traits::RangeType> sol;

	  debug_jochen << "  [HDF5Tools] Writing dataset " << dataset << std::endl;

	  typename Ax1Output2D<OutputTraits>::SolutionVectorInfo info;
	  Ax1Output2D<OutputTraits>::getSolutionVector(physics, gf, SUBSAMPLING_POINTS, pos, sol, info);

	  // Apply user-defined functor to solution vector
	  sol = f(sol);

    //debug_jochen << "sol.size() = " << sol.size() << std::endl;

    // Extract number of points in each coordinate direction
    std::vector<typename GF::Traits::DomainFieldType> xcoord, ycoord;
    Tools::extractComponentFromVector(pos,xcoord,0);
    Tools::extractComponentFromVector(pos,ycoord,1);

    // Global domain size (number of elements)
    const int dim = GF::Traits::GridViewType::dimension;
    std::vector<int> global_dims(dim);
    global_dims[0] = physics.getParams().X().size() - 1;
    global_dims[1] = physics.getParams().Y().size() - 1;

    // Offset of this processor partition with respect to full domain
    std::vector<int> offsets(dim);
    offsets[0] = physics.nOffset(0);
    offsets[1] = physics.nOffset(1);

    std::vector<int>& dims = info.dimensions;
    std::string path = "";
    // Domain data (values on every element)
    int gridDomain = physics.getDomain(gf);
    if(gridDomain == GridDomains::DOMAIN_ALL)
    {
      path = "/DOMAIN/";
    }
    // Membrane data (only one layer of elements in y-direction)
    if(gridDomain == GridDomains::DOMAIN_MEMB_INTERFACE)
    {
      path = "/MEMBRANE/";
      global_dims[1] = physics.getParams().nMembranes();
      offsets[1] = 0;
    }
    // 'Real' membrane element data (is there any?)
    if(gridDomain == GridDomains::DOMAIN_MEMB)
    {
      path = "/MEMBRANE/";
      global_dims[1] = physics.getParams().nMembraneElements() * physics.getParams().nMembranes();
      offsets[1] = 0;
    }
    // Electrolyte data (subtract one membrane element in y-direction)
    if(gridDomain == GridDomains::DOMAIN_ELEC)
    {
      path = "/ELEC/";
      if(physics.getParams().useMembrane())
      {
        global_dims[1]-= physics.getParams().nMembranes() * physics.getParams().nMembraneElements();
      }
    }

    //debug_jochen << "dims[0] = " << dims[0] << std::endl;
    //debug_jochen << "dims[1] = " << dims[1] << std::endl;
    assert(dims[0]*dims[1] == sol.size());

    // Multiply directional offsets by number of values per element
    int nValuesPerElement = info.nValuesPerElement[0] * info.nValuesPerElement[1];
    if(nValuesPerElement > 1)
    {
      //debug_jochen << "nValues X = " << info.nValuesPerElement[0]
      //  << ", nValuesY = " << info.nValuesPerElement[1] << std::endl;
      offsets[0] *= info.nValuesPerElement[0];
      offsets[1] *= info.nValuesPerElement[1];

      global_dims[0] *= info.nValuesPerElement[0];
      global_dims[1] *= info.nValuesPerElement[1];
    }

    // There seems to already be a path in the dataset name!
    if(dataset.find("/") != std::string::npos)
    {
      path = "";
    }

    // Add path to dataset
    dataset = path + dataset;

    double tElapsed1 = 0;
    double tElapsed2 = 0;
    double tHdf5 = 0;

    // Watch out! We got a badass vector over here
    if(sol[0].size() > 1)
    {
      for(int k=0; k<sol[0].size(); k++)
      {
        std::vector<typename GF::Traits::RangeFieldType> sol_comp;
        Tools::extractComponentFromVector(sol,sol_comp,k);

        std::stringstream dataset_comp;
        dataset_comp << dataset;
        dataset_comp << "_" << k;

        tElapsed1 = timer.elapsed();
#if USE_PARALLEL == 1
//        debug_jochen
//          << " local_dims[0] = " << dims[0]
//          << ", local_dims[1] = " << dims[1]
//          << ", offset[0] = " << offsets[0]
//          << ", offset[1] = " << offsets[1] << std::endl;
        HDF5Tools<OutputTraits>::writeVector_MPI(filename, dims, global_dims, offsets, sol_comp,
            dataset_comp.str(), physics.gridView().comm());

#else
        HDF5Tools<OutputTraits>::writeVector(filename, dims, sol_comp, dataset_comp.str());
#endif
        tElapsed2 = timer.elapsed();
        tHdf5 += (tElapsed2 - tElapsed1);
      }
    } else {
      tElapsed1 = timer.elapsed();
#if USE_PARALLEL == 1
//      debug_jochen
//        << " local_dims[0] = " << dims[0]
//        << ", local_dims[1] = " << dims[1]
//        << ", offset[0] = " << offsets[0]
//        << ", offset[1] = " << offsets[1] << std::endl;
      HDF5Tools<OutputTraits>::writeVector_MPI(filename, dims, global_dims, offsets, sol,
          dataset, physics.gridView().comm());
#else
      HDF5Tools<OutputTraits>::writeVector(filename, dims, sol, dataset);
#endif
      tElapsed2 = timer.elapsed();
      tHdf5 += (tElapsed2 - tElapsed1);
	  }

    tElapsed1 = timer.elapsed();
    //TODO Put coordinate vectors into subfolders according to their output strategy (conforming, nonconforming)
#if USE_PARALLEL == 1
    HDF5Tools<OutputTraits>::writeVector_MPI(filename, dims, global_dims, offsets, xcoord,
        path+"x", physics.gridView().comm());
    HDF5Tools<OutputTraits>::writeVector_MPI(filename, dims, global_dims, offsets, ycoord,
        path+"y", physics.gridView().comm());
#else
    HDF5Tools<OutputTraits>::writeVector(filename, dims, xcoord, path + "x");
    HDF5Tools<OutputTraits>::writeVector(filename, dims, ycoord, path + "y");
#endif
    tElapsed2 = timer.elapsed();
    tHdf5 += (tElapsed2 - tElapsed1);

    debug_jochen << "    - Output time: " << timer.elapsed() << "s"
        << " (HDF5 I/O part: " << tHdf5 << "s / getSolutionVector part: "
        << info.tElapsed << "s)" << std::endl;
	}

	template<typename V>
  static void writeVector(
    	const std::string& filename,
			const std::vector<int>& dimensions,
			const V& data,
			const std::string& data_name)
	{
	  // general HDF5 status return value
    herr_t status;

    // get the dimensionality of the data
    UINT dim=dimensions.size();

    //number of elements the data vector should have
    UINT n=1;
    for(UINT i=0;i<dim;i++)
      n*=dimensions[i];

    //the data make no sense!
    assert(n==data.size());

    // Open file with r/w access
    hid_t file_id = H5Fopen(
      filename.c_str(),     // name of the file to be opened
      H5F_ACC_RDWR,         // open with read/write access
      H5P_DEFAULT);         // default file access property list

    assert( file_id > -1 );

    /* Check if the dataset already exists */
    // Assume intermediate groups have been created before, check full path to the actual dataset!
    htri_t exists = H5Lexists(file_id, data_name.c_str(), H5P_DEFAULT);
    assert(exists > -1);

    if(exists > 0)
    {
      //debug_jochen << "[HDF5Tools::writeVector()] Dataset '" << data_name
      //    << "' already exists, aborting write process." << std::endl;

      // Don't forget to close the file!
      status = H5Fclose( file_id );
      assert( status > -1 );

      return;
    }

    // create the memspace
    hsize_t mdims[1];
    mdims[0] = n;

    /* Create the dataspace for the dataset.
     * The dataspace describes the dimensions of the dataset array.
     */
    hsize_t dims[ dim ];
    for(UINT i=0;i<dim;i++)
      dims[i]=dimensions[i];

    hid_t dataspace_id = H5Screate_simple(
                        dim     // number of dimensions = rank
                        , dims  // vector containing sizes per dimension
                        , NULL  // maxdims == dims
                         );
    assert( dataspace_id > -1 );


    hid_t plist_id = H5Pcreate(                    // The new property list is initialized with default values for the specified class.
                   H5P_DATASET_CREATE  // Properties for dataset creation
                             );
    assert( plist_id > -1 );

    /*
     * chunks
     * !!! not optimized. Might be better if larger! -> use the inputdata
     */
    hsize_t chunk_dims[ dim ];
    for(UINT i=0;i<dim;i++)
      chunk_dims[i]=dims[i];  //maybe to small for fast saving!

    // set the chunk size!
    status = H5Pset_chunk(
                plist_id
                , dim             // must be == rank of the dataset
                , chunk_dims      // The values of the check_dims array define the size of the chunks to store the dataset's raw data. The unit of measure for check_dims values is dataset elements.
                 );
    assert( status > -1 );

    status = H5Pset_shuffle( plist_id ); // Sets the shuffle filter, H5Z_FILTER_SHUFFLE, in the dataset creation property list. This re-orders data to simplify compression.
    assert( status > -1 );

    status = H5Pset_deflate( plist_id, 1 ); // Sets deflate (GNU gzip) compression method and compression level. ( 0 < level < 9, lower = faster, but less compression )
    assert( status > -1 );

    /* Create the dataset. */
    hid_t dataset_id = H5Dcreate(
                   file_id             // Location identifier: id of the file or the group within which to create the dataset
                   , data_name.c_str()  // Dataset name: may be either an absolute path in the file or a relative path from file_id naming the dataset.
                   , H5T_NATIVE_DOUBLE   // Datatype identifier, here: IEEE floating point 32-bit little-endian
                   , dataspace_id      // Dataspace identifier
                   , H5P_DEFAULT
                   , plist_id          // Dataset creation property list identifier
                   , H5P_DEFAULT
                    );
    assert( dataset_id > -1 );

    /* Write the dataset. */
    status = H5Dwrite(                         // Writes raw data from a buffer to a dataset.
              dataset_id               // dataset identifier
              , H5T_NATIVE_DOUBLE      // memory datatype id
              , dataspace_id            // specifies the memory dataspace and the selection within it
              , H5S_ALL                // specifies the selection within the file dataset's dataspace. H5S_ALL indicates that the entire file dataspace, as defined by the current dimensions of the dataset, is to be selected
              , H5P_DEFAULT            // Identifier of a transfer property list for this I/O operation. H5P_DEFAULT: The default data transfer properties are used.
              , &(data[0])    // application memory buffer
                           );
    assert( status > -1 );

    /*
     * close everything
     */
    status = H5Sclose( dataspace_id );
    assert( status > -1 );

    status = H5Dclose( dataset_id );
    assert( status > -1 );

    status = H5Pclose( plist_id );
    assert( status > -1 );

    status = H5Fclose( file_id );
    assert( status > -1 );
	}


  /** function to write a vector (parallel) to a hdf5 file
   *
   * \tparam gobal_dim the global dimension of the stored data (total size)
   * \tparam data data which will be written to the file
   * \tparam local_count give the size of the local data
   * \tparam local_offset the offset of the data (in each direction)
   * \tparam helper the DUNE MPIHelper for the MPI communication
   * \tparam data_name is the name/path where the data are stored in the HDF5 file
   * \tparam data_filename is the filename of the data file
   *
   *
   See
   http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
   for a documentation of the HDF5 API.
   *
   */
	template<typename V>
  static void writeVector_MPI(const std::string& filename,
                              const std::vector<int>& local_dims,
                              const std::vector<int>& global_dims,
                              const std::vector<int>& offsets,
                              const V& data,
                              const std::string& data_name,
                              MPI_Comm communicator)
  {
	  //debug_jochen << "[writeVector_MPI] data.size() = " << data.size() << std::endl;
	  //debug_jochen << "[writeVector_MPI] data_name = " << data_name << std::endl;

    //Info variable needed for the HDF5
    MPI_Info mpiInfo = MPI_INFO_NULL;
    herr_t status = -1;             /* Generic return value */

    //get the dimension of the problem -> no checking if the given data makes any sense!
    UINT dim = local_dims.size();

    // Set up file access property list with parallel I/O access
    hid_t plist_id= H5Pcreate( H5P_FILE_ACCESS );

    // Set up file access property list with parallel I/O access
    status=H5Pset_fapl_mpio(plist_id, communicator, mpiInfo);  //collective MPI!!! needs to be called on each processor of the communicator
    assert(plist_id>-1);

//    // Create a new file using default properties.
//    hid_t file_id= H5Fcreate(
//                             data_filename.c_str()    // name of the file to be created
//                             , H5F_ACC_TRUNC          // if you are trying to create a file that exists already, the existing file will be truncated, i.e., all data stored on the original file will be erased
//                             , H5P_DEFAULT            // default file creation property list
//                             , plist_id
//                              );
//    assert( file_id > -1 );
//    H5Pclose(plist_id);

    // Open file with r/w access
    hid_t file_id = H5Fopen(
      filename.c_str(),     // name of the file to be opened
      H5F_ACC_RDWR,         // open with read/write access
      plist_id);         // default file access property list

    assert( file_id > -1 );
    H5Pclose(plist_id);

    /* Check if the dataset already exists */
    // Assume intermediate groups have been created before, check full path to the actual dataset!
    htri_t exists = H5Lexists(file_id, data_name.c_str(), H5P_DEFAULT);
    assert(exists > -1);

    if(exists > 0)
    {
      //debug_jochen << "[HDF5Tools::writeVector_MPI()] Dataset '" << data_name
      //    << "' already exists, aborting write process." << std::endl;

      // TODO Check if existing data (x,y vectors) match the data to be written
      // (otherwise GFs with non-matching (x,y) vectors might be added here)

      // Don't forget to close the file!
      status = H5Fclose( file_id );
      assert( status > -1 );

      return;
    }

    // (No compression settings applied here. Parallel HDF5 does not support compression.)

    // set the global size of the grid into a vector of type hsize_t (needed for HDF5 routines)
    hsize_t global_dim_HDF5[ dim];
    for(UINT i=0; i<dim;i++)
    {
      //global_dim_HDF5[i] = global_dim[i];
      global_dim_HDF5[dim-i-1]=global_dims[i];

      //debug_jochen << "global_dim_HDF5[" << (dim-i-1) << "] = " << global_dim_HDF5[dim-i-1] << std::endl;
    }

    // set the count and offset in the different dimensions (determine the size of the hyperslab)
    // (in hsize_t format, needed for HDF5 routines)
    hsize_t count[dim], offset[dim];
    for(UINT i=0;i<dim;i++)
    {
      //count[i]=local_count[i];
      //offset[i]=local_offset[i];
      count[dim-i-1]=local_dims[i];
      offset[dim-i-1]=offsets[i];

      //debug_jochen << "count[" << (dim-i-1) << "] = " << count[dim-i-1] << std::endl;
      //debug_jochen << "offset[" << (dim-i-1) << "] = " << offset[dim-i-1] << std::endl;
    }

    //define the total size of the local data
    hsize_t nAllLocalCells = count[0] * count[1];

    //  Create the dataspace for the dataset.
    hid_t filespace = H5Screate_simple(dim, global_dim_HDF5, NULL);
    assert(filespace>-1);

    // Create the dataset with default properties and close filespace.
    hid_t dset_id = H5Dcreate(file_id, data_name.c_str(), H5T_NATIVE_DOUBLE, filespace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Sclose(filespace);
    assert(dset_id>-1);

    //get the memoryspace (but only if something needs to be written on this processor!)
    hid_t memspace_id = -1;
    if(nAllLocalCells!=0) // -> otherwise HDF5 warning, because of writing nothing!
    {
      memspace_id = H5Screate_simple(dim, count, NULL);
      assert(memspace_id>-1);
    }

    // Select hyperslab in the file.
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    //logger<< "write_parallel_to_HDF5_without_DUNE:  hyperslab selected!" << std::endl;

    // Create property list for collective dataset write.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); // all MPI processors of this communicator need to call this (if they write or not)

    // finally write the data to the disk
    // even if nothing should be written H5Dwrite needs to be called!!
    if(nAllLocalCells!=0) // -> otherwise HDF5 warning, because of writing nothing!
    {
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace,
                        plist_id, &( data[0] )  );
      assert(status>-1);
    } else{ // IMPORTANT. otherwise the H5Dwrite() blocks!!!
      status = H5Dwrite(dset_id,
                        H5T_NATIVE_DOUBLE,
                        H5S_ALL,
                        filespace,
                        plist_id,
                        &( data[0] )  );
      assert(status>-1);
    }

    // Close the property list;
    status=H5Pclose(plist_id);
    assert( status > -1 );

    // Close the filespace;
    status=H5Sclose(filespace);
    assert( status > -1 );

    //if something writen close the memspace
    if(nAllLocalCells!=0)
    {
      // Close the mem space;
      status=H5Sclose(memspace_id);
      assert( status > -1 );
    }

    // Close the dataset;
    status=H5Dclose(dset_id);
    assert( status > -1 );

    // Close the file.
    status = H5Fclose( file_id );
    assert( status > -1 );

    //propably not needed. because the H5Dwrite blocks anyway!!
    //MPI_Barrier(communicator);
  }

};

#endif /* DUNE_AX1_HDF5_TOOLS_HH */
