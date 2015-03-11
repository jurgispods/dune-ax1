/* 
 * File:   HDF5Tools.hh
 * Author: ngo
 *
 * Created on Jan 26, 2013
 */

#ifndef DUNE_GEOINVERSION_HDF5_TOOLS_HH
#define	DUNE_GEOINVERSION_HDF5_TOOLS_HH

#include <dune/pdelab/common/geometrywrapper.hh>
#include <assert.h>
#include <sstream>

#include "hdf5.h"

#define USE_HDF5

namespace Dune {
  namespace GeoInversion {


    class HDF5Tools{
      
    private:
      HDF5Tools(){};

    public:

#ifdef USE_HDF5
      /*
       * Note that in hdf5, all array indices are ordered the other way round!
       *
       * So far, we had this:
       *
       *   nCells[0] = number of cells in x-direction
       *   nCells[1] = number of cells in y-direction
       *   nCells[2] = number of cells in z-direction
       *
       * But storing data in hdf5 requires to define a dimensions vector 'dim' with 
       *
       *  dims[2] = (hsize_t) ( nCells[0] ); // x-direction
       *  dims[1] = (hsize_t) ( nCells[1] ); // y-direction
       *  dims[0] = (hsize_t) ( nCells[2] ); // z-direction
       * 
       * for 3d, and 
       *
       *
       * dims[1] = (hsize_t) ( nCells[0] ); // x-direction
       * dims[0] = (hsize_t) ( nCells[1] ); // y-direction
       * 
       * for 2d.
       *
       */


      /*
       * \tparam dims = numbers of cells in each dimension
       *
       *
       * See
       * http://www.xdmf.org/index.php/Main_Page
       * for more information on this xml format!
       */
      static void writeMetaFileForHDF5( const std::string& meta_filename
                                 , const std::string& data_filename
                                 , const std::string& data_name
                                 , const hsize_t* dims
                                 , const std::string& meta_title
                                 )
      {

        // save header:
#ifdef DIMENSION3
        const int dim = 3;
        hsize_t col = dims[2]; // x-direction
        hsize_t row = dims[1]; // y-direction
        hsize_t dep = dims[0]; // z-direction
#else
        const int dim = 2;
        hsize_t col = dims[1]; // x-direction
        hsize_t row = dims[0]; // y-direction
#endif

        std::ofstream outfile( meta_filename.c_str(), std::ios::out );
        if( outfile.is_open() )
          {

            // Strip filename from path:
            size_t found=data_filename.find_last_of("/\\");
	  
            outfile << "<?xml version=\"1.0\" ?>" << std::endl;
            outfile << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
            outfile << "<Xdmf Version=\"2.0\">" << std::endl;
            outfile << "  <Domain>" << std::endl;
            outfile << "    <Grid Name=\"" << meta_title.c_str() <<  "\" GridType=\"Uniform\">" << std::endl << std::endl;

#ifdef DIMENSION3
            outfile << "      <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\""
              //<< col+1 << " " << row+1 << " " << dep+1 << "\"/>" 
                    << dep+1 << " " << row+1 << " " << col+1 << "\"/>" 
                    << std::endl << std::endl;
#else
            // Take care: 2d-case must be imbedded into 3d-box with height=1 in order to be displayed properly inside paraview!
            outfile << "      <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"2 "
                    << row+1 << " " << col+1 << "\"/>" 
                    << std::endl << std::endl;
#endif


#ifdef _EXTRA_DATA_
#ifdef DIMENSION3
            outfile << "      <Geometry GeometryType=\"ORIGIN_DXDYDZ\">" << std::endl;
#else
            outfile << "      <Geometry GeometryType=\"ORIGIN_DXDYDZ\">" << std::endl;
#endif

            outfile << "        <DataItem ItemType=\"Uniform\" Format=\"HDF\" NumberType=\"Double\" Precision=\"4\" Dimensions=\"" << dim << "\">" << std::endl;
            outfile << data_filename.substr(found+1) << ":/origin" << std::endl;
            outfile << "        </DataItem>" << std::endl;

            outfile << "        <DataItem ItemType=\"Uniform\" Format=\"HDF\" NumberType=\"Double\" Precision=\"4\" Dimensions=\"" << dim << "\">" << std::endl;
            outfile << data_filename.substr(found+1) << ":/resolution" << std::endl;
            outfile << "        </DataItem>" << std::endl;

            outfile << "      </Geometry>" << std::endl;

#endif //_EXTRA_DATA_


            outfile << "      <Attribute Name=\"" << meta_title.c_str() << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;

#ifdef DIMENSION3
            outfile << "        <DataItem ItemType=\"Uniform\" Format=\"HDF\" NumberType=\"Double\" Precision=\"4\" Dimensions=\"" 
                    << dep << " " << row << " " << col 
                    << "\">" << std::endl;
#else
            // Take care: 2d-case must be imbedded into 3d-box with height=1 in order to be displayed properly inside paraview!
            outfile << "        <DataItem ItemType=\"Uniform\" Format=\"HDF\" NumberType=\"Double\" Precision=\"4\" Dimensions=\"1 " 
                    << row << " " << col 
                    << "\">" << std::endl;
#endif
            outfile << data_filename.substr(found+1) << ":" << data_name.c_str() << std::endl;

            outfile << "        </DataItem>" << std::endl;
            outfile << "      </Attribute>" << std::endl;
            outfile << "    </Grid>" << std::endl;
            outfile << "  </Domain>" << std::endl;
            outfile << "</Xdmf>" << std::endl;
  
            outfile.close();
          }

      }



      /** Function to create hdf5 file. It stores the data of the virtual Y-field vector in spatial format such that hyperslabs of data can be cut out later ( see function read_parallel_from_HDF5() )
       *
       \tparam data = vector that is been stored. (one processor only!!!)
       \tparam data_name = name/location where the data a stored in the hdf5-file
       \tparam nCells = numbers of virtual cells in each dimension
       \tparam nChunks = numbers of chunks (needed for hdf5) in each dimension
       \tparam gridsizes = virtual grid-sizes, needed only for extra data
       \tparam data_filename = name of the hdf5 file where all the data will be stored
       \tparam meta_title = display title for meta-file
       \tparam meta_filename = name of the meta file
       *
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */
      static void write_sequential_to_HDF5( 
                                    const Vector<REAL> &data
                                    , const std::string& data_name
                                    , const Vector<UINT> &nCells
                                    , const Vector<UINT> &nChunks
                                    , const Vector<CTYPE> &gridsizes
                                    , const std::string& data_filename
                                    , const std::string& meta_title="default"
#ifdef DIMENSION3
                                    , const std::string& meta_filename="BUFFER/3D/default.xmf"
#else
                                    , const std::string& meta_filename="BUFFER/2D/default.xmf"				   
#endif
                                     )
      {
        UINT dim = nCells.size();

        //std::cout << std::endl << "write_sequential_to_HDF5: dim = " << dim << std::endl;
        logger << "write_sequential_to_HDF5: data_name = " << data_name << std::endl;

        /* Create a new file using default properties. */
        hid_t file_id = H5Fcreate( 
                                  data_filename.c_str()    // name of the file to be created
                                  , H5F_ACC_TRUNC          // if you are trying to create a file that exists already, the existing file will be truncated, i.e., all data stored on the original file will be erased
                                  , H5P_DEFAULT            // default file creation property list
                                  , H5P_DEFAULT            // default file access property list
                                   );

        assert( file_id > -1 );
        //logger << "write_sequential_to_HDF5: h5-file created: " << data_filename.c_str() << std::endl;


        hsize_t mdims[1];
        mdims[0] = data.size();
        hid_t memspace_id = H5Screate_simple(           // H5Screate_simple creates a new simple dataspace and opens it for access, returning a dataset identifier.
                                             1          // rank=1 is the number of dimensions used in the dataspace
                                             , mdims    // mdims is a one-dimensional array of size rank specifying the size of each dimension of the dataset.
                                             , NULL     // maxdims is an array of the same size specifying the upper limit on the size of each dimension. maxdims may be NULL, in which case the maxdims == mdims. 
                                                        ); // must be released with H5Sclose or resource leaks will occur

        assert( memspace_id > -1 );
        //logger << "write_sequential_to_HDF5: memspace_id created." << std::endl;


        /* Create the dataspace for the dataset.
         * The dataspace describes the dimensions of the dataset array. */

        hsize_t mindims = 50000;

        hsize_t dims[ dim ];
        
        for(UINT i=0;i<dim;i++){
          dims[dim-i-1] = (hsize_t) ( nCells[i] );
          mindims = std::min( mindims, dims[dim-i-1] );
        }


        //logger << "dims [0] = "<< dims[0]<<std::endl;
        //logger << "dims [1] = "<< dims[1]<<std::endl;
        hid_t dataspace_id = H5Screate_simple( 
                                              dim     // number of dimensions = rank
                                              , dims  // vector containing sizes per dimension
                                              , NULL  // maxdims == dims
                                               );
        assert( dataspace_id > -1 );
        //logger << "write_sequential_to_HDF5: dataspace_id created." << std::endl;



        hid_t plist_id = H5Pcreate(                    // The new property list is initialized with default values for the specified class.
                                   H5P_DATASET_CREATE  // Properties for dataset creation
                                                       );
        assert( plist_id > -1 );
        //logger << "write_sequential_to_HDF5: plist_id created." << std::endl;


        herr_t status;
        
        hsize_t maxnchunks = 1;
        for(UINT i=0;i<dim;i++){
          maxnchunks = std::max( maxnchunks, (hsize_t) nChunks[i] );
        }


        hsize_t chunk_dims[ dim ];
        if( maxnchunks >= mindims ){
          logger << "Warning: maxnchunks > mindims, set chunk_dims[i]=0" << std::endl;
          for(UINT i=0;i<dim;i++)
            chunk_dims[i]=(hsize_t)1;
        }
        else{
          for(UINT i=0;i<dim;i++){
            chunk_dims[dim-i-1] = (hsize_t) nChunks[i];
          }
        }

        

        status = H5Pset_chunk( 
                              plist_id
                              , dim             // must be == rank of the dataset
                              , chunk_dims      // The values of the check_dims array define the size of the chunks to store the dataset's raw data. The unit of measure for check_dims values is dataset elements. 
                               );

        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Pset_chunk() o.k." << std::endl;


        status = H5Pset_shuffle( plist_id ); // Sets the shuffle filter, H5Z_FILTER_SHUFFLE, in the dataset creation property list. This re-orders data to simplify compression.
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Pset_shuffle() o.k." << std::endl;
  

        status = H5Pset_deflate( plist_id, 1 ); // Sets deflate (GNU gzip) compression method and compression level. ( 0 < level < 9, lower = faster, but less compression )
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Pset_deflate() o.k." << std::endl;



        /* Create the dataset. */
        hid_t dataset_id = H5Dcreate( 
                                     file_id             // Location identifier: id of the file or the group within which to create the dataset
                                     , data_name.c_str()  // Dataset name: may be either an absolute path in the file or a relative path from file_id naming the dataset. 
                                     , HDF5_DATA_TYPE   // Datatype identifier, here: IEEE floating point 32-bit little-endian
                                     , dataspace_id      // Dataspace identifier 
                                     , plist_id          // Dataset creation property list identifier
                                      );
        assert( dataset_id > -1 );
        //logger << "write_sequential_to_HDF5: dataset_id created." << std::endl;


        /* Write the dataset. */
        status = H5Dwrite(                         // Writes raw data from a buffer to a dataset.
                          dataset_id               // dataset identifier
                          , H5T_NATIVE_DOUBLE      // memory datatype id
                          , memspace_id            // specifies the memory dataspace and the selection within it
                          , H5S_ALL                // specifies the selection within the file dataset's dataspace. H5S_ALL indicates that the entire file dataspace, as defined by the current dimensions of the dataset, is to be selected
                          , H5P_DEFAULT            // Identifier of a transfer property list for this I/O operation. H5P_DEFAULT: The default data transfer properties are used.
                          , &(data[0])    // application memory buffer
                                                   );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Dwrite( " <<data_name.c_str()<<" ) o.k." << std::endl;

        status = H5Sclose( dataspace_id );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: dataspace_id closed." << std::endl;

        status = H5Sclose( memspace_id );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: memspace_id closed." << std::endl;

        status = H5Dclose( dataset_id );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: dataset_id closed." << std::endl;



#ifdef _EXTRA_DATA_

        /* Create the data space for the attribute. */
        hsize_t adims = dim;
        hid_t attribute_dataspace_id = H5Screate_simple( 
                                                        1
                                                        , &adims
                                                        , NULL 
                                                         );


        /* Create the dataset for the location of the origin (reference point) */
        CTYPE attr_data[ 3 ]; // Take care: Here, we use 3d data even for the 2d case! Reason: Data mismatch in paraview + xmf
        attr_data[0]  = 0.0; 
        attr_data[1]  = 0.0; 
        attr_data[2]  = 0.0;

        /* Write the dataset. */
        hid_t attribute_id = H5Dcreate( 
                                       file_id
                                       , "/origin"
                                       , HDF5_DATA_TYPE
                                       , attribute_dataspace_id
                                       , H5P_DEFAULT 
                                        );
        assert( attribute_id > -1 );
        //logger << "write_sequential_to_HDF5: H5Dcreate(/origin) done." << std::endl;


        status = H5Dwrite( 
                          attribute_id
                          , H5T_NATIVE_DOUBLE
                          , H5S_ALL
                          , H5S_ALL
                          , H5P_DEFAULT
                          , attr_data 
                           );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Dwrite(/origin) done." << std::endl;


        /* Close the attribute for the origin. */
        status = H5Dclose( attribute_id );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Dclose(/origin) done." << std::endl;


        /* Create the dataset for the resolution */
        attr_data[0]= gridsizes[0];
        attr_data[1]= gridsizes[1];
        attr_data[2]= gridsizes[2];

        /* Write the dataset. */
        attribute_id = H5Dcreate(
                                 file_id
                                 , "/resolution"
                                 , HDF5_DATA_TYPE
                                 , attribute_dataspace_id
                                 , H5P_DEFAULT
                                 );
        assert( attribute_id > -1 );
        //logger << "write_sequential_to_HDF5: H5Dcreate(/resolution) done." << std::endl;
  

        status = H5Dwrite( 
                          attribute_id
                          , H5T_NATIVE_DOUBLE
                          , H5S_ALL
                          , H5S_ALL
                          , H5P_DEFAULT
                          , attr_data 
                           );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Dwrite(/resolution) done." << std::endl;

        status = H5Sclose( attribute_dataspace_id );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: attribute_dataspace_id closed." << std::endl;

        /* Close the attribute for the resolution. */
        status = H5Dclose( attribute_id );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Dclose(/resolution) done." << std::endl;

#endif //_EXTRA_DATA_

        /* Close the property list. */
        status = H5Pclose( plist_id );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Pclose(plist_id) done." << std::endl;

        /* Close the file. */
        status = H5Fclose( file_id );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Fclose( file_id ) done." << std::endl;

        /* Write meta file for display in paraview */
        writeMetaFileForHDF5( meta_filename, data_filename, data_name, &(dims[0]), meta_title  );

      };

      /** function to write data (in parallel) to the hdf5 file
       *
       * \tparam gv is the gridview
       * \tparam inputdata is the inputdata object
       * \tparam data are the data which are written into the file (current hyperslab)
       * \tparam data_name is the name/path where the data are stored in the HDF5 file
       * \tparam nCells  the number of virtual(global) cells for each dimension
       * \tparam nChunks = numbers of chunks (needed for hdf5) in each dimension
       * \tparam gridsizes = virtual grid-sizes, needed only for extra data (not used so far)
       * \tparam data_filename = name of the hdf5 file where all the data will be stored
       * \tparam writemetafile = flag indicating if a meta file should be written
       * \tparam meta_title = display title for meta-file
       * \tparam meta_filename = name of the meta file
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */


      template<typename GV>
      static void write_parallel_to_HDF5(
                                         const GV& gv
                                         , const CInputData& inputdata
                                         , const Vector<REAL> &data
                                         , const std::string& data_name
                                         , const Vector<UINT> &nCells   // nKnotsPerDim
                                         , const Vector<UINT> &nChunks
                                         , const Vector<CTYPE> &gridsizes
                                         , const std::string& data_filename
                                         , const int blocksizePerDimension = 1
                                         , const bool writemetafile = false
                                         , const FEMType::Type femtype=FEMType::DG
                                         , const std::string& meta_title="default"
#ifdef DIMENSION3
                                         , const std::string& meta_filename="BUFFER/3D/default.xmf"
#else
                                         , const std::string& meta_filename="BUFFER/2D/default.xmf"				   
#endif
                                         )
      {
        logger << "write_parallel_to_HDF5: ... " << data_filename << std::endl;
        
        //get the communicator of the dune grid
        const Dune::MPIHelper::MPICommunicator &comm = gv.comm();
  
        //Info varibale need for the HDF5
        MPI_Info mpiInfo = MPI_INFO_NULL;
  
        //some needed variables and typdefs
        const UINT dim = GV::dimension;  // dimensionality of the problem (2-d or 3-d)

        typedef typename GV::template Codim<0>::template Partition<Dune::All_Partition>::Iterator LeafIterator;

        typedef typename GV::Grid GRIDTYPE;

        // make a mapper for codim 0 entities in the leaf grid
        Dune::LeafMultipleCodimMultipleGeomTypeMapper<GRIDTYPE,P0Layout> mapper(gv.grid()); // get the underlying hierarchical grid out ot the gridview

        // variable fot the size of the local data set.
        //hsize_t mappersize = mapper.size();
        //logger << "mapper-size=" << mappersize << std::endl;
  
        // get indices
        //const typename GV::IndexSet& indexset = gv.indexSet();
  
        //variables for defining the offset and the count in the different dimensions
        std::vector<UINT> max_index( dim, 0 );
        std::vector<UINT> min_index( dim, 100000 );

        /*
          1.) Loop over the current hyperslab (hyperslab = all the leaf grid-cells that belong to the current processor) and get the global coordinate 'x' of each grid-cell center.
          Remark: From now on, we are working solemnly with the virtual Y-field grid, using virtual indices. Later, in the evaluate2d() and evaluate3d() member functions of the class kfieldgenerator, we are doing the very same conversion from a global coordinate to the virtual global index. That's why it can work independently of the grid type!
          2.) Translate this global coordinate 'x' into its corresponding Yfield-tensor-index 'global_index'. This requires only the virtual gridsize of the virtual Yfield grid. So 'global_index' is just a virtual index.
          3.) 'min_index' corresponds to the virtual cell with the lowest coordinate values in the current hyperslab. It is the offset (= translational displacement) from the very first virtual cell (0,0).
          4.) 'local_count[i]' is the number of virutal cells in the dimension 'i' of the current hyperslab. 
          Remark: The volume of the current hyperslab is 'nAllLocalCells'. 'mappersize' == 'nAllLocalCells' only if the resolution of the yaspgrid == resolution of the virtual grid !

        */

        for (LeafIterator it = gv.template begin<0,Dune::All_Partition> ()
               ; it != gv.template end<0,Dune::All_Partition> ()
               ; ++it) 
          {
            //int id = indexset.index(*it);
            //logger << "map index = " << id << std::endl;

            //Dune::FieldVector<CTYPE, dim> currentpoint(0.0);	  
            // Get the local coordinate:
            Dune::FieldVector<CTYPE,dim> x = it->geometry().center();
            //Dune::FieldVector<CTYPE,dim> currentpoint_local = it->geometry().local( x );

	  
            std::vector<UINT> global_index(dim,0);
            for(UINT i=0; i<dim; i++ )
              {
                global_index[i] = (UINT) floor( x[i] / inputdata.domain_data.virtual_gridsizes[i] );
                //broken_global_index[i] = global_index[i] / ifactor;
                max_index[i] = std::max( global_index[i], max_index[i] );
                min_index[i] = std::min( global_index[i], min_index[i] );
              }
#ifdef DEBUG_LOG_LEVEL_2
            logger << "global-index = ("  << global_index[0] << "," << global_index[1] 
#ifdef DIMENSION3
                   << "," << global_index[2] 
#endif 
                   << ")" << std::endl;
#endif //DEBUG_LOG_LEVEL_2
            // Idee für UG später: 
            // lese direkt punktweise aus der h5-Datei, über den offset[] = global_index[], mit count = (1 1 1) ?
          } // end of loop over leaf elements
  
        // Set up file access property list with parallel I/O access
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, comm, mpiInfo);
        assert(plist_id>-1);
        //logger << "write_parallel_to_HDF5: plist_id" << std::endl;
  
        //  Create a new file collectively and release property list identifier.
        hid_t file_id = H5Fcreate(data_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);
        assert(file_id>-1);
        //logger << "write_parallel_to_HDF5: "<<data_filename.c_str()<<" file created!" << std::endl;
  
        // set the global size of the grid into a vector of type hsize_t (needed for HDF5 routines)
        hsize_t dims_global[ dim];
  
        for(UINT i=0; i<dim;i++)
          {
            dims_global[dim -i-1]=nCells[i];
            //dims_global[i]=nCells[i];
            //logger << "global dims[ "<<dim -i-1<<" ] = "<<dims_global[dim -i-1]<<"."<<std::endl;
          }
    
        //  Create the dataspace for the dataset.
        hid_t filespace = H5Screate_simple(dim, dims_global, NULL);
        assert(filespace>-1);
        //logger << "write_parallel_to_HDF5:  fileSPACE created!" << std::endl;
  
  
        // Create the dataset with default properties and close filespace.
        hid_t dset_id = H5Dcreate(file_id, data_name.c_str(), HDF5_DATA_TYPE, filespace,
                                  H5P_DEFAULT);
        H5Sclose(filespace);
        assert(dset_id>-1);
        //logger<< "write_parallel_to_HDF5:  dataset created!" << std::endl;
  
  
        // set the count in the different dimensions (determine the size of the hyperslab)
        hsize_t count[ dim ];


        for(UINT i=0; i<dim; i++ )
          {
            //logger << "write_parallel_to_HDF5: nCells[" << i << "] = " << nCells[i] << std::endl;
            //logger << "write_parallel_to_HDF5: max_index[" << i << "] = " << max_index[i] << std::endl;
            //logger << "write_parallel_to_HDF5: min_index[" << i << "] = " << min_index[i] << std::endl;

#ifdef USE_FEM
            // be careful!!! needed when vertex-based data instead of cell based data
            count[dim-i-1] = max_index[i] - min_index[i] + 1 + ( nCells[i] - inputdata.domain_data.nCells[i] );
            //count[i] = max_index[i] - min_index[i] + 1;
#else
            if(femtype==FEMType::CG){
              count[dim-i-1] = max_index[i] - min_index[i] + 1 + ( nCells[i] - inputdata.domain_data.nCells[i] );
            }
            else{ 
              count[dim-i-1] = max_index[i] - min_index[i] + 1;
              count[dim-i-1] *= blocksizePerDimension;
            }
#endif
            //logger << "write_parallel_to_HDF5: count[" << dim-i-1 << "] = " << count[dim-i-1] << std::endl;
          }

        //define the total size of the local data
  
        //hsize_t nAllLocalCells = count[0] * count[1];
#ifdef DIMENSION3
        //nAllLocalCells *= count[2];
#endif
  
        //set the offset of the data!
        hsize_t offset[ dim ];
        for(UINT i=0; i<dim; i++ )
          {
#ifdef USE_FEM
            offset[dim-i-1] = min_index[i];
#else
            if(femtype==FEMType::CG)
              offset[dim-i-1] = min_index[i];
            else
              offset[dim-i-1] = min_index[i] * blocksizePerDimension;
#endif
            //offset[i] = min_index[i];
            //logger << "offset[" << dim-i-1 << "] = " << offset[dim-i-1] << std::endl;
          }
  
        // Each process defines dataset in memory and writes it to the hyperslab in the file.
        hid_t memspace_id = H5Screate_simple(dim, count, NULL);
        assert(memspace_id>-1);
        //logger<< "write_parallel_to_HDF5:  memspace created!" << std::endl;
  
        // Select hyperslab in the file.
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
        //logger<< "write_parallel_to_HDF5:  hyperslab selected!" << std::endl;
  
        // Create property list for collective dataset write.
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        //logger<< "write_parallel_to_HDF5:  properties set!" << std::endl;
  
        // fineally write the data to the disk
        //logger<< "write_parallel_to_HDF5:  writing ... " << std::endl;
        herr_t status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace,
                                 plist_id, &( data[0] )  );
        assert(status>-1);
        //logger<< "write_parallel_to_HDF5:  .... done (writing) " << std::endl;
  
#ifdef _EXTRA_DATA_

        /* Create the data space for the attribute. */
        hsize_t adims = dim;
        hid_t attribute_dataspace_id = H5Screate_simple( 
                                                        1
                                                        , &adims
                                                        , NULL 
                                                         );


        /* Create the dataset for the location of the origin (reference point) */
        CTYPE attr_data[ 3 ]; // Take care: Here, we use 3d data even for the 2d case! Reason: Data mismatch in paraview + xmf
        attr_data[0]  = 0.0; 
        attr_data[1]  = 0.0; 
        attr_data[2]  = 0.0;

        /* Write the dataset. */
        hid_t attribute_id = H5Dcreate( 
                                       file_id
                                       , "/origin"
                                       , HDF5_DATA_TYPE
                                       , attribute_dataspace_id
                                       , H5P_DEFAULT 
                                        );
        assert( attribute_id > -1 );
        //logger << "write_sequential_to_HDF5: H5Dcreate(/origin) done." << std::endl;


        status = H5Dwrite( 
                          attribute_id
                          , H5T_NATIVE_DOUBLE
                          , H5S_ALL
                          , H5S_ALL
                          , H5P_DEFAULT
                          , attr_data 
                           );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Dwrite(/origin) done." << std::endl;


        /* Close the attribute for the origin. */
        status = H5Dclose( attribute_id );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Dclose(/origin) done." << std::endl;


        /* Create the dataset for the resolution */
        attr_data[0]= gridsizes[0];
        attr_data[1]= gridsizes[1];
        attr_data[2]= gridsizes[2];

        /* Write the dataset. */
        attribute_id = H5Dcreate(
                                 file_id
                                 , "/resolution"
                                 , HDF5_DATA_TYPE
                                 , attribute_dataspace_id
                                 , H5P_DEFAULT
                                 );
        assert( attribute_id > -1 );
        //logger << "write_sequential_to_HDF5: H5Dcreate(/resolution) done." << std::endl;
  

        status = H5Dwrite( 
                          attribute_id
                          , H5T_NATIVE_DOUBLE
                          , H5S_ALL
                          , H5S_ALL
                          , H5P_DEFAULT
                          , attr_data 
                           );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Dwrite(/resolution) done." << std::endl;

        status = H5Sclose( attribute_dataspace_id );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: attribute_dataspace_id closed." << std::endl;

        /* Close the attribute for the resolution. */
        status = H5Dclose( attribute_id );
        assert( status > -1 );
        //logger << "write_sequential_to_HDF5: H5Dclose(/resolution) done." << std::endl;

#endif //_EXTRA_DATA_
  
        // close the used identifyers
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace_id);
        H5Pclose(plist_id);
        H5Fclose(file_id);
  
        // write meta file if needed
        if(writemetafile)
          /* Write meta file for display in paraview */
          writeMetaFileForHDF5( meta_filename, data_filename, data_name, &(dims_global[0]), meta_title  );
  
      };

      /** function to retrieve data (sequential) from the hdf5 file
       *
       * \tparam local_data is a 'return value' which gets the data that belongs to the current processor (current hyperslab)
       * \tparam data_name is the name/path where the data are stored in the HDF5 file
       * \tparam local_count is a 'return value' which gets the number of virtual cells for each dimension
       * \tparam local_offset is a 'return value' which gets the distance of the current hyperslab from the origin.
       * \tparam data_filename is the filename of the data file
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */





      static void read_sequential_from_HDF5(
                                            Vector<REAL>& local_data
                                            , const std::string& data_name
                                            , Vector<UINT>& local_count
                                            , Vector<UINT>& local_offset
                                            , const std::string& data_filename
                                            , const int & iPart=1
                                            , const int & nParts=1
                                            )
      {
 
  
        // open the file for reading
        hid_t file_id=  H5Fopen (data_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);  
        assert( file_id > -1 );
        //logger << "read_sequential_from_HDF5(): h5-file open: " << data_filename.c_str() << std::endl;
  
        // open the dataset
        hid_t dataset_id = H5Dopen(file_id, data_name.c_str());
        assert( dataset_id > -1 );
        //logger << "read_sequential_from_HDF5(): H5Dopen dataset done." << std::endl;
    
        // gwt the dataspace
        hid_t dataspace_id=H5Dget_space (dataset_id); 
        assert(dataspace_id>-1);
        //logger << "read_sequential_from_HDF5(): H5Dget_space  done." << std::endl;

        // some needed variables 
        herr_t status;
        hsize_t dim,*dims;
  
        // get the dimension (2-d or 3-d)
        //logger << "Getting parameters for reading from data (read all)"<<std::endl; 
        dim = H5Sget_simple_extent_ndims( dataspace_id );
        //logger << "Got dim = " << dim << std::endl;
  
        //logger << "read_sequential_from_HDF5()! ( dim = "<<dim<< " )" << std::endl;
   
        // get the size of the problem
        dims=(hsize_t*)malloc(dim * sizeof (hsize_t));
        status = H5Sget_simple_extent_dims( dataspace_id , dims, 0);
        assert( status > -1 );
 
        // file the return value. the local_count
        local_count.resize(dim);
        //logger << "Got local_count from dataspace:" << std::endl;
        for( UINT i=0; i<dim; i++ )
          {
            local_count[dim-1-i]=dims[i];
            //logger << "local_count[" << dim-i-1 << "] = " << local_count[dim-1-i] << std::endl; 
      
          }
    
        // set offset(for reading) and offset_local(as return value) to 0 (read all the data, maybe later also used to read only parts, maybe we can kick it out first)
        // like this it is not nice but works, maybe we change it later!!! 
        // local_size is the size of the data
        hsize_t offset[dim];
        local_offset.resize(dim);
        hsize_t  local_size=1;
        for( UINT i=0; i<dim; i++ ){
          local_offset[i]=0;
          offset[i]=0;
        }

        if(nParts>1){
          local_count[1]/=nParts;
          local_offset[1]=local_count[1]*iPart;
          offset[1]=local_offset[1];
          if(iPart+1==nParts){
            local_count[1]=dims[1]-local_offset[1];
          }
          for( UINT i=0; i<dim; i++ ){
            dims[i]=local_count[dim-1-i];
          }
        }
        
        //logger << "local_count[" << 1 << "] = " << local_count[1] << std::endl; 
        hsize_t stride[ dim ];
        for( UINT i=0; i<dim; i++ ){
          local_size*=local_count[i];
          stride[i]=1;
        }
        // create the memory space
        hid_t memspace_id = H5Screate_simple (dim, dims, NULL);
        //hid_t memspace_id = H5Screate_simple (1, &local_size, NULL);
    
        //select the hyperslab-
        status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, stride,
                                      dims, NULL);
        assert(status>-1);
        //logger << "read_sequential_from_HDF5(): hyperslab selected." << std::endl;

        //resize the return data
        local_data.resize( local_size /*mappersize*/ );
  
        /* set up the collective transfer properties list */
        hid_t xferPropList = H5Pcreate (H5P_DATASET_XFER);
        assert( xferPropList > -1 );
        //logger << "xferPropList created." << std::endl;

        local_data.resize(local_size);
        // finally the reading from the file
        status = H5Dread(                     dataset_id
                                              , H5T_NATIVE_DOUBLE //image.DataType
                                              , memspace_id
                                              , dataspace_id // H5S_ALL // 
                                              , xferPropList //H5P_DEFAULT // 
                                              , &( local_data[0] )
                                              );
        assert(status>-1);
        //logger << "read_sequential_from_HDF5(): reading data done." << std::endl;

        // close the identifyers
        H5Dclose (dataset_id);
        H5Sclose (dataspace_id);
        H5Sclose (memspace_id);
        free(dims);
        status = H5Fclose( file_id );
        assert( status > -1 );
        //logger << "read_sequential_from_HDF5(): H5Fclose( file_id ) done." << std::endl;
      };










#ifdef PARALLEL
      /** function to retrieve data (in parallel) from the hdf5 file
       *
       * \tparam gv is the gridview
       * \tparam inputdata is the inputdata object
       * \tparam local_data is a 'return value' which gets the data that belongs to the current processor (current hyperslab)
       * \tparam data_name is the name/path where the data are stored in the HDF5 file
       * \tparam local_count is a 'return value' which gets the number of virtual cells for each dimension
       * \tparam local_offset is a 'return value' which gets the distance of the current hyperslab from the origin.
       * \tparam data_filename is the filename of the data file
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */

      template<typename GV>
      static void read_parallel_from_HDF5(
                                          const GV& gv
                                          , const CInputData& inputdata
                                          , Vector<REAL>& local_data
                                          , const std::string& data_name
                                          , Vector<UINT>& local_count
                                          , Vector<UINT>& local_offset
                                          , const std::string& data_filename
                                          , const int blocksizePerDimension=1
                                          , const int baselevel=0
                                          , const FEMType::Type femtype=FEMType::DG
                                          ) 
      {
        logger << "read_parallel_from_HDF5 ... " << data_filename << std::endl;

        const Dune::MPIHelper::MPICommunicator &comm = gv.comm();

        MPI_Info mpiInfo = MPI_INFO_NULL;
        
        const UINT dim = GV::dimension;
        typedef typename GV::template Codim<0>::template Partition<Dune::All_Partition>::Iterator LeafIterator;
        typedef typename GV::Grid GRIDTYPE;
  
        //logger << "read_parallel_from_HDF5()! ( dim = "<<dim<<" )" << std::endl;

        // make a mapper for codim 0 entities in the leaf grid
        // Dune::LeafMultipleCodimMultipleGeomTypeMapper<GRIDTYPE, P0Layout> mapper(gv.grid());
        // make a mapper for codim 0 entities in the base-level grid
        Dune::LevelMultipleCodimMultipleGeomTypeMapper<GRIDTYPE, P0Layout> mapper(gv.grid(),baselevel);

        hsize_t mappersize = mapper.size();
        logger << "read_parallel_from_HDF5(): mapper-size=" << mappersize << std::endl;
        assert( mappersize > 0 );
        if( ! mappersize > 0 ){
            DUNE_THROW(Dune::RangeError, "At least one partition has lost its share of the grid! Maybe there are too many processors for this "
                       "resolution?");
        }
            
        //const typename GV::IndexSet& indexset = gv.indexSet();
  
        std::vector<UINT> max_index( dim, 0 );
        std::vector<UINT> min_index( dim, 100000 );

        /*

          1.) Loop over the current hyperslab (hyperslab = all the leaf grid-cells that belong to the current processor) and get the global coordinate 'x' of each grid-cell center.
          Remark: Form now on, we are working solemnly with the virtual Y-field grid, using virtual indices. Later, in the evaluate2d() and evaluate3d() member functions of the class kfieldgenerator, we are doing the very same conversion from a global coordinate to the virtual global index. That's why it can work independently of the grid type!
          2.) Translate this global coordinate 'x' into its corresponding Yfield-tensor-index 'global_index'. This requires only the virtual gridsize of the virtual Yfield grid. So 'global_index' is just a virtual index.
          3.) 'min_index' corresponds to the virtual cell with the lowest coordinate values in the current hyperslab. It is the offset (= translational displacement) from the very first virtual cell (0,0).
          4.) 'local_count[i]' is the number of virutal cells in the dimension 'i' of the current hyperslab. 
          Remark: The volume of the current hyperslab is 'nAllLocalCells'. 'mappersize' == 'nAllLocalCells' only if the resolution of the yaspgrid == resolution of the virtual grid !

        */

        for (LeafIterator it = gv.template begin<0,Dune::All_Partition> ()
               ; it != gv.template end<0,Dune::All_Partition> ()
               ; ++it) 
          {
            //int id = indexset.index(*it);
            //logger << "map index = " << id << std::endl;

            //Dune::FieldVector<CTYPE, dim> currentpoint(0.0);	  
            // Get the local coordinate:
            Dune::FieldVector<CTYPE, dim> x = it->geometry().center();
            //Dune::FieldVector<CTYPE, dim> currentpoint_local = it->geometry().local( x );

	  
            std::vector<UINT> global_index(dim,0);
            for(UINT i=0; i<dim; i++ )
              {
                global_index[i] = (UINT) floor( x[i] / inputdata.domain_data.virtual_gridsizes[i] );
                //broken_global_index[i] = global_index[i] / ifactor;
                max_index[i] = std::max( global_index[i], max_index[i] );
                min_index[i] = std::min( global_index[i], min_index[i] );
              }
#ifdef DEBUG_LOG_LEVEL_2
            logger << "global-index = ("  << global_index[0] << "," << global_index[1] 
#ifdef DIMENSION3
                   << "," << global_index[2] 
#endif 
                   << ")" << std::endl;
#endif //DEBUG_LOG_LEVEL_2
            // Idee für UG später: 
            // lese direkt punktweise aus der h5-Datei, über den offset[] = global_index[], mit count = (1 1 1) ?
          } // end of loop over leaf elements

  


        /* Do some preparations for accessing the hdf5 file... 
         * This hyperslab selection works only for a structured grid, I guess!
         * For an unstructured grid, a pointwise hyperslab selection should still be possible.
         */
        //logger << "HDF5 file - read access..." << std::endl;

        hsize_t count[ dim ];
        local_count.resize( dim );
        local_offset.resize( dim );

        for(UINT i=0; i<dim; i++ )
          {
            local_count[i] = max_index[i] - min_index[i] + 1;
            //logger << "local_count[" << i << "] = " << local_count[i] << std::endl;
          }

        for(UINT i=0; i<dim; i++ )
          {
            count[i] = (hsize_t) local_count[ dim - i - 1 ];
            //logger << "count[" << i <<"] = " << count[i] << std::endl;
          }

        hsize_t offset[ dim ];
        for(UINT i=0; i<dim; i++ )
          {
#ifdef USE_FEM
            local_offset[i] = min_index[i];
#else
            if(femtype==FEMType::CG)
              local_offset[i] = min_index[i];
            else
              local_offset[i] = min_index[i] * blocksizePerDimension;              
#endif

            //logger << "local_offset[" << i << "] = " << local_offset[i] << std::endl;
          }

        for(UINT i=0; i<dim; i++ )
          {	  
            offset[i] = (hsize_t) local_offset[ dim - i - 1 ];
            //logger << "offset[" << i <<"] = " << offset[i] << std::endl;
          }
	

        /* setup file access template with parallel IO access. */
        hid_t access_pList = H5Pcreate( H5P_FILE_ACCESS );
        assert( access_pList > -1 );
        //logger << "access_pList created." << std::endl;

        herr_t status;             /* Generic return value */
        status = H5Pset_fapl_mpio( access_pList, comm, mpiInfo ); // Take care, this requires the MPI version of hdf5!
        assert( status > -1 );
        //logger << "access_pList set" << std::endl;


        /* open the file collectively */
        hid_t infileId = H5Fopen( data_filename.c_str(), H5F_ACC_RDONLY, access_pList );
        assert( infileId > -1 );
        //logger << "infileId opened." << std::endl;
  
 
  
  

        /* Release file-access template */
        status = H5Pclose( access_pList );
        assert( status > -1 );
        //logger << "access_pList closed." << std::endl;

        hid_t indatasetId = H5Dopen( infileId, data_name.c_str() );
        assert( indatasetId > -1 );
        //logger << "indatasetId opened." << std::endl;

        hid_t indataspaceId = H5Dget_space( indatasetId );
        assert( indataspaceId > -1 );
        //logger << "indataspaceId opened." << std::endl;

        //hsize_t rank = H5Sget_simple_extent_ndims( indataspaceId );
        //logger << "Got rank = " << rank << std::endl;

        hsize_t dims[ dim ];
        status = H5Sget_simple_extent_dims( indataspaceId, dims, 0);
        assert( status > -1 );
        //logger << "Got dims from dataspace:" << std::endl;
  
        /*
         * Be careful!!!
         * correction if vertex based data a read!
         */
#ifdef USE_FEM
        for(UINT ii=0;ii<dim;ii++){
          if(dims[dim-1-ii]-inputdata.domain_data.nCells[ii]){
            local_count[ii]+=1;
            count[dim-1-ii]+=1; 
          }
        }
#else
        if(femtype==FEMType::CG) {
          for(UINT ii=0;ii<dim;ii++){

            //logger << "dims["<<dim-1-ii<<"] = " << dims[dim-1-ii] << std::endl;
            //logger << "local_count["<<ii<<"] = " << local_count[ii] << std::endl;
            //logger << "count["<<dim-1-ii<<"] = " << count[dim-1-ii] << std::endl;

            //if(dims[dim-1-ii]-inputdata.domain_data.nCells[ii]){
              local_count[ii]+=1;
              count[dim-1-ii]+=1; 
            //}

            //logger << "local_count["<<ii<<"] = " << local_count[ii] << std::endl;
            //logger << "count["<<dim-1-ii<<"] = " << count[dim-1-ii] << std::endl;
          }
        }
        else{
          for(UINT ii=0;ii<dim;ii++){

            //logger << "dims["<<dim-1-ii<<"] = " << dims[dim-1-ii] << std::endl;

            //if(dims[dim-1-ii]-inputdata.domain_data.nCells[ii]){
              local_count[ii] = blocksizePerDimension * local_count[ii];
              count[dim-1-ii] = blocksizePerDimension * count[dim-1-ii];
            //}

            //logger << "local_count["<<ii<<"] = " << local_count[ii] << std::endl;
            //logger << "count["<<dim-1-ii<<"] = " << count[dim-1-ii] << std::endl;
          }
        }
#endif


        hsize_t nAllLocalCells = local_count[0] * local_count[1];
#ifdef DIMENSION3
        nAllLocalCells *= local_count[2];
#endif
        /*  
            for( UINT i=0; i<dim; i++ )
            {
            logger << "dims[" << i << "] = " << dims[i] << std::endl;
            }
        */

        //hid_t memspaceId = H5Screate_simple(1, &mappersize, NULL);
        hid_t memspaceId = H5Screate_simple(1, &nAllLocalCells, NULL);
        assert(memspaceId > -1 );

        hsize_t stride[ dim ];
        for( size_t i=0; i<dim; ++i )
          stride[i]=1;

        /*
          logger << "stride : " << stride[0] << stride[1] 
          #ifdef DIMENSION3
          << stride[2] 
          #endif
          << std::endl;
          logger << "offset : " << offset[0] << offset[1] 
          #ifdef DIMENSION3
          << offset[2] 
          #endif
          << std::endl;
          logger << "count : " << count[0] << count[1] 
          #ifdef DIMENSION3
          << count[2] 
          #endif
          << std::endl;
        */
        status = H5Sselect_hyperslab( indataspaceId
                                      , H5S_SELECT_SET
                                      , offset
                                      , stride
                                      , count
                                      , NULL
                                      );
        assert( status > -1 );

        //logger << "H5Sselect_hyperslab() done." << std::endl;

        /* set up the collective transfer properties list */
        hid_t xferPropList = H5Pcreate (H5P_DATASET_XFER);
        assert( xferPropList > -1 );
        //logger << "xferPropList created." << std::endl;

        //Vector<double> local_Yfield_vector;
        local_data.resize( nAllLocalCells /*mappersize*/ );

        status = H5Dread( indatasetId
                          , H5T_NATIVE_DOUBLE //image.DataType
                          , memspaceId 
                          , indataspaceId // H5S_ALL // 
                          , xferPropList // H5P_DEFAULT // 
                          , &( local_data[0] )
                          );
        assert( status > -1 );
        //logger << "H5Dread() done." << std::endl;


        status = H5Sclose(indataspaceId);
        assert( status > -1 );
        //logger << "indataspaceId closed." << std::endl;

        status = H5Sclose(memspaceId);
        assert( status > -1 );
        //logger << "memspaceId closed." << std::endl;
  
        status = H5Pclose( xferPropList );
        assert( status > -1 );
        //logger << "xferPropList closed." << std::endl;

        /* Close the dataset. */
        status = H5Dclose( indatasetId );
        assert( status > -1 );
        //logger << "indatasetId closed." << std::endl;

        /* Close the file. */
        status = H5Fclose( infileId );
        assert( status > -1 );
        //logger << "infileId closed." << std::endl;


        //logger << "parallel_input_from_HDF5() done." << std::endl;
  

        /*
	
          The following code is just for debugging purposes!
	
        */

#ifdef DEBUG_LOG_LEVEL_2

        logger << "Now map it to the leaf grid elements..." << std::endl;
        for (LeafIterator it = gv.template begin<0,Dune::All_Partition> ()
               ; it != gv.template end<0,Dune::All_Partition> ()
               ; ++it) 
          {
            //int local1d_index = indexset.index( *it );
            //int map_index = mapper.map( *it );

            Dune::FieldVector<CTYPE, dim> x = it->geometry().center();
	  
            std::vector<UINT> global_index(dim,0);
            std::vector<UINT> local_index( dim, 0 );

            for(UINT i=0; i<dim; i++ )
              {
                global_index[i] = (UINT) floor( x[i] / inputdata.domain_data.virtual_gridsizes[i] );
                local_index[i] = global_index[i] - local_offset[ i ];
              }


#ifdef DIMENSION3
            UINT l = indexconversion_3d_to_1d( local_index[2], local_index[1], local_index[0]
                                               , local_count[2], local_count[1], local_count[0] );
#else
            UINT l = indexconversion_2d_to_1d( local_index[1], local_index[0]
                                               , local_count[1], local_count[0] );
#endif

            logger << "global_index = (" 
                   << global_index[0] << "," << global_index[1] 
#ifdef DIMENSION3
                   << "," << global_index[2] 
#endif
                   << ")   ";
            logger << "value = " 
                   <<  local_Yfield_vector[ l ] 
                   << std::endl;
	  
          }
        logger << std::endl;

#endif// DEBUG_LOG_LEVEL_2


        return;
      }
#endif // ifdef PARALLEL

      /** function to write a vector (sequential) to a hdf5 file 
       * 
       * \tparam dimensions the dimensions of the written data vector
       * \tparam data this data will be written to the file
       * \tparam data_name is the name/path where the data are stored in the HDF5 file
       * \tparam data_filename is the filename of the data file
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */

      static void write_sequential_to_HDF5_without_DUNE(
                                                 const Vector<UINT> dimensions
                                                 , const Vector<REAL> &data
                                                 , const std::string& data_name
                                                 , const std::string& data_filename
                                                 )
      {
        // ge the dimensionality of the data
        UINT dim=dimensions.size();
  
        //number of elements the data vector should have
        UINT n=1;
        for(UINT i=0;i<dim;i++)
          n*=dimensions[i];
  
        //the data make no sense!
        assert(n==data.size());

        /* Create a new file using default properties. */
        hid_t file_id = H5Fcreate( 
                                  data_filename.c_str()    // name of the file to be created
                                  , H5F_ACC_TRUNC          // if you are trying to create a file that exists already, the existing file will be truncated, i.e., all data stored on the original file will be erased
                                  , H5P_DEFAULT            // default file creation property list
                                  , H5P_DEFAULT            // default file access property list
                                   );

        assert( file_id > -1 );

        // create the memspace
        hsize_t mdims[1];
        mdims[0] = n;
        hid_t memspace_id = H5Screate_simple(           // H5Screate_simple creates a new simple dataspace and opens it for access, returning a dataset identifier.
                                             1          // rank=1 is the number of dimensions used in the dataspace
                                             , mdims    // mdims is a one-dimensional array of size rank specifying the size of each dimension of the dataset.
                                             , NULL     // maxdims is an array of the same size specifying the upper limit on the size of each dimension. maxdims may be NULL, in which case the maxdims == mdims. 
                                                        ); // must be released with H5Sclose or resource leaks will occur

        assert( memspace_id > -1 );

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

        //needed variable for HDF5 operations
        herr_t status;


        /*
         * chunks
         * !!! not optimized. Might be better if larger! -> use the inputdata 
         */
        hsize_t chunk_dims[ dim ];
        for(UINT i=0;i<dim;i++)
          chunk_dims[i]=1;  //maybe to small for fast saving!

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
                                     , HDF5_DATA_TYPE   // Datatype identifier, here: IEEE floating point 32-bit little-endian
                                     , dataspace_id      // Dataspace identifier 
                                     , plist_id          // Dataset creation property list identifier
                                      );
        assert( dataset_id > -1 );

        /* Write the dataset. */
        status = H5Dwrite(                         // Writes raw data from a buffer to a dataset.
                          dataset_id               // dataset identifier
                          , H5T_NATIVE_DOUBLE      // memory datatype id
                          , memspace_id            // specifies the memory dataspace and the selection within it
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

        status = H5Sclose( memspace_id );
        assert( status > -1 );

        status = H5Dclose( dataset_id );
        assert( status > -1 );

        status = H5Pclose( plist_id );
        assert( status > -1 );

        status = H5Fclose( file_id );
        assert( status > -1 );
      };



      /** function to retrieve a vector (sequential) from the hdf5 file
       *
       * \tparam local_data is a 'return value' which gets the data that belongs to the current processor (current hyperslab)
       * \tparam data_name is the name/path where the data are stored in the HDF5 file
       * \tparam data_filename is the filename of the data file
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */

      static void read_sequential_from_HDF5_without_DUNE(
                                                  Vector<REAL>& local_data
                                                  , const std::string& data_name
                                                  , const std::string& data_filename
                                                  )
      {
 
  
        // open the file for reading
        hid_t file_id=  H5Fopen (data_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);  
        assert( file_id > -1 );
  
        // open the dataset
        hid_t dataset_id = H5Dopen(file_id, data_name.c_str());
        assert( dataset_id > -1 );
    
        // get the dataspace
        hid_t dataspace_id=H5Dget_space (dataset_id); 
        assert(dataspace_id>-1);

        // some needed variables 
        herr_t status;
        hsize_t dim,*dims;
  
        // get the dimension (2-d or 3-d)
        dim = H5Sget_simple_extent_ndims( dataspace_id );
   
        // get the size of the problem
        dims=(hsize_t*)malloc(dim * sizeof (hsize_t));
        status = H5Sget_simple_extent_dims( dataspace_id , dims, 0);
        assert( status > -1 );
 
        UINT local_size=1;
        hsize_t offset[dim];
        for( UINT i=0; i<dim; i++ )
          {
            local_size*=dims[i];
            offset[i]=0;
          }
  
        // create the memory space
        hid_t memspace_id = H5Screate_simple (dim, dims, NULL);
    
        //select the hyperslab-
        status = H5Sselect_hyperslab (memspace_id, H5S_SELECT_SET, offset, NULL,
                                      dims, NULL);
        assert(status>-1);

        //resize the return data
        local_data.resize( local_size );
  
        /* set up the collective transfer properties list */
        hid_t xferPropList = H5Pcreate (H5P_DATASET_XFER);
        assert( xferPropList > -1 );

        //resize the return vector!
        local_data.resize(local_size);
  
        // finally the reading from the file
        status = H5Dread(                     dataset_id
                                              , H5T_NATIVE_DOUBLE //image.DataType
                                              , memspace_id
                                              , dataspace_id // H5S_ALL // 
                                              , xferPropList //H5P_DEFAULT // 
                                              , &( local_data[0] )
                                              );
        assert(status>-1);
  
        // close the identifyers
        H5Dclose (dataset_id);
        H5Sclose (dataspace_id);
        H5Sclose (memspace_id);
        free(dims);
        status = H5Fclose( file_id );
        assert( status > -1 );
      };

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

      static void write_parallel_to_HDF5_without_DUNE( 
                                               const Vector<UINT>& global_dim
                                               , const Vector<REAL> &data
                                               , const Vector<UINT>& local_count
                                               , const Vector<UINT>& local_offset
                                               , MPI_Comm communicator
                                               //, Dune::MPIHelper& helper
                                               , const std::string& data_name
                                               , const std::string& data_filename
                                               , const bool writemetafile = false
                                               , const std::string& meta_title="default"
#ifdef DIMENSION3
                                               , const std::string& meta_filename="DATA/3D/default.xmf"
#else
                                               , const std::string& meta_filename="DATA/2D/default.xmf"				   
#endif
                                                )
      {
  
        //Info varibale need for the HDF5
        MPI_Info mpiInfo = MPI_INFO_NULL;
        herr_t status;             /* Generic return value */
  
        //get the dimension of the problem -> no cross checking ig the given data make any sense!
        UINT dim = local_count.size();

        // Set up file access property list with parallel I/O access
        hid_t plist_id= H5Pcreate( H5P_FILE_ACCESS );

  
        // Set up file access property list with parallel I/O access
        status=H5Pset_fapl_mpio(plist_id, communicator, mpiInfo);  //collcetive MPI!!! needs to be called on each processor of the communicator
        assert(plist_id>-1);
   
        // Create a new file using default properties.
        hid_t file_id= H5Fcreate( 
                                 data_filename.c_str()    // name of the file to be created
                                 , H5F_ACC_TRUNC          // if you are trying to create a file that exists already, the existing file will be truncated, i.e., all data stored on the original file will be erased
                                 , H5P_DEFAULT            // default file creation property list
                                 , plist_id
                                  );
        assert( file_id > -1 );
        H5Pclose(plist_id);

        // set the global size of the grid into a vector of type hsize_t (needed for HDF5 routines)
        hsize_t global_dim_HDF5[ dim];
        for(UINT i=0; i<dim;i++)
          global_dim_HDF5[dim -i-1]=global_dim[i];
    
        // set the count and offset in the different dimensions (determine the size of the hyperslab) (in hsize_t format, needed for HDF5 routines)
        hsize_t count[ dim ],offset[ dim ];
        for(UINT i=0;i<dim;i++){
          count[dim-i-1]=local_count[i];
          offset[dim-i-1]=local_offset[i];   
        }
  
        //define the total size of the local data
        hsize_t nAllLocalCells = count[0] * count[1];
#ifdef DIMENSION3
        nAllLocalCells *= count[2];
#endif
  
        //  Create the dataspace for the dataset.
        hid_t filespace = H5Screate_simple(dim, global_dim_HDF5, NULL);
        assert(filespace>-1);
  
        // Create the dataset with default properties and close filespace.
        hid_t dset_id = H5Dcreate(file_id, data_name.c_str(), HDF5_DATA_TYPE, filespace,
                                  H5P_DEFAULT);
        H5Sclose(filespace);
        assert(dset_id>-1);
 
        //get the memoryspace (but only if something needs to be written on this processor!)
        hid_t memspace_id;
        if(nAllLocalCells!=0){  // -> otherwise HDF5 warning, because of writing nothing!
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
    
        // fineally write the data to the disk
        // even if nothing should be written H5Dwrite needs to be called!!
        if(nAllLocalCells!=0){ // -> otherwise HDF5 warning, because of writing nothing!
          status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace,
                            plist_id, &( data[0] )  );
          assert(status>-1);
        }else{ // IMPORTANT. otherwise the H5Dwrite() blocks!!!
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
        if(nAllLocalCells!=0){   
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
        MPI_Barrier(communicator);


        // write meta file if needed
        if(writemetafile)
          /* Write meta file for display in paraview */
          writeMetaFileForHDF5( meta_filename, data_filename, data_name, &(global_dim_HDF5[0]), meta_title  );

      };



      /** function to retrieve a vector (sequential) from the hdf5 file
       *
       * \tparam local_data is a 'return value' which gets the data that belongs to the current processor (current hyperslab)
       * \tparam local_count give the size of the local data
       * \tparam local_offset the offset of the data (in each direction)
       *  \tparam helper the DUNE MPIHelper for the MPI communication
       * \tparam data_name is the name/path where the data are stored in the HDF5 file
       * \tparam data_filename is the filename of the data file
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */

      static void read_parallel_from_HDF5_without_DUNE(
                                                Vector<REAL>& local_data
                                                , const Vector<UINT>& local_count
                                                , const Vector<UINT>& local_offset
                                                , MPI_Comm communicator
                                                //, Dune::MPIHelper& helper
                                                , const std::string& data_name
                                                , const std::string& data_filename
                                                )
      {
        /* setup file access template with parallel IO access. */
        hid_t access_pList = H5Pcreate( H5P_FILE_ACCESS );
        assert( access_pList > -1 );

        herr_t status;             /* Generic return value */
        MPI_Info mpiInfo = MPI_INFO_NULL;
        status = H5Pset_fapl_mpio( access_pList, communicator, mpiInfo ); // Take care, this requires the MPI version of hdf5!
        assert( status > -1 );
  
        // open the file for reading
        hid_t file_id=  H5Fopen (data_filename.c_str(), H5F_ACC_RDONLY, access_pList);  
        assert( file_id > -1 );
  
        /* Release file-access template */
        status = H5Pclose( access_pList );
        assert( status > -1 );
  
        // open the dataset
        hid_t dataset_id = H5Dopen(file_id, data_name.c_str());
        assert( dataset_id > -1 );
    
        // get the dataspace
        hid_t dataspace_id=H5Dget_space (dataset_id); 
        assert(dataspace_id>-1);

        // some needed variables 
        hsize_t dim,*dims;
  
        // get the dimension (2-d or 3-d)
        dim = H5Sget_simple_extent_ndims( dataspace_id );
   
        // get the size of the data structure
        dims=(hsize_t*)malloc(dim * sizeof (hsize_t));
        status = H5Sget_simple_extent_dims( dataspace_id , dims, 0);
        assert( status > -1 );
  
        //set the local, offset, and count as hsize_t, which is needed by the HDF5 routines 
        hsize_t local_size=1;
        hsize_t offset[dim],count[dim];
        for( UINT i=0; i<dim; i++ )
          {
            local_size*=local_count[i];
            offset[dim-i-1]=local_offset[i];
            count[dim-i-1]=local_count[i];
          }
  
        // create the memory space, if something needes to be read on this processor
        hid_t memspace_id;
        if(local_size!=0)
          memspace_id = H5Screate_simple (1, &local_size, NULL);
    
        //select the hyperslab
        status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, NULL,
                                      count, NULL);
        assert(status>-1);

        //resize the return data
        local_data.resize( local_size );
  
        // set up the collective transfer properties list 
        hid_t xferPropList = H5Pcreate (H5P_DATASET_XFER);
        assert( xferPropList > -1 );
        //logger << "xferPropList created." << std::endl;

        // finally the reading from the file, only if something needes to be written
        if(local_size!=0){
          status = H5Dread(                     dataset_id
                                                , H5T_NATIVE_DOUBLE //image.DataType
                                                , memspace_id
                                                , dataspace_id // H5S_ALL // 
                                                , xferPropList //H5P_DEFAULT // 
                                                , &( local_data[0] )
                                                );
          assert(status>-1);
        }

        // close the identifyers
        H5Dclose (dataset_id);
        H5Sclose (dataspace_id);
        if(local_size!=0) //this identifier only exists if somethings needs to be read
          H5Sclose (memspace_id);
        free(dims);
        status = H5Fclose( file_id );
        assert( status > -1 );
      };

      /** function to write a vector (parallel) to an EXISTING hdf5 file
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

      static void write_parallel_to_existing_HDF5_without_DUNE( 
                                                        const Vector<REAL> &data
                                                        , const Vector<UINT>& local_count
                                                        , const Vector<UINT>& local_offset
                                                        , MPI_Comm communicator
                                                        //, Dune::MPIHelper& helper
                                                        , const std::string& data_name
                                                        , const std::string& data_filename
                                                         )
      {
  
        // setup file access template with parallel IO access. 
        hid_t access_pList = H5Pcreate( H5P_FILE_ACCESS );
        assert( access_pList > -1 );

        herr_t status;             // Generic return value 
        MPI_Info mpiInfo = MPI_INFO_NULL;
        status = H5Pset_fapl_mpio( access_pList, communicator, mpiInfo ); // Take care, this requires the MPI version of hdf5!
        assert( status > -1 );
  
        // open the file for reading
        hid_t file_id=  H5Fopen (data_filename.c_str(), H5F_ACC_RDWR, access_pList);  
        assert( file_id > -1 );
        H5Pclose(access_pList);
  
        // open the dataset
        hid_t dataset_id = H5Dopen(file_id, data_name.c_str());
        assert( dataset_id > -1 );
    
        // get the dataspace
        hid_t dataspace_id=H5Dget_space (dataset_id); 
        assert(dataspace_id>-1);
 

        // some needed variables 
        hsize_t dim,*dims;
 
        // get the dimension (2-d or 3-d)
        dim = H5Sget_simple_extent_ndims( dataspace_id );
   
        // get the size of the data structure
        dims=(hsize_t*)malloc(dim * sizeof (hsize_t));
        status = H5Sget_simple_extent_dims( dataspace_id , dims, 0);
        assert( status > -1 );
  
  
        // set the count and offset in the different dimensions (determine the size of the hyperslab) (in hsize_t format, needed for HDF5 routines)
        hsize_t count[ dim ],offset[ dim ];
        for(UINT i=0;i<dim;i++){
          count[dim-i-1]=local_count[i];
          offset[dim-i-1]=local_offset[i];   
        }
  
  
        //define the total size of the local data
        hsize_t nAllLocalCells = count[0] * count[1];
#ifdef DIMENSION3
        nAllLocalCells *= count[2];
#endif
  
        //get the memoryspace (but only if something needs to be written on this processor!)
        hid_t memspace_id;
        if(nAllLocalCells!=0){  // -> otherwise HDF5 warning, because of writing nothing!
          memspace_id = H5Screate_simple(dim, count, NULL);
          assert(memspace_id>-1);
        }
  
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
  
        access_pList = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(access_pList, H5FD_MPIO_COLLECTIVE); // all MPI processors of this communicator need to call this (if they write or not)
  
        // fineally write the data to the disk
        // even if nothing should be written H5Dwrite needs to be called!!
        if(nAllLocalCells!=0){ // -> otherwise HDF5 warning, because of writing nothing!
          status = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,
                             access_pList, &( data[0] )  );
          assert(status>-1);
        }else{ // IMPORTANT. otherwise the H5Dwrite() blocks!!!
          status = H5Dwrite( dataset_id, 
                             H5T_NATIVE_DOUBLE, 
                             H5S_ALL, 
                             dataspace_id,
                             access_pList, 
                             &( data[0] )  );
          assert(status>-1);
        }
  
        // close the identifyers
        H5Dclose (dataset_id);
        H5Sclose (dataspace_id);
        // Release file-access template 
        status = H5Pclose( access_pList );
        assert( status > -1 );
        if(nAllLocalCells!=0) //this identifier only exists if somethings needs to be read
          H5Sclose (memspace_id);
        free(dims);
        status = H5Fclose( file_id );
        assert( status > -1 );
  
        return;
      };





      
      template<typename NUM,typename GFS>
      static void write_vector_to_HDF5( std::vector<NUM>& vData,
                                        const GFS &gfs,
                                        const std::string & filename, 
                                        const std::string & datasetname
                                        ){

        const typename GFS::Traits::GridViewType& gv = gfs.gridView();

        UINT local_dofs = gfs.size();
        UINT total_dofs = gfs.size();
        total_dofs = gv.comm().sum( total_dofs );

        UINT global_dofs = gfs.size();
        global_dofs = gv.comm().max( global_dofs );

        logger << "local_dofs = " << local_dofs << std::endl;
        logger << "total_dofs = " << total_dofs << std::endl;
        
        int mpi_size = gv.comm().size();
        int mpi_rank = gv.comm().rank();
        std::vector<int> v1(mpi_size,0);
        v1[ mpi_rank ] = local_dofs;

        std::vector<int> v2(mpi_size,0);
        
        MPI_Allreduce( &(v1[0]), &(v2[0]), mpi_size, MPI_INT, MPI_SUM, gv.comm() );
        int my_offset = 0;
        for(int i=0;i<mpi_rank;i++){
          my_offset += v2[i];
        }
        logger << "my_offset = " << my_offset << std::endl;

        const int h5rank = 2;
        hid_t plist_id;          /* property list identifier */
        hid_t file_id;
        hid_t dset_id;           /* file and dataset identifiers */
        hid_t filespace_id;
        hid_t memspace_id;       /* file and memory dataspace identifiers */
        hsize_t dimsf[h5rank];   /* dataset dimensions */
        hsize_t	count[h5rank];	 /* hyperslab selection parameters */
        hsize_t	offset[h5rank];
        double *data;            /* pointer to data buffer to write */
        herr_t status;

  
        //Info varibale need for the HDF5
        MPI_Info mpiInfo = MPI_INFO_NULL;
  
        //some needed variables and typdefs
        //const UINT dim = GV::dimension;  // dimensionality of the problem (2-d or 3-d)

        /* 
         * Set up file access property list with parallel I/O access
         */
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id,gv.comm(),mpiInfo);

        /*
         * Create a new file collectively and release property list identifier.
         */
        file_id = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id );
        H5Pclose(plist_id);

        /*
         * Create the dataspace for the dataset.
         */
        dimsf[0] = total_dofs;
        dimsf[1] = 1;
        filespace_id = H5Screate_simple(h5rank,dimsf,NULL);
        /*
         * Create the dataset with default properties and close filespace.
         */
        dset_id = H5Dcreate(file_id,
                            datasetname.c_str(),
                            H5T_NATIVE_DOUBLE,
                            filespace_id,
                            H5P_DEFAULT);
        H5Sclose(filespace_id);
        /* 
         * Each process defines dataset in memory and writes it to the hyperslab
         * in the file.
         */
        count[0] = local_dofs;
        count[1] = 1;
        offset[0] = my_offset;
        offset[1] = 0;
        memspace_id = H5Screate_simple(h5rank,count,NULL);
        /*
         * Select hyperslab in the file.
         */
        filespace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace_id,H5S_SELECT_SET,offset,NULL,count,NULL);
        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  
        status = H5Dwrite( dset_id,
                           H5T_NATIVE_DOUBLE,
                           memspace_id,
                           filespace_id,
                           plist_id,
                           &(vData[0]) );
        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace_id);
        H5Sclose(memspace_id);
        H5Pclose(plist_id);
        H5Fclose(file_id);
      };





      // Save BackendVector to hdf5 (This implementation works not only for Q1 finite elements!)
      // read_BackendVector_from_HDF5()
      template<typename GFS>
      static void write_BackendVector_to_HDF5( 
                                              const GFS& gfs
                                              , const CInputData& inputdata
                                              , const std::string & filename 
                                              , const std::string & data_name
                                              , const typename Dune::PDELab::BackendVectorSelector<GFS,REAL>::Type& backend_vector
                                               ){

        std::vector<REAL> standard_vector;
        backend_vector.std_copy_to( standard_vector );

        write_vector_to_HDF5( standard_vector, 
                              gfs,
                              filename, 
                              data_name 
                              );

      };










      template<typename NUM,typename GFS>
      static void read_vector_from_HDF5( std::vector<NUM>& vData,
                                         const GFS& gfs,
                                         const std::string & filename, 
                                         const std::string & datasetname
                                         ){

        const typename GFS::Traits::GridViewType& gv = gfs.gridView();

        UINT local_dofs = gfs.size();
        UINT total_dofs = gfs.size();
        total_dofs = gv.comm().sum( total_dofs );

        UINT global_dofs = gfs.size();
        global_dofs = gv.comm().max( global_dofs );

        logger << "local_dofs = " << local_dofs << std::endl;
        logger << "total_dofs = " << total_dofs << std::endl;
        
        int mpi_size = gv.comm().size();
        int mpi_rank = gv.comm().rank();
        std::vector<int> v1(mpi_size,0);
        v1[ mpi_rank ] = local_dofs;

        std::vector<int> v2(mpi_size,0);
        
        MPI_Allreduce( &(v1[0]), &(v2[0]), mpi_size, MPI_INT, MPI_SUM, gv.comm() );
        int my_offset = 0;
        for(int i=0;i<mpi_rank;i++){
          my_offset += v2[i];
        }
        logger << "my_offset = " << my_offset << std::endl;

        
        const int h5rank = 2;
        hid_t plist_id;          /* property list identifier */
        hid_t file_id;
        hid_t dset_id;           /* file and dataset identifiers */
        hid_t filespace_id;
        hid_t memspace_id;       /* file and memory dataspace identifiers */
        hsize_t dimsf[h5rank];   /* dataset dimensions */
        hsize_t	count[h5rank];	 /* hyperslab selection parameters */
        hsize_t	offset[h5rank];
        double *data;            /* pointer to data buffer to write */
        herr_t status;

  
        //Info varibale need for the HDF5
        MPI_Info mpiInfo = MPI_INFO_NULL;

        /* 
         * Set up file access property list with parallel I/O access
         */
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id,gv.comm(),mpiInfo);

        /* open the file collectively */
        file_id = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, plist_id );

        H5Pclose(plist_id);

        dset_id = H5Dopen( file_id, datasetname.c_str() );
    
        hid_t dspace_id = H5Dget_space( dset_id );
    
        hsize_t dims[ h5rank ];
        status = H5Sget_simple_extent_dims( dspace_id, dims, 0 );


        /* 
         * Each process defines dataset in memory and reads it from the hyperslab
         */
        count[0] = local_dofs;
        count[1] = 1;
        offset[0] = my_offset;
        offset[1] = 0;

        memspace_id = H5Screate_simple( h5rank, count, NULL );

        status = H5Sselect_hyperslab( dspace_id
                                      , H5S_SELECT_SET
                                      , offset
                                      , NULL
                                      , count
                                      , NULL
                                      );
        
        vData.resize( local_dofs );
        status = H5Dread( dset_id
                          , H5T_NATIVE_DOUBLE
                          , memspace_id
                          , dspace_id
                          , H5P_DEFAULT
                          , &( vData[0] )
                          );
        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(dspace_id);
        H5Sclose(memspace_id);
        H5Fclose(file_id);

        //vData.resize( gfs.size() ); // cut vector back to original size

      };





      // Read BackendVector from hdf5 (This implementation works only for Q1 finite elements!)
      // <---> write_BackendVector_to_HDF5()
      template<typename GFS>
      static void read_BackendVector_from_HDF5( const GFS& gfs
                                                , const CInputData& inputdata
                                                , const std::string & filename
                                                , const std::string & data_name
                                                , typename Dune::PDELab::BackendVectorSelector<GFS,REAL>::Type& backend_vector
                                                ){

        std::vector<REAL> read_local_data;
        read_vector_from_HDF5( read_local_data, gfs, filename, data_name );
        backend_vector.std_copy_from( read_local_data );

        logger << "Read backend-vector from " << filename.c_str() << std::endl;
        logger << "standard vector size = " << read_local_data.size() << std::endl;

      };


#endif // USE_HDF5

    }; // class HDF5Tools


  } // GeoInversion
} // Dune








#endif	/* DUNE_GEOINVERSION_HDF5_TOOLS_HH */

