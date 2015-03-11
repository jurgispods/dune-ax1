/*
 * output.hh
 *
 *  Created on: May 11, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_OUTPUT_HH
#define DUNE_AX1_OUTPUT_HH

#include <iostream>
#include <iomanip>
#include <string>
#include <valarray>
#include <vector>
#include <fstream>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/ios_state.hh>

#include <dune/pdelab/backend/istlvectorbackend.hh>

#include <dune/ax1/common/constants.hh>


class Output
{
  public:

    static const int prec = 12;
    static constexpr double eps = 1e-8;

    static void gnuplotInitialize(const std::string& filename,
                                  const std::string& prefix = "",
                                  const std::string& header = "")
    {
      std::ofstream file;
      file.open((prefix + filename).c_str() );
      if (header != "") file << header << std::endl;
      file.close();
    }

    // add xy data to file

    template<typename T>
    static void gnuplotAppend( const std::string& filename,
                               const T& x,
                               const T& y,
                               const std::string& prefix = "",
                               const std::string& infoString = "")
    {
      std::ofstream file;
      file.open((prefix + filename).c_str(), std::ofstream::app);

      if(infoString != "")
      {
        file << std::endl;
        file << std::endl;
        file << "# " << infoString << std::endl;
      }

      file << std::setprecision(prec) << std::scientific << x << " " << y << std::endl;
      file.close();
    }
    
    // add xyz data to file
    
    template<typename T>
    static void gnuplotDoubleAppend( const std::string& filename,
                                     const T& x,
                                     const T& y,
                                     const T& z,
                                     const std::string& prefix = "",
                                     const std::string& infoString = "" )
    {
      std::ofstream file;
      file.open((prefix + filename).c_str(), std::ofstream::app);

      if(infoString != "")
      {
        file << std::endl;
        file << std::endl;
        file << "# " << infoString << std::endl;
      }

      file << std::setprecision(prec) << std::scientific << x << " " << y << " " << z << std::endl;
      file.close();
    }
    
    // add a row of diagnostic values to file

    template<typename T>
    static void gnuplotMultiAppend( const std::string& filename,
                                    const T& x,
    																const std::vector<T>& y,
    																const std::string& prefix = "",
    																const std::string& infoString = "")
		{
			std::ofstream file;
			file.open((prefix + filename).c_str(), std::ofstream::app);

			if(infoString != "")
			{
			  file << std::endl;
        file << std::endl;
			  file << "# " << infoString << std::endl;
			}

			file << std::setprecision(prec) << std::scientific << x << " ";
			for(int i=0; i<y.size(); ++i)
			{
			  file << std::scientific << y[i] << " ";
			}
			file << std::endl;
			file.close();
		}

    // add mutiple diagnostic arrays to a gnuplot file
    template<typename T>
    static void gnuplotAppendMultiArray( const std::string& filename,
                                    const std::vector<std::vector<T> >& diagValues,
                                    const std::string& infoString = "",
                                    const std::string& prefix = "" )
    {
      std::ofstream file;
      file.open((prefix + filename).c_str(), std::ofstream::app);

      file << std::endl;
      file << std::endl;

      if(infoString != "")
      {
        file << "# " << infoString << std::endl;
      }

      for(int j=0; j<diagValues[0].size(); ++j)
      {
        // Print iteration number
        file << std::fixed << std::setw(3) << j << " ";
        for(int i=0; i<diagValues.size(); ++i)
        {
          file << std::setprecision(prec) << std::scientific << diagValues[i][j] << " ";
        }
        file << std::endl;
      }
      file.close();
    }

    // write two valarrays to file, deleting all previous entries

    template<typename T>
    static void gnuplotArray(const std::string& filename, const std::valarray<T>& x, const std::valarray<T>& y,
        const std::string& prefix = "")
    {
      std::ofstream file;
      file.open((prefix + filename).c_str());
      if (! file.is_open())
        DUNE_THROW(Dune::IOError, "Could not write to file " << filename);

      for(int i=0; i<x.size(); ++i)
      {
        file << std::setprecision(prec) << std::scientific << x[i] << " " << y[i] << std::endl;
      }

      file.close();
    }
    
    // write two or more valarrays to file, deleting all previous entries
    
    template<typename T>
    static void gnuplotMultiArray(const std::string& filename,
                                  const std::valarray<T>& x,
                                  const std::vector<std::valarray<T> >& y,
                                  const std::string& prefix = "")
    {
      std::ofstream file;
      file.open((prefix + filename).c_str());
      if (! file.is_open())
        DUNE_THROW(Dune::IOError, "Could not write to file " << filename);

      for(int i=0; i<x.size(); ++i)
      {
        file << std::setprecision(prec) << std::scientific << x[i];
        for (int j=0; j<y.size(); j++)
        {
          file << " " << y[j][i];
        }
        file << std::endl;
      }

      file.close();
    }
    
    // add valarrays to file, separated by two blank lines from previous entry (accessible with
    // "index" in gnuplot), optionally with a info line starting with "#" (to be skipped in gnuplot)
    // (for DG: DG elements are plotted "discontinuous")
    template<typename T>
    static void gnuplotAppendArray(const std::string& filename,
                                   	 const std::valarray<T>& x,
                                   	 const std::valarray<T>& y,
                                   	 const std::string& infoString = "",
                                   	 const std::string& prefix = "")
    {
      std::ofstream file;
      file.open((prefix + filename).c_str(), std::ofstream::app);
      if (! file.is_open())
        DUNE_THROW(Dune::IOError, "Could not write to file " << filename);

      file << "" << std::endl;
      file << "" << std::endl;
      if ( infoString != "") file << "# " << infoString << std::endl;

      for(int i=0; i<x.size(); ++i)
      {
        file << std::setprecision(prec) << std::scientific << x[i] << " " << y[i] << std::endl;
        //if (i>0 and std::abs(x[i+1]-x[i]) < eps)
        //{
      	//  file << "" << std::endl;
        //}
      }

      file.close();
    }

    // add valarrays to file, seperated by two blank lines from previous entry (accesible with
		// "index" in gnuplot), optionally with a info line starting with "#" (to be skipped in gnuplot)
		// (for DG: DG elements are plotted "discontinuous")
		template<typename T>
		static void gnuplotAppendDoubleArray(const std::string& filename,
																		     const std::valarray<T>& x,
																		     const std::valarray<T>& y,
																		     const std::valarray<T>& z,
																		     const std::string& infoString = "",
																		     const std::string& prefix = "")
		{
			std::ofstream file;
			file.open((prefix + filename).c_str(), std::ofstream::app);
			if (! file.is_open())
				DUNE_THROW(Dune::IOError, "Could not write to file " << filename);

			file << "" << std::endl;
			file << "" << std::endl;
			if ( infoString != "") file << "# " << infoString << std::endl;

			for(int i=0; i<x.size(); ++i)
			{
				file << std::setprecision(prec) << std::scientific << x[i] << " " << y[i] << " " << z[i] << std::endl;
				//if (i>0 and std::abs(x[i+1]-x[i]) < eps)
				//{
				//	file << "" << std::endl;
				//}
			}

			file.close();
		}

    template<typename T>
    static void printVector(const std::vector<T>& v, const std::string prefix = "")
    {
      debug_info << prefix << "[ ";
      for(int i=0; i<v.size(); ++i)
      {
        debug_info << v[i] << " ";
      }
      debug_info << "]" << std::endl;
    }

    template<typename Iterator>
    static void printVectorInterval(Iterator start, Iterator end,
        const std::string prefix = "")
    {
      debug_info << prefix << "[ ";
      for(Iterator it = start; it != end; ++it)
      {
        debug_info << *it << " ";
      }
      debug_info << "]" << std::endl;
    }

    template<typename T>
    static void printValarray(const std::valarray<T>& v)
    {
      debug_info << "[ ";
      for(int i=0; i<v.size(); ++i)
      {
        debug_info << v[i] << " ";
      }
      debug_info << "]" << std::endl;
    }

    template<typename U>
    static void printRawCoefficientVector(const U& u, const char* name, const int precision = prec)
    {
      Dune::ios_base_all_saver restorer(std::cout); // store old ios flags

      debug_info << "------------ " << name  << " -------" << std::endl;
      debug_info << std::scientific;
      std::size_t i = 0;
      for (typename U::const_iterator it = u.begin(); it != u.end(); ++it)
      {
        debug_info << "u[x_" << i << "]:  |";
        debug_info << std::setprecision(precision) << std::setfill(' ') << *it  << "  |";
        debug_info << std::endl;
        i++;
      }
      debug_info << "---------------------------" << std::endl;
    }

    //! \brief Print multiple coefficient vectors side by side for direct comparison
    template<typename U>
    static void printRawCoefficientVector(const std::vector<U>& u, const std::vector<std::string>& name,
        const int precision = prec)
    {
      Dune::ios_base_all_saver restorer(std::cout); // store old ios flags

      std::vector<typename U::const_iterator> it;
      // Initialize all iterators
      for(int j=0; j<u.size(); j++)
      {
        it.push_back(u[j].begin());
        debug_info << "----- " << name[j];
      }
      debug_info << " -----" << std::endl;
      debug_info << std::scientific;

      std::size_t i = 0;

      // Now loop over all coefficient vectors and print components side by side
      for (; it[0] != u[0].end(); ++it[0])
      {
        debug_info << "u[x_" << i << "]:  |";
        for(int j=0; j<u.size(); j++)
        {
          debug_info << std::setprecision(precision) << std::setfill(' ') << *(it[j])  << "  |";
          if(j > 0)
          {
            ++it[j];
          }
        }
        debug_info << std::endl;
        i++;
      }
      debug_info << "---------------------------" << std::endl;
    }

    template<typename U>
    static void printPermutedCoefficientVector(const U& u, const char* name, const std::vector<int>& permutation)
    {
      assert(u.flatsize() == permutation.size());

      Dune::ios_base_all_saver restorer(std::cout); // store old ios flags

      debug_info << "------------ " << name  << " -------" << std::endl;
      debug_info << std::scientific;
      std::size_t i = 0;
      std::vector<double> raw(u.flatsize());
      for (typename U::const_iterator it = u.begin(); it != u.end(); ++it)
      {
        raw[i] = *it;
//        std::cout << "u[x_" << i << "]:  |";
//        std::cout << std::setprecision(prec) << std::setfill(' ') << raw[i] << " (" << *it << ")  |";
//        std::cout << std::endl;
        i++;
      }
      for(i=0; i<raw.size(); i++)
      {
        debug_info << "u[x_" << i << "]:  |";
        debug_info << std::setprecision(prec) << std::setfill(' ') << raw[permutation[i]]  << "  |";
        debug_info << std::endl;
      }
      debug_info << "---------------------------" << std::endl;
    }


    /**
     * Using this method for numEquations > 1 only makes sense when lexicographic ordering is
     * used to order the coefficient vector u. Otherwise use numEquations = 1 for printing the
     * plain vector 'as is'
     *
     * @param u
     * @param numEquations
     */
    template<typename U>
    static void printCoefficientVector(const U& u, int numEquations)
    {
      //typedef typename U::Backend Backend;

      Dune::ios_base_all_saver restorer(std::cout); // store old ios flags

      // Don't try to group coeff vector equation-wise when blocksize > 1
      if(u.flatsize() != u.N()) numEquations = 1;

      assert(u.flatsize() % numEquations == 0);

      std::size_t size = u.flatsize() / numEquations;

      // First copy everything to a temporary std::vector;
      // TODO choose a more performant implementation here
      std::vector<typename U::field_type> stdVec(u.flatsize());
      std::size_t i = 0;
      for (typename U::const_iterator it = u.begin(); it != u.end(); ++it)
      {
        stdVec[i] = *it;
        i++;
      }

      Dune::ios_base_all_saver ohrensessel(std::cout);

      //debug_info << std::fixed << std::setprecision(4);
      debug_info << std::scientific;

      for(i=0; i<size; ++i)
      {
        debug_info << "u[x_" << i << "]:  |";
        for(std::size_t j=0; j<numEquations; ++j)
        {
          //typename U::field_type value = Backend::access(u,j*size+i);
          typename U::field_type value = stdVec[j*size+i];
          //double value = u[j*size+i];
          debug_info << std::setprecision(prec) << std::setfill(' ') << value  << "  |";
        }
        debug_info << std::endl;
      }
    }

    template<typename U>
    static void printSingleCoefficientVector(const U& u, const char* name)
    {
      debug_info << "------------ " << name  << " -------" << std::endl;
      Output::printCoefficientVector(u, 1);
      debug_info << "---------------------------" << std::endl;
    }

    template<typename U>
    static void printMultipleComponentCoefficientVector(const U& u, int numEquations)
    {
      debug_info << "--------------- ";
      for(int j=0; (j<numEquations && j<ION_NAMES.size()); ++j)
      {
        debug_info << ION_NAMES[j] << " ----------- ";
      }
      if(ION_NAMES.size() < numEquations)
      {
        debug_info << "pot ------------- ";
      }
      debug_info << std::endl;
      Output::printCoefficientVector(u, numEquations);
      debug_info << "-----------------------------------------------------------------------------" << std::endl;
    }

    // Helper method for printSingle... / printMultiple...
    template<typename U>
    static void printCoefficientVectorDG(const U& u, int numEquations)
    {
      Dune::ios_base_all_saver restorer(std::cout); // store old ios flags

      assert(u.N() % numEquations == 0);
      // Assume block size of 1
      assert(u.N() == u.flatsize());

      std::size_t size = u.N() / numEquations;
      std::size_t i;

      Dune::ios_base_all_saver camembert(std::cout);
      //debug_info << std::fixed << std::setprecision(4);
      debug_info << std::scientific;

      for(i=0; i<size; i=i+2)
      {
        debug_info << "u[x_" << i << "]:  |";
        for(std::size_t j=0; j<numEquations; ++j)
        {
          typename U::field_type value = u[j*size+i];
          //double value = u[j*size+i];
          debug_info << std::setprecision(prec) << std::setfill(' ') << value  << "  |";
        }
        debug_info << std::endl;
      }
    }

    template<typename U>
    static void printSingleCoefficientVectorDG(const U& u, const char* name)
    {
      debug_info << "------------ " << name  << " -------" << std::endl;
      Output::printCoefficientVectorDG(u, 1);
      debug_info << "---------------------------" << std::endl;
    }

    template<typename U>
    static void printMultipleComponentCoefficientVectorDG(const U& u, int numEquations)
    {
      debug_info << "--------------- ";
      for(int j=0; (j<numEquations && j<ION_NAMES.size()); ++j)
      {
        debug_info << ION_NAMES[j] << " ----------- ";
      }
      if(ION_NAMES.size() < numEquations)
      {
        debug_info << "pot ------------- ";
      }
      debug_info << std::endl;
      Output::printCoefficientVectorDG(u, numEquations);
      debug_info << "-----------------------------------------------------------------------------" << std::endl;
    }

};

#endif /* DUNE_AX1_OUTPUT_HH */
