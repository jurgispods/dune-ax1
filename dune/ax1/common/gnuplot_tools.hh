/*
 * gnuplot_tools.hh
 *
 *  Created on: May 16, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_GNUPLOT_TOOLS_HH
#define DUNE_AX1_GNUPLOT_TOOLS_HH

#include <valarray>
#include <vector>
#include <algorithm>
#include <iomanip>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/ax1_outputfunctors.hh>
#include <dune/ax1/common/ax1_output2d.hh>
#include <dune/ax1/common/ax1_parallelhelper.hh>

namespace
{
  enum GnuplotFunctions { GNUPLOT_ADD_BLOCK = 0, GNUPLOT_ADD_LINE = 1,
                        GNUPLOT_ADD_MULTI_COLUMN_BLOCK = 2};

  // Default output strategy
  template<typename AcmeOutputTraits, typename GF>
  struct OutputFunction
  {
   enum { value = GNUPLOT_ADD_BLOCK };
  };

  // Specialization for membrane potential
  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::GF_MEMB_POT>
  {
   enum { value = GnuplotFunctions::GNUPLOT_ADD_LINE };
  };

  // Specialization for membrane potential
  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::GF_MEMB_POT_MD>
  {
   enum { value = GnuplotFunctions::GNUPLOT_ADD_LINE };
  };

  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_FLUX>
  {
   enum { value = GnuplotFunctions::GNUPLOT_ADD_LINE };
  };

  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_FLUX_DIFF>
  {
   enum { value = GnuplotFunctions::GNUPLOT_ADD_LINE };
  };

  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_FLUX_DRIFT>
  {
   enum { value = GnuplotFunctions::GNUPLOT_ADD_LINE };
  };

  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_FLUX_LEAK>
  {
   enum { value = GnuplotFunctions::GNUPLOT_ADD_LINE };
  };

  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_FLUX_VOLTAGE_GATED>
  {
   enum { value = GnuplotFunctions::GNUPLOT_ADD_LINE };
  };

  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_FLUX_MORI>
  {
   enum { value = GnuplotFunctions::GNUPLOT_ADD_LINE };
  };

  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_MEMB_CURRENT>
  {
   enum { value = GnuplotFunctions::GNUPLOT_ADD_LINE };
  };

  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::GF_SINGLE_PARAM_MEMB_FLUX>
  {
   enum { value = GnuplotFunctions::GNUPLOT_ADD_LINE };
  };

  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::BGF_PARAM_POT_NEUMANN>
  {
   enum { value = GnuplotFunctions::GNUPLOT_ADD_LINE };
  };

  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::GF_CHANNEL>
  {
   enum { value = GnuplotFunctions::GNUPLOT_ADD_LINE };
  };

  template<typename AcmeOutputTraits>
  struct OutputFunction<AcmeOutputTraits, typename AcmeOutputTraits::GF_MEMB_GROUPS>
  {
    enum { value = GnuplotFunctions::GNUPLOT_ADD_BLOCK };
  };


}

template<typename OutputTraits>
class GnuplotTools2D
{
  public:
    static const int prec = 12;
    static constexpr double eps = 1e-8;

    template<typename PHYSICS, typename GF>
    static void writeGF(const PHYSICS& physics, const GF& gf,
        const int SUBSAMPLING_POINTS, const std::string filename, const double time,
        const std::string& infoStr = "", const bool forceWrite = false)
    {
      Ax1OutputFunctors::IdentityFunctor f;
      GnuplotTools2D<OutputTraits>::writeGF(physics,gf,SUBSAMPLING_POINTS,filename,time,infoStr,forceWrite,f);
    }

    template<typename PHYSICS, typename GF, typename FUNC>
    static void writeGF(const PHYSICS& physics, const GF& gf,
        const int SUBSAMPLING_POINTS, const std::string filename, const double time,
        const std::string& infoStr, const bool forceWrite, FUNC& f)
    {
      if(physics.getParams().general.get("suppressAllOutput",false))
        return;

      if(forceWrite || physics.getParams().general.get("gnuplotOutput",true))
      {
        Dune::Timer timer;

        std::vector<typename PHYSICS::Element::Geometry::GlobalCoordinate> pos;
        std::vector<typename GF::Traits::RangeType> sol;

        debug_jochen << "  [GnuplotTools2D] Writing gridfunction to file " << filename << "..." << std::endl;

        typename Ax1Output2D<OutputTraits>::SolutionVectorInfo info;
        Ax1Output2D<OutputTraits>::getSolutionVector(physics, gf, SUBSAMPLING_POINTS, pos, sol, info);

        // Apply user-defined functor to solution vector
        sol = f(sol);

        //debug_jochen << "info.dimensions[0] = " << info.dimensions[0] << std::endl;
        //debug_jochen << "info.dimensions[1] = " << info.dimensions[1] << std::endl;

#if USE_PARALLEL == 1
        // Gather partial vectors on the root node
        const int rank = gf.getGridView().comm().rank();
        const int size = gf.getGridView().comm().size();

        // Commented out communication for every GF, rather use static information from Physics! (see below)
  //      // Gather x dimensions of all processors on root node
  //      std::vector<int> all_dimensions_x(size);
  //      gf.getGridView().comm().gather(&info.dimensions[0], &all_dimensions_x[0], 1, 0);
  //
  //      // Calculate displacements from all_dimensions_x
  //      std::vector<int> displacements(size);
  //      for(int i=1; i<displacements.size(); i++)
  //      {
  //        displacements[i] = displacements[i-1]+all_dimensions_x[i-1];
  //      }

        // Old version to determine nValuesPerElement, use new version provided below
  //      // This should work as long as info.dimensions[0] is always a multiple of the number of elements in x-direction on this process;
  //      // otherwise use the more general (but more expensive) version above!
  //      int nValuesPerElement = sol.size() / physics.nElementsForGF(gf);
  ////      debug_jochen << "  sol.size(): " << sol.size() << std::endl;
  ////      debug_jochen << "  physics.nElementsForGF(gf): " << physics.nElementsForGF(gf) << std::endl;

        int nValuesPerElement = info.nValuesPerElement[0] * info.nValuesPerElement[1];

        assert(nValuesPerElement * physics.nElementsForGF(gf) == sol.size()); // Make sure integer division leaves no remainder

        std::vector<int> all_dimensions_x = physics.nElements_AllProcessors(0); // copy, must not be const in MPI_Gatherv
        std::vector<int> displacements_x = physics.nOffset_AllProcessors(0); // copy, must not be const in MPI_Gatherv

        std::vector<int> all_dimensions(all_dimensions_x.size());
        std::vector<int> displacements(all_dimensions_x.size());
        for(int i=0; i<all_dimensions_x.size(); i++)
        {
          all_dimensions_x[i] *= info.nValuesPerElement[0];
          all_dimensions[i] = all_dimensions_x[i] * info.dimensions[1];
          //debug_jochen << "all_dimensions[" << i << "] = " << all_dimensions[i] << std::endl;

          displacements_x[i] *= info.nValuesPerElement[0];
          displacements[i] = (i == 0 ? 0 : displacements[i-1] + all_dimensions[i-1]);

          //debug_jochen << "displacements[" << i << "] = " << displacements[i] << std::endl;
        }

        int total_size_x = std::accumulate(all_dimensions_x.begin(),all_dimensions_x.end(),0);
        int total_size_y = info.dimensions[1]; // assume info.dimensions[1] to be the same on each process!
        int totalSize = total_size_x * total_size_y;

        //debug_jochen << "total_size_x = " << total_size_x << std::endl;
        //debug_jochen << "totalSize = " << totalSize << std::endl;

        std::vector<typename PHYSICS::Element::Geometry::GlobalCoordinate> pos_all(rank == 0 ? totalSize : 0);
        std::vector<typename GF::Traits::RangeType> sol_all(rank == 0 ? totalSize : 0);

        std::vector<typename PHYSICS::Element::Geometry::GlobalCoordinate> pos_all_temp(rank == 0 ? totalSize : 0);
        std::vector<typename GF::Traits::RangeType> sol_all_temp(rank == 0 ? totalSize : 0);

        // New communication method: Communicate all at once and sort everything later in order to avoid
        // lots of separate communication calls which conderably slow down performance on distributed systems
        // like Helics
        if(gf.getGridView().comm().size() > 1)
        {
          /* MPI_Gatherv arguments:
           * 1 - starting address of my send buffer
           * 2 - number of elements I am sending
           * 3 - type of elements I am sending
           * 4 - starting address of receive buffer (root node only)
           * 5 - array with number of elements each node is sending (root only)
           * 6 - array with displacements relative to receive buffer starting address where to place
           *     incoming data (root only)
           * 7 - type of elements received (root only)
           * 8 - root node (receiver)
           * 9 - communicator
           */
          int status = MPI_Gatherv(&pos[0], pos.size(),
              Dune::MPITraits<typename PHYSICS::Element::Geometry::GlobalCoordinate>::getType(),
              &pos_all_temp[0], &all_dimensions[0], &displacements[0],
              Dune::MPITraits<typename PHYSICS::Element::Geometry::GlobalCoordinate>::getType(),
              0, gf.getGridView().comm());

          status = MPI_Gatherv(&sol[0], sol.size(),
              Dune::MPITraits<typename GF::Traits::RangeType>::getType(),
              &sol_all_temp[0], &all_dimensions[0], &displacements[0],
              Dune::MPITraits<typename GF::Traits::RangeType>::getType(),
              0, gf.getGridView().comm());

          assert(status == 0);


          // Reorder entries (root node only)
          if(gf.getGridView().comm().rank() == 0)
          {
            int index = 0;
            for(int p=0; p<gf.getGridView().comm().size(); p++)
            {
              for(int i=0; i<info.dimensions[1]; i++)
              {
                for(int j=0; j<all_dimensions_x[p]; j++)
                {
                  // Running index (local within this processor block)
                  int running_index = j + i*all_dimensions_x[p];

                  // Don't ask me what this does. It does the right thing, that's all I know. I'm going to bed now.
                  int new_index = displacements_x[p] + i*(total_size_x - all_dimensions_x[p])
                      + running_index;

                  //debug_jochen << "index " << index << " -> " << new_index << std::endl;
                  pos_all[new_index] = pos_all_temp[index];
                  sol_all[new_index] = sol_all_temp[index];
                  index++;
                }

              }
            }
          }

        } else {
          pos_all = pos;
          sol_all = sol;
        }
#else
        const std::vector<typename PHYSICS::Element::Geometry::GlobalCoordinate>& pos_all(pos);
        const std::vector<typename GF::Traits::RangeType>& sol_all(sol);
#endif

        if(gf.getGridView().comm().rank() == 0)
        {
          // Note to self: Previously, here was a switch/case block which dynamically determined the GF output method at
          // runtime; however, the compiler still tries to generate code for each and every of the cases, for each gridfunction.
          // This is 1) unnecessarily increasing compile time and 2) will yield compiler errors of the GF does not comply with
          // the interface expected by the different output methods.
          // The new way to handle this is via a TMP, which chooses the appropriate (and only that one!) method to be compiled
          // for each GF separately.
          std::string prefix("");
          strategy<PHYSICS,GF,OutputFunction<OutputTraits,GF>::value>::write(gf,filename,time,pos_all,sol_all,prefix,infoStr);
        }
        debug_jochen << "    - Output time: " << timer.elapsed() << " s"
            << " (getSolutionVector part: " << info.tElapsed << "s)" << std::endl;
      }
    }

    //! struct for choosing the matching gnuplot output function at compile-time; default version throws an exception,
    //! has to be spcialized for each gnuplot output method
    template<typename PHYSICS, typename DGF, int>
    struct strategy
    {
      static void write(const DGF& gf,
          const std::string& filename,
          const double time,
          const std::vector<typename PHYSICS::Element::Geometry::GlobalCoordinate>& pos,
          const std::vector<typename DGF::Traits::RangeType>& sol,
          const std::string& prefix,
          const std::string& infoString)
      {
        DUNE_THROW(Dune::Exception, "Could not find an output strategy for GF type "
                  << Tools::getTypeName(gf));
      }

    };

    template<typename PHYSICS, typename DGF>
    struct strategy<PHYSICS,DGF,GnuplotFunctions::GNUPLOT_ADD_BLOCK>
    {
      static void write(const DGF& gf,
          const std::string& filename,
          const double time,
          const std::vector<typename PHYSICS::Element::Geometry::GlobalCoordinate>& pos,
          const std::vector<typename DGF::Traits::RangeType>& sol,
          const std::string& prefix,
          const std::string& infoString)
      {
        //debug_jochen << "GNUPLOT_ADD_BLOCK" << std::endl;
        GnuplotTools2D<OutputTraits>::gnuplotAddBlock(filename,pos,sol,prefix,infoString);
      }
    };

    template<typename PHYSICS, typename DGF>
    struct strategy<PHYSICS,DGF,GnuplotFunctions::GNUPLOT_ADD_LINE>
    {
      static void write(const DGF& gf,
          const std::string& filename,
          const double time,
          const std::vector<typename PHYSICS::Element::Geometry::GlobalCoordinate>& pos,
          const std::vector<typename DGF::Traits::RangeType>& sol,
          const std::string& prefix,
          const std::string& infoString)
      {
        //debug_jochen << "GNUPLOT_ADD_LINE" << std::endl;
        GnuplotTools2D<OutputTraits>::gnuplotAddLine(filename,time,sol,prefix,infoString);
      }
    };

    template<typename PHYSICS, typename DGF>
    struct strategy<PHYSICS,DGF,GnuplotFunctions::GNUPLOT_ADD_MULTI_COLUMN_BLOCK>
    {
      static void write(const DGF& gf,
          const std::string& filename,
          const double time,
          const std::vector<typename PHYSICS::Element::Geometry::GlobalCoordinate>& pos,
          const std::vector<typename DGF::Traits::RangeType>& sol,
          const std::string& prefix,
          const std::string& infoString)
      {
        //debug_jochen << "GNUPLOT_ADD_MULTI_COLUMN_BLOCK" << std::endl;
        //GnuplotTools2D<OutputTraits>::gnuplotAddMultiColumnBlock(filename,pos,sol,prefix,infoString);
        DUNE_THROW(Dune::NotImplemented, "Huhn-Hahn!");
      }
    };


    template<typename DT, typename RT>
    static void gnuplotAddBlock(const std::string& filename,
                                 const std::vector<DT>& x,
                                 const std::vector<RT>& y,
                                 const std::string& prefix = "",
                                 const std::string& infoString = "",
                                 bool insertNewlines = true)
    {
      std::ofstream file;
      file.open((prefix + filename).c_str(), std::ofstream::app);
      if (! file.good())
        DUNE_THROW(Dune::IOError, "Could not write to file " << filename);

      file << "# " << infoString << std::endl;

      file << std::setprecision(prec) << std::scientific;

      DT x0 = x[0];
      for(int i=0; i<x.size(); ++i)
      {
        file << x[i] << " " << y[i] << std::endl;

        if(i<x.size()-1)
        {
          DT xDiff = x[i+1];
          xDiff -= x[i];

          // Insert newline delimiter in case of different y values => new subblock
          // Second criterion: x is x0 => new subblock
          // (useful for distinct block with identical y values)
          if (insertNewlines && (std::abs(xDiff[1]) > eps || std::abs(x0[0]-x[i+1][0]) < eps) )
          {
            file << std::endl;
          }
        }

      }

      // Finalize this block with two newlines
      file << std::endl << std::endl;
      file.close();
    }

    template<typename DT, typename... TupleElements>
    static void gnuplotAddBlock(const std::string& filename,
                                 const std::vector<DT>& x,
                                 const std::vector<std::tuple<TupleElements...> >& y,
                                 const std::string& prefix = "",
                                 const std::string& infoString = "",
                                 bool insertNewlines = true)
    {
      std::ofstream file;
      file.open((prefix + filename).c_str(), std::ofstream::app);
      if (! file.good())
        DUNE_THROW(Dune::IOError, "Could not write to file " << filename);

      file << "# " << infoString << std::endl;

      file << std::setprecision(prec) << std::scientific;

      DT x0 = x[0];
      for(int i=0; i<x.size(); ++i)
      {
        // coordinate
        file << x[i] << " ";

        // Print tuple
        typedef std::tuple<TupleElements...> Tuple;
        Dune::ForLoop<write_to_ostream,0,std::tuple_size<Tuple>::value-1>::apply(file, y[i]);

        // newline
        file << std::endl;

        if(i<x.size()-1)
        {
          DT xDiff = x[i+1];
          xDiff -= x[i];

          // Insert newline delimiter in case of different y values => new subblock
          // Second criterion: x is x0 => new subblock
          // (useful for distinct block with identical y values)
          if (insertNewlines && (std::abs(xDiff[1]) > eps || std::abs(x0[0]-x[i+1][0]) < eps) )
          {
            file << std::endl;
          }
        }

      }

      // Finalize this block with two newlines
      file << std::endl << std::endl;
      file.close();
    }

    template<typename DT, typename RT>
    static void gnuplotAddMultiColumnBlock(const std::string& filename,
                                 const std::vector<DT>& x,
                                 const std::vector<std::vector<RT> >& y,
                                 const std::string& prefix = "",
                                 const std::string& infoString = "",
                                 bool insertNewlines = true)
    {
      std::ofstream file;
      file.open((prefix + filename).c_str(), std::ofstream::app);
      if (! file.good())
        DUNE_THROW(Dune::IOError, "Could not write to file " << filename);

      file << "# " << infoString << std::endl;

      file << std::setprecision(prec) << std::scientific;

      DT x0 = x[0];
      for(int i=0; i<x.size(); ++i)
      {
        file << x[i];
        for(int j=0; j<y.size(); j++)
        {
          file << " " << y[j][i];
        }
        file << std::endl;

        if(i<x.size()-1)
        {
          DT xDiff = x[i+1];
          xDiff -= x[i];

          // Insert newline delimiter in case of different y values => new subblock
          // Second criterion: x is x0 => new subblock
          // (useful for distinct block with identical y values)
          if (insertNewlines && (std::abs(xDiff[1]) > eps || std::abs(x0[0]-x[i+1][0]) < eps) )
          {
            file << std::endl;
          }
        }
      }

      // Finalize this block with two newlines
      file << std::endl << std::endl;
      file.close();
    }

    template<typename DT, typename RT>
    static void gnuplotAddLine(const std::string& filename,
                                 const DT& x,
                                 const std::vector<RT>& y,
                                 const std::string& prefix = "",
                                 const std::string& infoString = "")
    {
      std::ofstream file;
      file.open((prefix + filename).c_str(), std::ofstream::app);
      if (! file.good())
        DUNE_THROW(Dune::IOError, "Could not write to file " << filename);

      if(infoString != "")
      {
        file << "# " << infoString << std::endl;
      }

      file << std::setprecision(prec) << std::scientific << x << " ";
      for(int i=0; i<y.size(); ++i)
      {
        file << y[i] << " ";
      }
      file << std::endl;

      file.close();
    }


    template<typename DT, typename RT>
    static void gnuplotAddLine(const std::string& filename,
                                 const DT& x,
                                 const RT& y,
                                 const std::string& prefix = "",
                                 const std::string& infoString = "")
    {
      std::ofstream file;
      file.open((prefix + filename).c_str(), std::ofstream::app);
      if (! file.good())
        DUNE_THROW(Dune::IOError, "Could not write to file " << filename);

      if(infoString != "")
      {
        file << "# " << infoString << std::endl;
      }

      file << std::setprecision(prec) << std::scientific << x << " " << y << std::endl;

      file.close();
    }

    template<int i>
    struct write_to_ostream
    {
      template<typename Tuple>
      static void apply(std::ostream& stream, const Tuple& tuple)
      {
        stream << std::get<i>(tuple) << " ";
      }

    };

};



#endif /* DUNE_AX1_GNUPLOT_TOOLS_HH */
