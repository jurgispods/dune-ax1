/*
 * ax1_parallelhelper.hh
 *
 *  Created on: May 5, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_PARALLELHELPER_HH
#define DUNE_AX1_PARALLELHELPER_HH

#include <dune/common/forloop.hh>

#if HAVE_MPI

namespace Dune {


  template<typename... Elements>
  struct MPITraits<std::tuple<Elements...> >
  {
    typedef std::tuple<Elements...> TupleType;

    static MPI_Datatype datatype;

    static inline MPI_Datatype getType()
    {
     if(datatype==MPI_DATATYPE_NULL) {
       const int n = sizeof...(Elements);

       MPI_Datatype oldtypes[n];
       int blockcounts[n];
       MPI_Aint offsets[n], extent;

       // Initialize data types for elements i=1,...,n-1
       Dune::ForLoop<create_tuple_type, 0, n-1>::apply(offsets, blockcounts, oldtypes);

//       std::cout << std::endl;
//       std::cout << "--------------" << std::endl;
//       for(int i=0; i<n; i++)
//       {
//         std::cout << "offsets[" << i << "] = " << offsets[i] << std::endl;
//         std::cout << "blockcounts[" << i << "] = " << blockcounts[i] << std::endl;
//       }
//       std::cout << "--------------" << std::endl;
//       std::cout << std::endl;

       MPI_Type_struct(n, blockcounts, offsets, oldtypes, &datatype);
       MPI_Type_commit(&datatype);
     }

     return datatype;
    }



    template<int i>
    struct create_tuple_type {

      static void apply(MPI_Aint* offsets, int* blockcounts, MPI_Datatype* oldtypes)
      {
        TupleType tuple;
        MPI_Aint base;
        MPI_Aint displ;
        MPI_Address(&tuple, &base);
        MPI_Address(&(std::get<i>(tuple)), &displ);
        displ -= base;

        offsets[i] = displ;
        //offsets[i] = offsets[i-1] + blockcounts[i-1] * previous_extent;
        oldtypes[i] = MPITraits<typename std::tuple_element<i,TupleType>::type>::getType();
        blockcounts[i] = 1;

//        std::cout << "offsets[" << i << "] = " << offsets[i] << std::endl;
//        std::cout << "oldtypes[" << i << "] = " << oldtypes[i] << std::endl;
//        std::cout << "blockcounts[" << i << "] = " << blockcounts[i] << std::endl;
      }
    };

  };

  template<typename... Elements>
  MPI_Datatype MPITraits<std::tuple<Elements...> >::datatype = MPI_DATATYPE_NULL;

}

#endif



class Ax1LogTag
{
public:
  Ax1LogTag(const int rank_)
  : rank(rank_)
  {}

  std::ostream& operator() (std::ostream &s) const
  {
    return s << "[p" << rank << "]";
  }

private:
  const int rank;
};

template<typename DebugStream>
class Ax1ParallelDebugStream : public DebugStream
{
public:
  typedef Ax1ParallelDebugStream<DebugStream> This;

  Ax1ParallelDebugStream(std::ostream& out = std::cerr)
  : DebugStream(out),
    logtag(""),
    insertTag(true)
  {}

  Ax1ParallelDebugStream(Dune::DebugStreamState& master,
      std::ostream& fallback = std::cerr, std::string logtag_ = "")
  : DebugStream(master, fallback),
    logtag(logtag_),
    insertTag(true)
  {}

  //! \brief Generic types are passed on to current output stream
  template <class T>
  This& operator<<(const T data)
  {
    if(insertTag)
    {
      insertTag = false;
      return (This&) DebugStream::operator<<(logtag) << data;
    }
    else return (This&) DebugStream::operator<<(data);

    //return DebugStream::operator<<(Dune::PDELab::logtag) << data;
  }

  This& operator<<(const int data)
  {
    if(insertTag)
    {
      insertTag = false;
      return (This&) DebugStream::operator<<(logtag) << data;
    } else
      return (This&) DebugStream::operator<<(data);

    insertTag = false;
  }

  This& operator<<(std::ostream& (*f)(std::ostream&))
  {
    typedef std::ostream& (*manip_t)(std::ostream&);

    //if(f == std::endl)
    if(f == static_cast<manip_t>(std::endl))
    {
      //std::cerr << "std::endl detected!" << std::endl;
      insertTag = true;
    }

    return (This&) DebugStream::operator<<(f);
  }

  void setLogTag(std::string tag)
  {
    logtag = tag;
  }

private:
 std::string logtag;
 bool insertTag;

};


#endif /* AX1_PARALLELHELPER_HH_ */
