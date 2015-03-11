/*
 * ax1_lfs_tools.hh
 *
 *  Created on: Apr 4, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_LFS_TOOLS_HH
#define DUNE_AX1_LFS_TOOLS_HH

class CoordinateEvaluation
  {
  public:
    // store the coordinate to evaluate
    CoordinateEvaluation (int i_) : i(i_) {}

    // eval coordinate i
    template<typename DT, typename RT>
    inline void evaluate (const DT& x, RT& y) const
    {
      y = x[i];
      return;
    }

  private:
    int i;
  };



class Ax1LFSTools
{
  public:
    template<typename V, typename FV>
    static int getNodeIndex(const V& v, const FV& pos)
    {
      for(int i=0; i<v.size(); i++)
      {
        if(Ax1LFSTools::fieldVectorEqual(v[i],pos))
          return i;
      }
      return -1;
    }

    template<typename FV>
    static bool fieldVectorEqual(const FV& v1, const FV& v2)
    {
      bool equal = true;
      for(int i=0; i<v1.size(); i++)
      {
        equal = equal && (std::abs(v1[i] -v2[i]) < 1e-6);
      }
      return equal;
    }

};

struct field_vector_equal
{
  template<typename FV>
  inline bool operator() (const FV& v1, const FV& v2) const
  {
    bool equal = true;
    for(int i=0; i<v1.size(); i++)
    {
      equal = equal && (std::abs(v1[i] -v2[i]) < 1e-6);
    }
    return equal;
  }
};

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

//! Create std::hash specialization for Dune::FieldVector
namespace std
{
  template<typename RF, int dim> struct hash<Dune::FieldVector<RF,dim> >
  {
    inline size_t operator()(const Dune::FieldVector<RF,dim> & v) const
    {
      size_t seed = 0;
      for(int i=0; i<v.size(); i++)
        ::hash_combine(seed, v[i]);

      return seed;
    }
  };
}


#endif /* DUNE_AX1_LFS_TOOLS_HH */
