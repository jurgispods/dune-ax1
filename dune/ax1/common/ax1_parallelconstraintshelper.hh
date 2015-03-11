/*
 * ax1_parallelconstraintshelper.hh
 *
 *  Created on: Dec 3, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_PARALLELCONSTRAINTSHELPER_HH
#define DUNE_AX1_PARALLELCONSTRAINTSHELPER_HH

//! \brief This class is derived from OverlappingConformingDirichletConstraints and overwrites the method
// processor(), which is important in the parallel case. It contains a flag 'overlapIsDirichlet', which by
// default is true and simply calls the father class method. If it is set to false, it does nothing, which
// means that a contraints containers constructed from this class will NOT consider overlap entities to
// be Dirichlet constraints, which is useful for a very special case when interpolating time-dependent
// Dirichlet boundary values in parallel.
class Ax1OverlappingDirichletContraints : public Dune::PDELab::OverlappingConformingDirichletConstraints
{
  public:
    typedef Dune::PDELab::OverlappingConformingDirichletConstraints BaseT;

    Ax1OverlappingDirichletContraints(bool overlapIsDirichlet_ = true)
    : overlapIsDirichlet(overlapIsDirichlet_)
    {}

    enum { doProcessor = true };

    //! processor constraints
    /**
     * \tparam IG  intersection geometry
     * \tparam LFS local function space
     * \tparam T   TransformationType
     */
    template<typename IG, typename LFS, typename T>
    void processor (const IG& ig, const LFS& lfs, T& trafo) const
    {
      if(overlapIsDirichlet)
        BaseT::processor(ig, lfs, trafo);

      // Do nothing when flag is not set (=> overlap is not a considered Dirichlet boundary)
    }

    void setOverlapIsDirichlet(bool overlapIsDirichlet_)
    {
      overlapIsDirichlet =  overlapIsDirichlet_;
    }

  private:
    bool overlapIsDirichlet;
};

#endif /* DUNE_AX1_PARALLELCONSTRAINTSHELPER_HH */
