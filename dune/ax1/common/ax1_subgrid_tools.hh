/*
 * ax1_subgrid_tools.hh
 *
 *  Created on: Jan 27, 2012
 *      Author: jpods
 */

#ifndef DUNE_AX1_SUBGRID_TOOLS_HH
#define DUNE_AX1_SUBGRID_TOOLS_HH

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

//! \todo Please doc me !
template<typename GV, typename SubGV, typename FULL_GFS, typename SUB_GFS, typename UFULL, typename USUB>
class SubGridCoefficientVectorRestrictor
{
  public:
    typedef typename SubGV::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GV::template Codim<0>::EntityPointer HostElementPointer;


    /** \brief Restrict hostgrid coefficient vector to subgrid coefficient vector
     * The underlying local function spaces are assumed to have identical DOF orderings per entity!
     *
     * Also works if FULL_GFS and SUB_GFS are PowerGFS's!
     *
     * This is a shortcut for creating an (temporary) restrictor object and
     * handing it to the transfer method in the subgrid
     */
    static void restrict(const GV& gv_, const SubGV& subGV_, const FULL_GFS& fullGfs_,
        const SUB_GFS& subGfs_, const UFULL& uFull_, USUB& uSub_)
    {
      dune_static_assert(FULL_GFS::CHILDREN == SUB_GFS::CHILDREN,
          "Power grid function spaces must have the same number of children!");

      typedef SubGridCoefficientVectorRestrictor<GV,SubGV,FULL_GFS,SUB_GFS,UFULL,USUB> MyType;

      MyType restrictor(gv_, subGV_, fullGfs_, subGfs_, uFull_, uSub_);
      subGV_.grid().template transfer<MyType> (restrictor);
    }

    //! \todo Please doc me !
    SubGridCoefficientVectorRestrictor(const GV& gv_, const SubGV& subGV_, const FULL_GFS& fullGfs_,
        const SUB_GFS& subGfs_, const UFULL& uFull_, USUB& uSub_)
    : gv(gv_),
      subGV(subGV_),
      fullGfs(fullGfs_),
      subGfs(subGfs_),
      uFull(uFull_),
      uSub(uSub_)
    {}

    //! Prepare subgrid vector
    void pre()
    {
      //uSub.resize(grid.size(dim));
      uSub = 0.0;
    }

    //! \todo Please doc me !
    void transfer(const ElementPointer& element, const HostElementPointer& hostElement)
    {
      //typedef Dune::FiniteElementInterfaceSwitch<
      //    typename Dune::PDELab::LocalFunctionSpace<FULL_GFS>::Traits::FiniteElementType> FESwitch;

      // TODO get real RangeFieldType
      typedef double RangeFieldType;
      std::vector<RangeFieldType> localContainer(fullGfs.maxLocalSize());

      Dune::PDELab::LocalFunctionSpace<FULL_GFS> fullLfs(fullGfs);
      //Bind to local function space of host entity
      fullLfs.bind(*hostElement);

      // Read local coefficient vector from global coefficient vector
      fullLfs.vread(uFull,localContainer);

      //debug_verb << "= Host entity @" << hostElement->geometry().center() << std::endl;
      //for(int i=0; i<localContainer.size(); i++)
      //{
      //  debug_verb << " u_local[" << i << "] = " << localContainer[i] << std::endl;
      //}

      // Now transfer the coefficients; assume that numbering is identical for both local function spaces!
      Dune::PDELab::LocalFunctionSpace<SUB_GFS> subLfs(subGfs);
      //Bind to local function space of sub entity
      subLfs.bind(*element);

      // Fill global coefficient vector by reading from local coefficient vector
      subLfs.vwrite(localContainer,uSub);
    }

    //! does nothing
    void post()
    {
    }

  private:

    const GV& gv;
    const SubGV& subGV;

    const FULL_GFS& fullGfs;
    const SUB_GFS& subGfs;

    const UFULL& uFull;
    USUB& uSub;
};


//! \todo Please doc me !
template<typename SubGV, typename GV, typename SUB_GFS, typename FULL_GFS, typename USUB, typename UFULL>
class SubGridCoefficientVectorInterpolator
{
  public:
    typedef typename SubGV::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GV::template Codim<0>::EntityPointer HostElementPointer;

    /** \brief Interpolate subgrid coefficient vector onto host grid coefficient vector
     * The underlying local function spaces are assumed to have identical DOF orderings per entity!
     *
     * Also works if FULL_GFS and SUB_GFS are PowerGFS's!
     *
     * This is a shortcut for creating an (temporary) interpolator object and
     * handing it to the transfer method in the subgrid
     */
    static void interpolate(const SubGV& subGV_, const GV& gv_, const SUB_GFS& subGfs_,
        const FULL_GFS& fullGfs_, const USUB& uSub_, UFULL& uFull_)
    {
      dune_static_assert(FULL_GFS::CHILDREN == SUB_GFS::CHILDREN,
          "Power grid function spaces must have the same number of children!");

      typedef SubGridCoefficientVectorInterpolator<SubGV,GV,SUB_GFS,FULL_GFS,USUB,UFULL> MyType;

      MyType interpolator(subGV_, gv_, subGfs_, fullGfs_, uSub_, uFull_);
      subGV_.grid().template transfer<MyType> (interpolator);
    }

    //! \todo Please doc me !
    SubGridCoefficientVectorInterpolator(const SubGV& subGV_, const GV& gv_, const SUB_GFS& subGfs_,
        const FULL_GFS& fullGfs_, const USUB& uSub_, UFULL& uFull_)
    : subGV(subGV_),
      gv(gv_),
      subGfs(subGfs_),
      fullGfs(fullGfs_),
      uSub(uSub_),
      uFull(uFull_)
    {}

    //! Does nothing
    void pre()
    {
    }

    //! \todo Please doc me !
    void transfer(const ElementPointer& element, const HostElementPointer& hostElement)
    {
      //typedef Dune::FiniteElementInterfaceSwitch<
      //    typename Dune::PDELab::LocalFunctionSpace<FULL_GFS>::Traits::FiniteElementType> FESwitch;

      // TODO get real RangeFieldType
      typedef double RangeFieldType;
      std::vector<RangeFieldType> localContainer(subGfs.maxLocalSize());

      Dune::PDELab::LocalFunctionSpace<SUB_GFS> subLfs(subGfs);
      //Bind to local function space of subgrid entity
      subLfs.bind(*element);

      // Read local coefficient vector from global coefficient vector
      subLfs.vread(uSub,localContainer);

      //debug_verb << "= Host entity @" << hostElement->geometry().center() << std::endl;
      //for(int i=0; i<localContainer.size(); i++)
      //{
      //  debug_verb << " u_local[" << i << "] = " << localContainer[i] << std::endl;
      //}

      // Now transfer the coefficients; assume that numbering is identical for both local function spaces!
      Dune::PDELab::LocalFunctionSpace<FULL_GFS> fullLfs(fullGfs);
      //Bind to local function space of host grid entity
      fullLfs.bind(*hostElement);

      // Fill global coefficient vector by reading from local coefficient vector
      fullLfs.vwrite(localContainer,uFull);
    }

    //! does nothing
    void post()
    {
    }

  private:

    const SubGV& subGV;
    const GV& gv;

    const SUB_GFS& subGfs;
    const FULL_GFS& fullGfs;

    const USUB& uSub;
    UFULL& uFull;
};

#endif /* DUNE_AX1_SUBGRID_TOOLS_HH */
