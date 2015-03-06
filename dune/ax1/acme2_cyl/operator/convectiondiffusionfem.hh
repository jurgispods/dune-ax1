// -*- tab-width: 4; indent-tabs-mode: nil -*-

/*
 * This file was copied from <dune/pdelab/localoperator/convectiondiffusionfem.hh>
 * so that it can be modified for testing purposes without affecting the original one
 */
#ifndef DUNE_PDELAB_CONVECTIONDIFFUSIONFEM_HH
#define DUNE_PDELAB_CONVECTIONDIFFUSIONFEM_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/type.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include<dune/pdelab/localoperator/callswitch.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <dune/ax1/common/ax1_lfs_tools.hh>

namespace Dune {
  namespace PDELab {

    /** a local operator for solving convection-diffusion equation with standard FEM
     *  
     * \f{align*}{
     *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \\
     *                                              u &=& g \mbox{ on } \partial\Omega_D \\
     *                (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
     *                        -(A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_O
     * \f}
     * Note:
     *  - This formulation is valid for velocity fields which are non-divergence free.
     *  - Outflow boundary conditions should only be set on the outflow boundary
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T, typename FiniteElementMap, bool useMembraneContributions=true>
    class ConvectionDiffusionFEM : 
      public Dune::PDELab::NumericalJacobianApplyVolume<ConvectionDiffusionFEM<T,FiniteElementMap,useMembraneContributions> >,
      public Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionFEM<T,FiniteElementMap,useMembraneContributions> >,
      public Dune::PDELab::NumericalJacobianVolume<ConvectionDiffusionFEM<T,FiniteElementMap,useMembraneContributions> >,
      public Dune::PDELab::NumericalJacobianBoundary<ConvectionDiffusionFEM<T,FiniteElementMap,useMembraneContributions> >,
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags,
      public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>
    {
      enum { dim = T::Traits::GridViewType::dimension };
      typedef typename T::Traits::RangeFieldType Real;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doAlphaBoundary = true };

      ConvectionDiffusionFEM (T& param_, const int lfsPotSize, int intorderadd_=0)
        : param(param_),
          intorderadd(intorderadd_),
          POTENTIAL_RESIDUAL_SCALE(param.getPhysics().getParams().getPotentialResidualScalingFactor()),
          lfsuPotSize(lfsPotSize),
          doVolumeScaling(param.getPhysics().getParams().doVolumeScaling()),
          volumeScale(lfsuPotSize, 1.0),
          local_position(lfsuPotSize),
          local_volume(lfsuPotSize),
          c(lfsuPotSize),
          useJacobianVolume(param.getPhysics().getParams().general.get("useMembOperatorJacobianVolume",false))
      {
        debug_jochen << "ConvectionDiffusionFEM(), useMembraneContributions="
            << useMembraneContributions << std::endl;
      }

      // volume integral depending on test and ansatz functions
      template<typename EG_ORIG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG_ORIG& eg_orig, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        if(!useMembraneContributions && param.getPhysics().isMembrane(eg_orig.entity()))
        {
//          debug_jochen << "Skipping intersection coming from membrane element!" << std::endl;
          return;
        }

//        debug_jochen << "ConvectionDiffusionFEM::alpha_volume @ "<< eg_orig.geometry().center() << std::endl;
//        debug_jochen << "  useMembraneContributions: " << useMembraneContributions << std::endl;
        typedef typename Acme2CylGeometrySwitch::ElementGeometrySwitch<EG_ORIG>::type EG;
        const EG& eg(eg_orig);

        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // ============================================================================================
        // Determine Lagrange points for each basis (=test) function
        // => get node position the test functions belong to
        if(doVolumeScaling)
        {
          for (int k=0; k<dim; k++)
          {
            CoordinateEvaluation f(k);
            lfsu.finiteElement().localInterpolation().interpolate(f,c);
            for (size_type i=0; i<lfsu.size(); i++)
              local_position[i][k] = c[i];
          }
          for (size_type i=0; i<lfsu.size(); i++)
          {
            local_volume[i] = param.getPhysics().getNodeVolume(eg.geometry().global(local_position[i]));
            volumeScale[i] = 1./ local_volume[i];
          }
        }
        // ============================================================================================

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const int intorder = intorderadd+2*lfsu.finiteElement().localBasis().order();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
#if USE_CACHE
            const std::vector<RangeType>& phi = cache.evaluateFunction(it->position(),lfsu.finiteElement().localBasis());
#else
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);
#endif


            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // evaluate diffusion tensor at quadrature point
            typename T::Traits::PermTensorType tensor;
            tensor = param.A(eg.entity(),it->position());


            //debug_verb << "[ConvectionDiffusionFEM] u = " << u << std::endl;

            //Dune::ios_base_all_saver atzbert(std::cout);
            //if(physics.getElementIndex(eg.entity()) == 0)
            //{
            //  debug_jochen << "u  @ " << eg.entity().geometry().global(it->position()) << " : "
            //      << std::scientific << std::setprecision(16)<< u << std::endl;
            //}

            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
#if USE_CACHE
            const std::vector<JacobianType>& js = cache.evaluateJacobian(it->position(),lfsu.finiteElement().localBasis());
#else
            std::vector<JacobianType> js(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);
#endif


            // transform gradients of shape functions to real element
            const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            for (size_type i=0; i<lfsu.size(); i++)
              jac.mv(js[i][0],gradphi[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (size_type i=0; i<lfsu.size(); i++)
              gradu.axpy(x(lfsu,i),gradphi[i]);

            //debug_verb << "[ConvectionDiffusionFEM] gradu = " << gradu << std::endl;

            // compute A * gradient of u
            Dune::FieldVector<RF,dim> Agradu(0.0);
            tensor.umv(gradu,Agradu);

            // evaluate velocity field, sink term and source te
            typename T::Traits::RangeType b = param.b(eg.entity(),it->position());
            typename T::Traits::RangeFieldType c = param.c(eg.entity(),it->position());

            typename T::Traits::RangeFieldType f = param.f(eg.entity(),it->position());
            //typename T::Traits::RangeFieldType f = param.f(eg.entity(),it->position(),u);

            // integrate (A grad u)*grad phi_i - u b*grad phi_i + c*u*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            //debug_verb << "acc2" << std::endl;
            for (size_type i=0; i<lfsu.size(); i++)
            {
//              debug_jochen << "[LOP] pot @" << lfsu.dofIndex(i)
//                  << (( Agradu*gradphi[i] - u*(b*gradphi[i]) + (c*u-f)*phi[i] )
//                  *factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE) << std::endl;
              r.accumulate(lfsu,i,( Agradu*gradphi[i] - u*(b*gradphi[i]) + (c*u-f)*phi[i] )
                  *factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE);
            }
          }
      }

#ifdef AX1_USE_JACOBIAN_METHODS_IN_OPERATOR
      // jacobian of volume term
      template<typename EG_ORIG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG_ORIG& eg_orig, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        if(! useJacobianVolume)
        {
          Dune::PDELab::NumericalJacobianVolume<ConvectionDiffusionFEM<T,FiniteElementMap,
            useMembraneContributions> >::jacobian_volume(eg_orig, lfsu, x, lfsv,mat);
          return;
        }

        typedef typename Acme2CylGeometrySwitch::ElementGeometrySwitch<EG_ORIG>::type EG;
        const EG& eg(eg_orig);

        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // ============================================================================================
        // Determine Lagrange points for each basis (=test) function
        // => get node position the test functions belong to
        if(doVolumeScaling)
        {
          for (int k=0; k<dim; k++)
          {
            CoordinateEvaluation f(k);
            lfsu.finiteElement().localInterpolation().interpolate(f,c);
            for (size_type i=0; i<lfsu.size(); i++)
              local_position[i][k] = c[i];
          }
          for (size_type i=0; i<lfsu.size(); i++)
          {
            local_volume[i] = param.getPhysics().getNodeVolume(eg.geometry().global(local_position[i]));
            volumeScale[i] = 1./ local_volume[i];
          }
        }
        // ============================================================================================

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const int intorder = intorderadd+2*lfsu.finiteElement().localBasis().order();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        typename T::Traits::PermTensorType tensor;
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
        tensor = param.A(eg.entity(),localcenter);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        {
          // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
#if USE_CACHE
          const std::vector<JacobianType>& js = cache.evaluateJacobian(it->position(),lfsu.finiteElement().localBasis());
#else
          std::vector<JacobianType> js(lfsu.size());
          lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);
#endif

          // transform gradient to real element
          const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
          std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
          std::vector<Dune::FieldVector<RF,dim> > Agradphi(lfsu.size());
          for (size_type i=0; i<lfsu.size(); i++)
          {
            jac.mv(js[i][0],gradphi[i]);
            tensor.mv(gradphi[i],Agradphi[i]);
          }

          // evaluate basis functions
#if USE_CACHE
          const std::vector<RangeType>& phi = cache.evaluateFunction(it->position(),lfsu.finiteElement().localBasis());
#else
          std::vector<RangeType> phi(lfsu.size());
          lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);
#endif

          // evaluate velocity field, sink term and source term
          typename T::Traits::RangeType b = param.b(eg.entity(),it->position());
          typename T::Traits::RangeFieldType c = param.c(eg.entity(),it->position());

          // integrate (A grad phi_j)*grad phi_i - phi_j b*grad phi_i + c*phi_j*phi_i
          RF factor = it->weight() * eg.geometry().integrationElement(it->position());
          for (size_type j=0; j<lfsu.size(); j++)
          {
            for (size_type i=0; i<lfsu.size(); i++)
            {
              mat.accumulate(lfsu,i,lfsu,j,( Agradphi[j]*gradphi[i]-phi[j]*(b*gradphi[i])+c*phi[j]*phi[i] )
                  *factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE);
            }
          }
        }
      }
#endif

      // boundary integral
      template<typename IG_ORIG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG_ORIG& ig_orig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        if(!useMembraneContributions && param.getPhysics().isMembraneInterface(ig_orig.intersection())
            && param.getPhysics().isMembrane(*ig_orig.inside()))
        {
          //debug_jochen << "Skipping intersection coming from membrane element!" << std::endl;
          return;
        }

//        debug_jochen << "ConvectionDiffusionFEM::alpha_boundary @ "<< ig_orig.geometry().center() << std::endl;
//        debug_jochen << "  Inside element: "<< ig_orig.inside()->geometry().center() << std::endl;
//        debug_jochen << "  isMembrane: " << param.getPhysics().isMembrane(*ig_orig.inside()) << std::endl;
//        debug_jochen << "  useMembraneContributions: " << useMembraneContributions << std::endl;

        typedef typename Acme2CylGeometrySwitch::IntersectionGeometrySwitch<IG_ORIG>::type IG;
        const IG& ig(ig_orig);

        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        
        // evaluate boundary condition type
        Dune::GeometryType gtface = ig.geometryInInside().type();
        Dune::FieldVector<DF,dim-1> facecenterlocal = Dune::ReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        ConvectionDiffusionBoundaryConditions::Type bctype;
        bctype = param.bctype(ig.intersection(),facecenterlocal);

//        debug_jochen << "  bctype: " << bctype << std::endl;
 
        // skip rest if we are on Dirichlet boundary
        if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet) return;

        // Special case for dune-multidomain
        if (bctype==ConvectionDiffusionBoundaryConditions::None) return;

        // ============================================================================================
        // Determine Lagrange points for each basis (=test) function
        // => get node position the test functions belong to
        if(doVolumeScaling)
        {
          for (int k=0; k<dim; k++)
          {
            CoordinateEvaluation f(k);
            lfsu_s.finiteElement().localInterpolation().interpolate(f,c);
            for (size_type i=0; i<lfsu_s.size(); i++)
              local_position[i][k] = c[i];
          }
          for (size_type i=0; i<lfsu_s.size(); i++)
          {
            local_volume[i] = param.getPhysics().getNodeVolume(
                (*ig.inside()).geometry().global(local_position[i]) );
            volumeScale[i] = 1./local_volume[i];
          }
        }
        // ============================================================================================


        // select quadrature rule
        const int intorder = intorderadd+2*lfsu_s.finiteElement().localBasis().order();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

            // evaluate shape functions (assume Galerkin method) 
#if USE_CACHE
            const std::vector<RangeType>& phi = cache.evaluateFunction(local,lfsu_s.finiteElement().localBasis());
#else
            std::vector<RangeType> phi(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(local,phi);
#endif

            if (bctype==ConvectionDiffusionBoundaryConditions::Neumann)
              {
                // evaluate flux boundary condition
                typename T::Traits::RangeFieldType j = param.j(ig.intersection(),it->position());
            
                // integrate j
                RF factor = it->weight()*ig.geometry().integrationElement(it->position());

                // TODO Dirty hack, remove this
#ifdef USE_LAPLACE_OPERATOR
                // Switch off cylinder geometry transformation on the lower boundary
                if(param.doCollapseToLineSource()
                    && ig.geometry().global(it->position())[1]-1e-6 < param.getPhysics().getParams().yMin())
                {
                  factor = it->weight()*ig_orig.geometry().integrationElement(it->position());
                  if(param.getTime() < 20.0)
                    debug_jochen << "Switching off geometry transformation, factor =" << factor << std::endl;
                }
#endif

                for (size_type i=0; i<lfsu_s.size(); i++)
                {
//                  if(ig.geometry().global(it->position())[0] < 100e3
//                      && std::abs(ig.geometry().global(it->position())[1] - 505) < 1e-6)
//                  {
//                    debug_verb << "[LOP] pot @" << lfsu_s.dofIndex(i) << ": "
//                      << std::scientific << std::setprecision(16)
//                      << ", phi[i] " << phi[i]
//                      << ", factor " << factor
//                      << ", volumeScale[i] " << volumeScale[i]
//                      << ", POTENTIAL_RESIDUAL_SCALE" << POTENTIAL_RESIDUAL_SCALE
//                      << " -> " << (j*phi[i]*factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE) << std::endl;
//                  }


                  r_s.accumulate(lfsu_s,i,j*phi[i]*factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE);
                }
              }

            if (bctype==ConvectionDiffusionBoundaryConditions::Outflow)
              {
                // evaluate u
                RF u=0.0;
                for (size_type i=0; i<lfsu_s.size(); i++)
                  u += x_s(lfsu_s,i)*phi[i];

                // evaluate velocity field and outer unit normal
                typename T::Traits::RangeType b = param.b(*(ig.inside()),local);
                const Dune::FieldVector<DF,dim> n = ig.unitOuterNormal(it->position());

                // evaluate outflow boundary condition
                typename T::Traits::RangeFieldType o = param.o(ig.intersection(),it->position());
            
                // integrate o
                RF factor = it->weight()*ig.geometry().integrationElement(it->position());
                for (size_type i=0; i<lfsu_s.size(); i++)
                {
                  r_s.accumulate(lfsu_s,i,( (b*n)*u + o)*phi[i]*factor*volumeScale[i]*POTENTIAL_RESIDUAL_SCALE);
                }
              }
          }
      }


      /*
      // jacobian contribution from boundary
      template<typename IG_ORIG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_boundary (const IG_ORIG& ig_orig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              M& mat_s) const
      {
        typedef typename Acme2CylGeometrySwitch::IntersectionGeometrySwitch<IG_ORIG>::type IG;
        const IG& ig(ig_orig);

        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        
        // evaluate boundary condition type
        Dune::GeometryType gtface = ig.geometryInInside().type();
        Dune::FieldVector<DF,dim-1> facecenterlocal = Dune::ReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        ConvectionDiffusionBoundaryConditions::Type bctype;
        bctype = param.bctype(ig.intersection(),facecenterlocal);
 
        // skip rest if we are on Dirichlet boundary
        if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet) return;
        if (bctype==ConvectionDiffusionBoundaryConditions::Neumann) return;

        // Special case for dune-multidomain
        if (bctype==ConvectionDiffusionBoundaryConditions::None) return;

        // select quadrature rule
        const int intorder = intorderadd+2*lfsu_s.finiteElement().localBasis().order();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

            // evaluate shape functions (assume Galerkin method)
#if USE_CACHE
            const std::vector<RangeType>& phi = cache.evaluateFunction(local,lfsu_s.finiteElement().localBasis());
#else
            std::vector<RangeType> phi(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(local,phi);
#endif

            // evaluate velocity field and outer unit normal
            typename T::Traits::RangeType b = param.b(*(ig.inside()),local);
            const Dune::FieldVector<DF,dim> n = ig.unitOuterNormal(it->position());
        
            // integrate 
            RF factor = it->weight()*ig.geometry().integrationElement(it->position());
            for (size_type j=0; j<lfsu_s.size(); j++)
              for (size_type i=0; i<lfsu_s.size(); i++)
                mat_s.accumulate(lfsu_s,i,lfsu_s,j,(b*n)*phi[j]*phi[i]*factor);
          }
      }
      */



      //! set time in parameter class
      void setTime (double t)
      {
        param.setTime(t);
      }

      std::string getName()
      {
        return "ConvectionDiffusionFEM";
      }

    private:
      T& param;
      int intorderadd;
      const double POTENTIAL_RESIDUAL_SCALE;
      const int lfsuPotSize;
      const bool doVolumeScaling;
      mutable std::vector<double> volumeScale;
      mutable std::vector<Dune::FieldVector<Real,dim> > local_position;
      mutable std::vector<double> local_volume;
      mutable std::vector<Real> c;
      const bool useJacobianVolume;

      typedef typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;
      Dune::PDELab::LocalBasisCache<LocalBasisType> cache;
    };



    /** a local operator for residual-based error estimation
     *  
     * A call to residual() of a grid operator space will assemble
     * the quantity \eta_T^2 for each cell. Note that the squares
     * of the cell indicator \eta_T is stored. To compute the global
     * error estimate sum up all values and take the square root.
     *
     * Assumptions and limitations:
     * - Assumes that LFSU is P_1/Q_1 finite element space
     *   and LFSV is a P_0 finite element space (one value per cell).
     * - Convection term is ignored (but reaction term is included)
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T>
    class ConvectionDiffusionFEMResidualEstimator 
      : public Dune::PDELab::LocalOperatorDefaultFlags
    {
      enum { dim = T::Traits::GridViewType::dimension };
 
      typedef typename T::Traits::RangeFieldType Real;
      typedef typename ConvectionDiffusionBoundaryConditions::Type BCType;

    public:
      // pattern assembly flags
      enum { doPatternVolume = false };
      enum { doPatternSkeleton = false };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };

      //! constructor: pass parameter object
      ConvectionDiffusionFEMResidualEstimator (T& param_) 
        : param(param_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG_ORIG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG_ORIG& eg_orig, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        typedef typename Acme2CylGeometrySwitch::ElementGeometrySwitch<EG_ORIG>::type EG;
        const EG& eg(eg_orig);

        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int intorder = 2*lfsu.finiteElement().localBasis().order();

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        RF sum(0.0);
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // evaluate reaction term
            typename T::Traits::RangeFieldType c = param.c(eg.entity(),it->position());

            // evaluate right hand side parameter function
            typename T::Traits::RangeFieldType f = param.f(eg.entity(),it->position());

            // integrate f^2
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            sum += (f*f-c*c*u*u)*factor;
          }

        // accumulate cell indicator 
        DF h_T = diameter(eg.geometry());
        r.accumulate(lfsv,0,h_T*h_T*sum);
      }


      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG_ORIG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG_ORIG& ig_orig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           R& r_s, R& r_n) const
      {
        typedef typename Acme2CylGeometrySwitch::IntersectionGeometrySwitch<IG_ORIG>::type IG;
        const IG& ig(ig_orig);

        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>& 
          outside_local = Dune::ReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);
        typename T::Traits::PermTensorType A_s, A_n;
        A_s = param.A(*(ig.inside()),inside_local);
        A_n = param.A(*(ig.outside()),outside_local);

        // select quadrature rule
        const int intorder = 2*lfsu_s.finiteElement().localBasis().order();
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

        // tensor times normal
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        Dune::FieldVector<RF,dim> An_F_n;
        A_n.mv(n_F,An_F_n);

        // loop over quadrature points and integrate normal flux
        RF sum(0.0);
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // position of quadrature point in local coordinates of elements 
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());
            Dune::FieldVector<DF,dim> iplocal_n = ig.geometryInOutside().global(it->position());

            // evaluate gradient of basis functions
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
            std::vector<JacobianType> gradphi_n(lfsu_n.size());
            lfsu_n.finiteElement().localBasis().evaluateJacobian(iplocal_n,gradphi_n);

            // transform gradients of shape functions to real element
            jac = ig.inside()->geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
            jac = ig.outside()->geometry().jacobianInverseTransposed(iplocal_n);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_n(lfsu_n.size());
            for (size_type i=0; i<lfsu_n.size(); i++) jac.mv(gradphi_n[i][0],tgradphi_n[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu_s(0.0);
            for (size_type i=0; i<lfsu_s.size(); i++)
              gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);
            Dune::FieldVector<RF,dim> gradu_n(0.0);
            for (size_type i=0; i<lfsu_n.size(); i++)
              gradu_n.axpy(x_n(lfsu_n,i),tgradphi_n[i]);

            // integrate
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());
            RF jump = (An_F_s*gradu_s)-(An_F_n*gradu_n);
            sum += 0.25*jump*jump*factor;
          }

        // accumulate indicator
        // DF h_T = diameter(ig.geometry());
        DF h_T = std::max(diameter(ig.inside()->geometry()),diameter(ig.outside()->geometry()));
        r_s.accumulate(lfsv_s,0,h_T*sum);
        r_n.accumulate(lfsv_n,0,h_T*sum);
      }


      // boundary integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG_ORIG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG_ORIG& ig_orig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        typedef typename Acme2CylGeometrySwitch::IntersectionGeometrySwitch<IG_ORIG>::type IG;
        const IG& ig(ig_orig);

        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        typename T::Traits::PermTensorType A_s;
        A_s = param.A(*(ig.inside()),inside_local);
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);

        // select quadrature rule
        const int intorder = 2*lfsu_s.finiteElement().localBasis().order();
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

        // evaluate boundary condition
        const Dune::FieldVector<DF,dim-1> 
          face_local = Dune::ReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        BCType bctype = param.bctype(ig.intersection(),face_local);
        if (bctype != ConvectionDiffusionBoundaryConditions::Neumann)
          return;

        // loop over quadrature points and integrate normal flux
        RF sum(0.0);
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // position of quadrature point in local coordinates of elements 
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());

            // evaluate gradient of basis functions
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);

            // transform gradients of shape functions to real element
            jac = ig.inside()->geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu_s(0.0);
            for (size_type i=0; i<lfsu_s.size(); i++)
              gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);

            // evaluate flux boundary condition
            RF j = param.j(ig.intersection(),it->position());
                
            // integrate
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());
            RF jump = j+(An_F_s*gradu_s);
            sum += jump*jump*factor;
          }

        // accumulate indicator
        //DF h_T = diameter(ig.geometry());
        DF h_T = diameter(ig.inside()->geometry());
        r_s.accumulate(lfsv_s,0,h_T*sum);
      }

    private:
      T& param;  // two phase parameter class

      template<class GEO>
      typename GEO::ctype diameter (const GEO& geo) const
      {
        typedef typename GEO::ctype DF;
        DF hmax = -1.0E00;
        const int dim = GEO::coorddimension;
        for (int i=0; i<geo.corners(); i++)
          {
            Dune::FieldVector<DF,dim> xi = geo.corner(i);
            for (int j=i+1; j<geo.corners(); j++)
              {
                Dune::FieldVector<DF,dim> xj = geo.corner(j);
                xj -= xi;
                hmax = std::max(hmax,xj.two_norm());
              }
          }
        return hmax;
      }

    };
  }
}
#endif
