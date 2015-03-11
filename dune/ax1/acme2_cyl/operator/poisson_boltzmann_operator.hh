
#ifndef DUNE_AX1_POISSON_BOLTZMANN_OPERATOR_HH
#define DUNE_AX1_POISSON_BOLTZMANN_OPERATOR_HH

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>


template<typename PHYSICS, typename PARAMS, typename INITIAL_CON>
class PoissonBoltzmannLocalOperator :
public Dune::PDELab::NumericalJacobianApplyVolume<PoissonBoltzmannLocalOperator<PHYSICS,PARAMS,INITIAL_CON> >,
public Dune::PDELab::NumericalJacobianVolume<PoissonBoltzmannLocalOperator<PHYSICS,PARAMS,INITIAL_CON> >,
public Dune::PDELab::NumericalJacobianApplyBoundary<PoissonBoltzmannLocalOperator<PHYSICS,PARAMS,INITIAL_CON> >,
public Dune::PDELab::NumericalJacobianBoundary<PoissonBoltzmannLocalOperator<PHYSICS,PARAMS,INITIAL_CON> >,
public Dune::PDELab::FullVolumePattern,
public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
	// pattern assembly flags
	enum { doPatternVolume = true };
	
	// residual assembly flags
	enum { doAlphaVolume = true };
	
	static constexpr double EPS = 1e-8;

	PoissonBoltzmannLocalOperator (PHYSICS& physics_, PARAMS& param_, INITIAL_CON& gfInitialCon_,
	    int intorderadd_ = 0)
    : Dune::PDELab::NumericalJacobianApplyVolume<PoissonBoltzmannLocalOperator<PHYSICS,PARAMS,INITIAL_CON> >(EPS),
      Dune::PDELab::NumericalJacobianVolume<PoissonBoltzmannLocalOperator<PHYSICS,PARAMS,INITIAL_CON> >(EPS),
      Dune::PDELab::NumericalJacobianApplyBoundary<PoissonBoltzmannLocalOperator<PHYSICS,PARAMS,INITIAL_CON> >(EPS),
      Dune::PDELab::NumericalJacobianBoundary<PoissonBoltzmannLocalOperator<PHYSICS,PARAMS,INITIAL_CON> >(EPS),
      physics(physics_),
      param(param_),
      gfInitialCon(gfInitialCon_),
      intorderadd(intorderadd_),
      poissonConstant(physics_.getPoissonConstant())
	{}
	
	// volume integral depending on test and ansatz functions ########################################
  
	template<typename EG_ORIG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume (const EG_ORIG& eg_orig, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	{
	  typedef typename Acme2CylGeometrySwitch::ElementGeometrySwitch<EG_ORIG>::type EG;
    const EG& eg(eg_orig);

		// extract some types
		typedef typename LFSU::Traits::FiniteElementType::
		Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::
		Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::
		Traits::LocalBasisType::Traits::JacobianType JacobianType;
		typedef typename LFSU::Traits::FiniteElementType::
		Traits::LocalBasisType::Traits::RangeType RangeType;
		typedef typename LFSU::Traits::SizeType size_type;

		/*
		Dune::ios_base_all_saver nutt(std::cout);
    for (size_type i=0; i<lfsu.size(); i++)
    {
      debug_jochen << std::scientific << std::setprecision(16) << "x[i] = " << x(lfsu,i) << std::endl;
    }
    */
        
		// dimensions
		const int dim = EG::Geometry::dimension;
		const int dimw = EG::Geometry::dimensionworld;
		
		const int intorder = 2*lfsu.finiteElement().localBasis().order() + intorderadd;

		int elemIndex = physics.getElementIndex(eg.entity());

		// select quadrature rule
		Dune::GeometryType gt = eg.geometry().type();
		const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
		
		// loop over quadrature points
		for (typename Dune::QuadratureRule<DF,dim>::const_iterator 
			 it=rule.begin(); it!=rule.end(); ++it)
		{
			// evaluate basis functions on reference element
			std::vector<RangeType> phi(lfsu.size());
			lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);
			
			// compute u at integration point
			RF u=0.0;
			for (size_type i=0; i<lfsu.size(); i++)
				u += x(lfsu,i)*phi[i];
			

			debug_verb << std::scientific << std::setprecision(16) << "u  = " << u << std::endl;

			// evaluate gradient of basis functions on reference element
			std::vector<JacobianType> js(lfsu.size());
			lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);
			
			// transform gradients from reference element to real element
			const Dune::FieldMatrix<DF,dimw,dim> 
        jac = eg.geometry().jacobianInverseTransposed(it->position());
			std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
			for (size_type i=0; i<lfsu.size(); i++)
				jac.mv(js[i][0],gradphi[i]);
			
			// compute gradient of u
			Dune::FieldVector<RF,dim> gradu(0.0);
			for (size_type i=0; i<lfsu.size(); i++)
        gradu.axpy(x(lfsu,i),gradphi[i]);
			
      RF epsilon = physics.getPermittivity(elemIndex);

      // Calculate RHS of Poisson-Boltzmann equation
      RF f = 0.0;
      if(not physics.isMembrane(eg.entity()))
      {
        typename PHYSICS::FieldType poissonConstant = physics.getPoissonConstant();

        typename INITIAL_CON::Traits::RangeType conc(0.0);
        gfInitialCon.evaluate(eg.entity(), it->position(), conc);

        typename INITIAL_CON::Traits::RangeFieldType sum = 0.0;
        for (int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          int valence = physics.getValence(j);
          Dune::ios_base_all_saver huhn(std::cout);
          //debug_jochen << std::scientific << std::setprecision(6)
          //  << "u = " << u << " -> "
          //  << "conc[" << j << "] = " << (conc[j] * std::exp(-valence * u))
          //  << std::endl;
          sum += valence * conc[j] * std::exp(-valence * u);
        }
        //debug_jochen << "pot: " << u << " -> cd: " << sum << std::endl;
        f = - poissonConstant * sum;
        //debug_jochen << "poisson const = " << poissonConstant << std::endl;
        //debug_jochen << "f= " << f << std::endl;
        //debug_jochen << "cd = " << cd << ", exp-factor = " << (sum/cd) << ", sum = " << sum << std::endl;
      }
      
			RF a =  0.0; 
			
			// integrate grad u * grad phi_i + a*u*phi_i - f phi_i
			RF factor = it->weight()*eg.geometry().integrationElement(it->position());
			//debug_verb << "acc1" << std::endl;
			for (size_type i=0; i<lfsu.size(); i++)
			{
				//r.accumulate(lfsu,i,(epsilon*gradu*gradphi[i] + a*u*phi[i] - f*phi[i])*factor);
			  r.accumulate(lfsu,i,(epsilon*gradu*gradphi[i] + a*u*phi[i] + f*phi[i])*factor);
			}
		}
	}

  // boundary integral #############################################################################
  
  template<typename IG_ORIG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG_ORIG& ig_orig, const LFSU& lfsu_s, const X& x_s,
                       const LFSV& lfsv_s, R& r_s) const
  {
    typedef typename Acme2CylGeometrySwitch::IntersectionGeometrySwitch<IG_ORIG>::type IG;
    const IG& ig(ig_orig);

    // some types
    typedef typename LFSV::Traits::FiniteElementType::
    Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSV::Traits::FiniteElementType::
    Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSV::Traits::FiniteElementType::
    Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename LFSV::Traits::SizeType size_type;
    
    // dimensions
    const int dim = IG::dimension;
    
    const int intorder = 2*lfsu_s.finiteElement().localBasis().order();
    int elemIndex = physics.getElementIndex(*ig.inside());

    // select quadrature rule for face
    Dune::GeometryType gtface = ig.geometryInInside().type();
    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);
    
    // evaluate boundary condition type
    Dune::FieldVector<DF,dim-1> facecenterlocal = Dune::ReferenceElements<DF,dim-1>::general(gtface).position(0,0);
    Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type bctype;
    bctype = param.bctype(ig.intersection(),facecenterlocal);

    // skip rest if we are on Dirichlet boundary
    if (bctype==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet) return;

    // loop over quadrature points and integrate normal flux
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); 
         it!=rule.end(); ++it)
    {
      // position of quadrature point in local coordinates of element 
      Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());
      
      // evaluate basis functions at integration point
      std::vector<RangeType> phi(lfsv_s.size());
      lfsu_s.finiteElement().localBasis().evaluateFunction(local,phi);
      
      if (bctype==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann)
      {
        // evaluate flux boundary condition
        typename PARAMS::Traits::RangeFieldType j = param.j(ig.intersection(),it->position());

        // evaluate u (e.g. flux may depend on u)
        //RF u=0.0;
        //for (size_type i=0; i<lfsu_s.size(); i++)
        //  u += x_s(lfsu_s,i)*phi[i];

        // integrate j
        RF factor = it->weight()*ig.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsv_s.size(); i++)
          r_s.accumulate(lfsv_s,i,j*phi[i]*factor);
      }
    }
  }
  
  std::string getName()
  {
    return "PoissonBoltzmannOperator";
  }

private:
  PHYSICS& physics;
  PARAMS& param;
  INITIAL_CON& gfInitialCon;
  const int intorderadd;

  typename PHYSICS::FieldType poissonConstant;
};

#endif /* DUNE_AX1_POISSON_BOLTZMANN_OPERATOR_HH */
