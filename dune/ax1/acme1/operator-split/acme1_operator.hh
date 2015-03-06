#ifndef DUNE_AX1_ACME1_OPERATOR_HH
#define DUNE_AX1_ACME1_OPERATOR_HH

#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

#include <dune/ax1/common/constants.hh>

template<class BCT, class Physics, class DGF_POT_GRAD>
class NernstPlanckOperator: public Dune::PDELab::NumericalJacobianApplyVolume<
		NernstPlanckOperator<BCT, Physics, DGF_POT_GRAD> >,
		public Dune::PDELab::NumericalJacobianVolume<NernstPlanckOperator<BCT,
				Physics, DGF_POT_GRAD> >,
		public Dune::PDELab::NumericalJacobianApplyBoundary<NernstPlanckOperator<
				BCT, Physics, DGF_POT_GRAD> >,
		public Dune::PDELab::NumericalJacobianBoundary<NernstPlanckOperator<BCT,
				Physics, DGF_POT_GRAD> >,
		public Dune::PDELab::FullVolumePattern,
		public Dune::PDELab::LocalOperatorDefaultFlags,
		public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:

	// pattern assembly flags
	enum
	{
		doPatternVolume = true
	};

	// residual assembly flags
	enum
	{
		doAlphaVolume = true
	};
	enum
	{
		doAlphaBoundary = true
	};

	// epsilon for numerical differentiation (default 1e-7)
	static const double eps = 1e-7;

	// constructor stores parameters
	NernstPlanckOperator(const BCT& bct_, Physics& physics_,
			DGF_POT_GRAD& dgfGradPot_, unsigned int intorder_ = 2) :
				Dune::PDELab::NumericalJacobianApplyVolume<NernstPlanckOperator<BCT,
						Physics, DGF_POT_GRAD> >(eps),
				Dune::PDELab::NumericalJacobianVolume<NernstPlanckOperator<BCT,
						Physics, DGF_POT_GRAD> >(eps),
				Dune::PDELab::NumericalJacobianApplyBoundary<NernstPlanckOperator<BCT,
						Physics, DGF_POT_GRAD> >(eps),
				Dune::PDELab::NumericalJacobianBoundary<NernstPlanckOperator<BCT,
						Physics, DGF_POT_GRAD> >(eps), bct(bct_), physics(physics_),
				dgfGradPot(dgfGradPot_), intorder(intorder_)
	{
	}

	// volume integral depending on test and ansatz functions
	template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	DUNE_DEPRECATED void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x,
			const LFSV& lfsv, R& r) const
	{
		// local function space for the NUMBER_OF_SPECIES components (assume Galerkin scheme U=V)
		typedef typename LFSU::template Child<0>::Type LFSU_SINGLE_CON;

		// Numer of local power function space must equal number of ion species
		assert(NUMBER_OF_SPECIES == physics.numOfSpecies());

		// domain and range field type (assume both components have same RF)
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType
				DF;
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType
				RF;
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType
				JacobianType;
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
				RangeType;
		typedef typename LFSU::Traits::SizeType size_type;

		// dimensions
		const int dim = EG::Geometry::dimension;
		const int dimw = EG::Geometry::dimensionworld;

		// select quadrature rule
		Dune::GeometryType gt = eg.geometry().type();
		const Dune::QuadratureRule<DF, dim>& rule =
				Dune::QuadratureRules<DF, dim>::rule(gt, intorder);

		// get element index
		int elemIndex = physics.getElementIndex(eg.entity());

		// get group index
		int subdomainIndex = physics.getSubdomainIndex(elemIndex);

		// current source test
		RF current = 0.0;

		bool useLogScaling = physics.getParams().useLogScaling();

		// loop over quadrature points #################################################################
		for (typename Dune::QuadratureRule<DF, dim>::const_iterator it =
				rule.begin(); it != rule.end(); ++it)
		{

			const Dune::FieldMatrix<DF, dimw, dim> jac =
					eg.geometry().jacobianInverseTransposed(it->position());

			RF factor = it->weight() * eg.geometry().integrationElement(
					it->position());

			// Evaluate gradPot DGF
			typename DGF_POT_GRAD::Traits::RangeType gradPotTemp;
			dgfGradPot.evaluate(eg.entity(), it->position(), gradPotTemp);
			RF gradPot = gradPotTemp[0];
			//debug_verb << "entity #" << elemIndex << ", gradPot @ " << it->position() << " = " << gradPot << std::endl;

			// *************** Concentration part *************************
			std::vector < LFSU_SINGLE_CON > v_lfsuSingleCon;
			std::vector < RF > v_singleCon;
			std::vector<Dune::FieldVector<RF, dim> > v_gradSingleCon;
			std::vector < std::vector<RangeType> > v_phiSingleCon;
			std::vector < std::vector<Dune::FieldVector<RF, dim> >
					> v_gradPhiSingleCon;

			for (int j = 0; j < NUMBER_OF_SPECIES; ++j)
			{
				const LFSU_SINGLE_CON& lfsuSingleCon = lfsu.child(j);

				// evaluate basis functions on reference element
				std::vector < RangeType > phiSingleCon(lfsuSingleCon.size());
				lfsuSingleCon.finiteElement().localBasis().evaluateFunction(
						it->position(), phiSingleCon);

				// compute Con at integration point
				RF singleCon = 0.0;
				// localIndex() maps dof within leaf space to all dofs within given element
				for (size_type i = 0; i < lfsuSingleCon.size(); i++)
					singleCon += x(lfsuSingleCon, i) * phiSingleCon[i];

				// evaluate gradient of basis functions on reference element
				std::vector < JacobianType > jsSingleCon(lfsuSingleCon.size());
				lfsuSingleCon.finiteElement().localBasis().evaluateJacobian(
						it->position(), jsSingleCon);

				// transform gradients from reference element to real element
				std::vector<Dune::FieldVector<RF, dim> > gradPhiSingleCon(
						lfsuSingleCon.size());
				for (size_type i = 0; i < lfsuSingleCon.size(); i++)
					jac.mv(jsSingleCon[i][0], gradPhiSingleCon[i]);

				// compute gradient of Con
				Dune::FieldVector<RF, dim> gradSingleCon(0.0);
				for (size_type i = 0; i < lfsuSingleCon.size(); i++)
					gradSingleCon.axpy(x(lfsuSingleCon, i), gradPhiSingleCon[i]);

				v_lfsuSingleCon.push_back(lfsuSingleCon);
				v_singleCon.push_back(singleCon);
				v_gradSingleCon.push_back(gradSingleCon);
				v_phiSingleCon.push_back(phiSingleCon);
				v_gradPhiSingleCon.push_back(gradPhiSingleCon);

			}

			// integrate components ****************************************************************
			// eq. 0: Nernst-Planck

			// To get the diffusion coefficients of the membrane, we first need to solve
			// the coupled system of ODEs for all the channels on this element
			double dt = physics.getTimeStep();

			// TODO Hack, calculate real deltaPot from gradPot DGF!
			double deltaPot = 0.;
			physics.updateDiffCoeffs(subdomainIndex, elemIndex, dt, deltaPot);

			for (int j = 0; j < NUMBER_OF_SPECIES; ++j)
			{
				RF diffCoeff = physics.getDiffCoeff(subdomainIndex, j, elemIndex);
				RF valence = physics.getValence(j);

				for (size_type i = 0; i < v_lfsuSingleCon[j].size(); ++i)
				{
					//debug_verb << "j=" << ION_NAMES[j] << ", i=" << i << ", gradCon = " << v_gradSingleCon[j]
					//  << ", gradPot = " << gradPot << std::endl;

					if (not useLogScaling)
					{
						r.accumulate(
								v_lfsuSingleCon[j],
								i,
								diffCoeff
										* (v_gradSingleCon[j] + valence * v_singleCon[j] * gradPot)
										* v_gradPhiSingleCon[j][i] - current * v_phiSingleCon[j][i]
										* factor);
					} else
					{
						r.accumulate(
								v_lfsuSingleCon[j],
								i,
								diffCoeff * std::exp(v_singleCon[j])
										* (v_gradSingleCon[j] + valence * gradPot)
										* v_gradPhiSingleCon[j][i] - current * v_phiSingleCon[j][i]
										* factor);
					}
				}
			}
		}
	}

	// boundary integral #############################################################################

	template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_boundary(const IG& ig, const LFSU& lfsu, const X& x,
			const LFSV& lfsv, R& r) const
	{
		// local function space for the NUMBER_OF_SPECIES components (assume Galerkin scheme U=V)
		typedef typename LFSU::template Child<0>::Type LFSU_SINGLE_CON;
		typedef typename BCT::template Child<0>::Type BCT_SINGLE_CON;

		// some types
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType
				DF;
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType
				RF;
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType
				JacobianType;
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
				RangeType;
		typedef typename LFSU::Traits::SizeType size_type;

		// dimensions
		const int dim = IG::dimension;
		const int dimw = IG::Geometry::dimensionworld;

		// select quadrature rule for face
		Dune::GeometryType gtface = ig.geometryInInside().type();
		const Dune::QuadratureRule<DF, dim - 1>& rule = Dune::QuadratureRules<DF,
				dim - 1>::rule(gtface, intorder);

		// loop over quadrature points and integrate normal flux
		for (typename Dune::QuadratureRule<DF, dim - 1>::const_iterator it =
				rule.begin(); it != rule.end(); ++it)
		{
			// Concentrations ////////////////////////////////////////////////////////////////////////////
			for (int j = 0; j < NUMBER_OF_SPECIES; ++j)
			{
				const BCT_SINGLE_CON& bctSingleCon = bct.child(j);

				if (not bctSingleCon.isDirichlet(ig, it->position()))
				{
					if (bctSingleCon.nonZeroBoundaryFlux(ig, it->position(), j))
					{
						const LFSU_SINGLE_CON& lfsuSingleCon = lfsu.child(j);

						// position of quadrature point in local coordinates of element
						Dune::FieldVector<DF, dim> local = ig.geometryInInside().global(
								it->position());

						// evaluate basis functions on reference element
						std::vector < RangeType > phiSingleCon(lfsuSingleCon.size());
						lfsuSingleCon.finiteElement().localBasis().evaluateFunction(local,
								phiSingleCon);

						// evaluate flux boundary condition
						//Dune::FieldVector<RF,dim> globalpos = ig.geometry().global(it->position());
						RF J = 0.1;

						// integrate J
						RF factor = it->weight() * ig.geometry().integrationElement(
								it->position());
						for (size_type i = 0; i < lfsuSingleCon.size(); i++)
							r.accumulate(lfsuSingleCon, i, -J * phiSingleCon[i] * factor);
					}
				}
			}
		}
	}


private:
	const BCT& bct;
	Physics& physics;
	DGF_POT_GRAD& dgfGradPot;
	unsigned int intorder;
};

template<class BCT, class Physics, class DGF_CON>
class PoissonOperator: public Dune::PDELab::NumericalJacobianApplyVolume<
		PoissonOperator<BCT, Physics, DGF_CON> >,
		public Dune::PDELab::NumericalJacobianVolume<PoissonOperator<BCT, Physics,
				DGF_CON> >,
		public Dune::PDELab::NumericalJacobianApplyBoundary<PoissonOperator<BCT,
				Physics, DGF_CON> >,
		public Dune::PDELab::NumericalJacobianBoundary<PoissonOperator<BCT,
				Physics, DGF_CON> >,
		public Dune::PDELab::FullVolumePattern,
		public Dune::PDELab::LocalOperatorDefaultFlags,
		public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
	// pattern assembly flags
	enum
	{
		doPatternVolume = true
	};

	// residual assembly flags
	enum
	{
		doAlphaVolume = true
	};
	enum
	{
		doAlphaBoundary = true
	};

	static const double eps = 1e-5;

	// constructor stores parameters
	PoissonOperator(const BCT& bct_, Physics& physics_, DGF_CON& dgfCon_,
			unsigned int intorder_ = 2) :
				Dune::PDELab::NumericalJacobianApplyVolume<PoissonOperator<BCT,
						Physics, DGF_CON> >(eps),
				Dune::PDELab::NumericalJacobianVolume<PoissonOperator<BCT, Physics,
						DGF_CON> >(eps),
				Dune::PDELab::NumericalJacobianApplyBoundary<PoissonOperator<BCT,
						Physics, DGF_CON> >(eps),
				Dune::PDELab::NumericalJacobianBoundary<PoissonOperator<BCT, Physics,
						DGF_CON> >(eps), bct(bct_), physics(physics_), dgfCon(dgfCon_),
				intorder(intorder_)
	{
	}

	// volume integral depending on test and ansatz functions
	template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x,
			const LFSV& lfsv, R& r) const
	{
		// Numer of local power function space must equal number of ion species
		assert(NUMBER_OF_SPECIES == physics.numOfSpecies());

		// domain and range field type (assume both components have same RF)
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType
				DF;
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType
				RF;
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType
				JacobianType;
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
				RangeType;
		typedef typename LFSU::Traits::SizeType size_type;

		// dimensions
		const int dim = EG::Geometry::dimension;
		const int dimw = EG::Geometry::dimensionworld;

		// select quadrature rule
		Dune::GeometryType gt = eg.geometry().type();
		const Dune::QuadratureRule<DF, dim>& rule =
				Dune::QuadratureRules<DF, dim>::rule(gt, intorder);

		// get element index
		int elemIndex = physics.getElementIndex(eg.entity());

		// get group index
		int subdomainIndex = physics.getSubdomainIndex(elemIndex);

		// get permittivity of the element
		RF permittivity = physics.getPermittivity(subdomainIndex);

		//TODO remove this, deltaPot needs to be calculated in NernstPlanckOperator from gradPot DGF instead
		// get potential jump through membrane
		RF deltaPot = 0.0;
		if (physics.isMembrane(subdomainIndex)) // membrane
		{
			RF potMin = 0.0;
			RF potMax = 0.0;
			for (typename Dune::QuadratureRule<DF, dim>::const_iterator it =
					rule.begin(); it != rule.end(); ++it)
			{
				std::vector < RangeType > phiPot(lfsu.size());
				lfsu.finiteElement().localBasis().evaluateFunction(it->position(),
						phiPot);
				RF Pot = 0.0;
				for (size_type i = 0; i < lfsu.size(); i++)
					Pot += x(lfsu, i) * phiPot[i];
				if (Pot < potMin)
					potMin = Pot;
				if (Pot > potMax)
					potMax = Pot;
			}
			deltaPot = std::abs(potMax - potMin);
			//debug_verb << "deltaPot = " << deltaPot << std::endl;
		}

		// loop over quadrature points #################################################################
		for (typename Dune::QuadratureRule<DF, dim>::const_iterator it =
				rule.begin(); it != rule.end(); ++it)
		{
			const Dune::FieldMatrix<DF, dimw, dim> jac =
					eg.geometry().jacobianInverseTransposed(it->position());

			RF factor = it->weight() * eg.geometry().integrationElement(
					it->position());

			// Calculate charge density from ion concentrations
			RF chargeDensity = 0.0;
			typename DGF_CON::Traits::RangeType tempCon(0.0);

			dgfCon.evaluate(eg.entity(), it->position(), tempCon);

			for (int ion = 0; ion < NUMBER_OF_SPECIES; ++ion)
			{
				chargeDensity += (physics.getElectrolyte().getValence(ion) * tempCon[ion]);
			}
			//debug_verb << "chargeDensity for element #" << elemIndex << " @ "
			//    << it->position() << " = "<< chargeDensity << std::endl;


			// *************** Potential part *************************
			// evaluate basis functions on reference element
			std::vector < RangeType > phiPot(lfsu.size());
			lfsu.finiteElement().localBasis().evaluateFunction(it->position(), phiPot);

			// compute Pot at integration point
			RF Pot = 0.0;
			// localIndex() maps dof within leaf space to all dofs within given element
			for (size_type i = 0; i < lfsu.size(); i++)
				Pot += x(lfsu, i) * phiPot[i];

			// evaluate gradient of basis functions on reference element
			std::vector < JacobianType > jsPot(lfsu.size());
			lfsu.finiteElement().localBasis().evaluateJacobian(it->position(), jsPot);

			// transform gradients from reference element to real element
			std::vector<Dune::FieldVector<RF, dim> > gradPhiPot(lfsu.size());
			for (size_type i = 0; i < lfsu.size(); i++)
				jac.mv(jsPot[i][0], gradPhiPot[i]);

			// compute gradient of Pot
			Dune::FieldVector<RF, dim> gradPot(0.0);
			for (size_type i = 0; i < lfsu.size(); i++)
				gradPot.axpy(x(lfsu, i), gradPhiPot[i]);

			// integrate components ****************************************************************
			// eq. 1: Poisson
			RF poissonConstant = physics.getPoissonConstant();

			for (size_type i = 0; i < lfsu.size(); ++i)
			{
				r.accumulate(
						lfsu,
						i,
						(-gradPot * gradPhiPot[i] + chargeDensity / permittivity
								* poissonConstant * phiPot[i])
								* factor);
			}
		}
	}

	// boundary integral #############################################################################

	template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_boundary(const IG& ig, const LFSU& lfsu, const X& x,
			const LFSV& lfsv, R& r) const
	{
		// some types
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType
				DF;
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType
				RF;
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType
				JacobianType;
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
				RangeType;
		typedef typename LFSU::Traits::SizeType size_type;

		// dimensions
		const int dim = IG::dimension;
		const int dimw = IG::Geometry::dimensionworld;

		// select quadrature rule for face
		Dune::GeometryType gtface = ig.geometryInInside().type();
		const Dune::QuadratureRule<DF, dim - 1>& rule = Dune::QuadratureRules<DF,
				dim - 1>::rule(gtface, intorder);

		// loop over quadrature points and integrate normal flux
		for (typename Dune::QuadratureRule<DF, dim - 1>::const_iterator it =
				rule.begin(); it != rule.end(); ++it)
		{
			// Potential /////////////////////////////////////////////////////////////////////////////////
			if (not bct.isDirichlet(ig, it->position()))
			{
				// position of quadrature point in local coordinates of element
				Dune::FieldVector<DF, dim> local = ig.geometryInInside().global(
						it->position());

				// evaluate basis functions at integration point
				std::vector < RangeType > phiPot(lfsu.size());
				lfsu.finiteElement().localBasis().evaluateFunction(local, phiPot);

				// evaluate flux boundary condition
				//Dune::FieldVector<RF,dim> globalpos = ig.geometry().global(it->position());
				RF J = 0.0;

				// integrate J
				RF factor = it->weight() * ig.geometry().integrationElement(
						it->position());
				for (size_type i = 0; i < lfsu.size(); i++)
					r.accumulate(lfsu, i, J * phiPot[i] * factor);
			}
		}
	}

private:
	const BCT& bct;
	Physics& physics;
	DGF_CON& dgfCon;
	unsigned int intorder;
};

#endif /* DUNE_AX1_ACME1_OPERATOR_HH */
