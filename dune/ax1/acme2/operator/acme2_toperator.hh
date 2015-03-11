#ifndef DUNE_AX1_ACME1MD_TOPERATOR_HH
#define DUNE_AX1_ACME1MD_TOPERATOR_HH

#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>


#include <dune/ax1/common/constants.hh>

template<typename PHYSICS, typename FiniteElementMap>
class NernstPlanckTimeLocalOperator: public Dune::PDELab::NumericalJacobianApplyVolume<
		NernstPlanckTimeLocalOperator<PHYSICS> >,
		public Dune::PDELab::NumericalJacobianVolume<NernstPlanckTimeLocalOperator<PHYSICS> >,
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

	// constructor remembers parameters
	NernstPlanckTimeLocalOperator(PHYSICS& physics_) :
		physics(physics_)
	{
	}

	// volume integral depending on test and ansatz functions
	template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x,
			const LFSV& lfsv, R& r) const
	{
		// select only concentration component
		//typedef typename LFSU::template Child<0>::Type PLFSU_CON;
		//const PLFSU_CON& plfsuCon = lfsu.template getChild<0>();

		typedef typename LFSU::template Child<0>::Type LFSU_SINGLE_CON;

		// Numer of local power function space must equal number of ion species
		//assert(NUMBER_OF_SPECIES == physics.numOfSpecies());

		// domain and range field type
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

		// Determine FE and quadrature order
		unsigned int maxFEOrder = 0;
		for(int i=0; i<NUMBER_OF_SPECIES; i++)
		{
		  maxFEOrder = std::max(maxFEOrder, lfsu.child(i).finiteElement().localBasis().order());
		}
		const unsigned int intorder = 2 * maxFEOrder;

		// select quadrature rule
		Dune::GeometryType gt = eg.geometry().type();
		const Dune::QuadratureRule<DF, dim>& rule =
				Dune::QuadratureRules<DF, dim>::rule(gt, intorder);

		// loop over quadrature points
		for (typename Dune::QuadratureRule<DF, dim>::const_iterator it = rule.begin(); it != rule.end(); ++it)
		{
			RF factor = it->weight() * eg.geometry().integrationElement(it->position());

			// *************** Concentration part ************************
			for (int j = 0; j < NUMBER_OF_SPECIES; ++j)
			{
				const LFSU_SINGLE_CON& lfsuCon = lfsu.child(j);

				// evaluate basis functions on reference element
				std::vector<RangeType> phiCon(lfsuCon.size());
				lfsuCon.finiteElement().localBasis().evaluateFunction(it->position(),phiCon);

				// compute con at integration point
				RF con = 0.0;
				for (size_type i=0; i<lfsuCon.size(); ++i)
					con += x(lfsuCon, i)*phiCon[i];

				// integration
				for (size_type i = 0; i < lfsuCon.size(); ++i)
				{
					r.accumulate(lfsuCon, i, con * phiCon[i] * factor);
				}
			}
		}
	}
private:
	PHYSICS& physics;
	typedef typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;
  Dune::PDELab::LocalBasisCache<LocalBasisType> cache;
};

#endif /* DUNE_AX1_ACME1MD_TOPERATOR_HH */
