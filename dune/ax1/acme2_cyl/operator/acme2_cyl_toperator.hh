#ifndef DUNE_AX1_ACME1MD_TOPERATOR_HH
#define DUNE_AX1_ACME1MD_TOPERATOR_HH

#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/ax1_lfs_tools.hh>

template<typename PHYSICS, typename FiniteElementMap>
class NernstPlanckTimeLocalOperator: public Dune::PDELab::NumericalJacobianApplyVolume<
		NernstPlanckTimeLocalOperator<PHYSICS,FiniteElementMap> >,
		public Dune::PDELab::NumericalJacobianVolume<NernstPlanckTimeLocalOperator<PHYSICS,FiniteElementMap> >,
		public Dune::PDELab::FullVolumePattern,
		public Dune::PDELab::LocalOperatorDefaultFlags,
		public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  typedef typename PHYSICS::FieldType Real;

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

	enum { dim = PHYSICS::GridView::dimension };

	// constructor remembers parameters
	NernstPlanckTimeLocalOperator(PHYSICS& physics_, int intorderadd_ = 0) :
		physics(physics_),
		intorderadd(intorderadd_),
		lfsuPotSize(4), // TODO
		doVolumeScaling(physics.getParams().doVolumeScaling()),
    volumeScale(lfsuPotSize, 1.0),
    local_position(lfsuPotSize),
    local_volume(lfsuPotSize),
    c(lfsuPotSize),
    useJacobianVolume(physics.getParams().general.get("useTimeOperatorJacobianVolume",false))
	{
	}

	// volume integral depending on test and ansatz functions
	template<typename EG_ORIG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume(const EG_ORIG& eg_orig, const LFSU& lfsu, const X& x,
			const LFSV& lfsv, R& r) const
	{
	  //debug_jochen << "NernstPlanckTimeLocalOperator::alpha_volume @ "<< eg_orig.geometry().center() << std::endl;
	  typedef typename Acme2CylGeometrySwitch::ElementGeometrySwitch<EG_ORIG>::type EG;
    const EG& eg(eg_orig);

		// select only concentration component
		//typedef typename LFSU::template Child<0>::Type PLFSU_CON;
		//const PLFSU_CON& plfsuCon = lfsu.template getChild<0>();

		typedef typename LFSU::template Child<0>::Type LFSU_SINGLE_CON;

		// Number of local power function space must equal number of ion species
		//assert(NUMBER_OF_SPECIES == physics.numOfSpecies());

		// domain and range field type
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType
				DF;
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType
				RF;
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType
				JacobianType;
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainType
        DomainType;
		typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
				RangeType;
		typedef typename LFSU::Traits::SizeType size_type;

		// dimensions
		const int dim = EG::Geometry::dimension;

		const LFSU_SINGLE_CON& lfsu0 = lfsu.child(0);

    // Assume equal FE order
    unsigned int feOrder = lfsu0.finiteElement().localBasis().order();
    const unsigned int intorder = 2 * feOrder + intorderadd;

		std::vector<double> volumeScale(lfsu0.size(), 1.0);
    // ============================================================================================
    // Determine Lagrange points for each basis (=test) function
    // => get node position the test functions belong to
    if(doVolumeScaling)
    {
      for (int k=0; k<dim; k++)
      {
        CoordinateEvaluation f(k);
        lfsu0.finiteElement().localInterpolation().interpolate(f,c);
        for (size_type i=0; i<lfsu0.size(); i++)
          local_position[i][k] = c[i];
      }
      for (size_type i=0; i<lfsu0.size(); i++)
      {
        local_volume[i] = physics.getNodeVolume(eg.geometry().global(local_position[i]));
        volumeScale[i] = 1./local_volume[i];
      }
    }
    // ============================================================================================


		// select quadrature rule
		Dune::GeometryType gt = eg.geometry().type();
		const Dune::QuadratureRule<DF, dim>& rule =
				Dune::QuadratureRules<DF, dim>::rule(gt, intorder);

		// loop over quadrature points
		for (typename Dune::QuadratureRule<DF, dim>::const_iterator it = rule.begin(); it != rule.end(); ++it)
		{
			RF factor = it->weight() * eg.geometry().integrationElement(it->position());

      // evaluate basis functions on reference element
#if USE_CACHE
      const std::vector<RangeType>& phi = cache.evaluateFunction(it->position(),lfsu0.finiteElement().localBasis());
#else
      std::vector<RangeType> phi(lfsu0.size());
      lfsu0.finiteElement().localBasis().evaluateFunction(it->position(),phi);
#endif

			// *************** Concentration part ************************
			for (int j = 0; j < NUMBER_OF_SPECIES; ++j)
			{
				const LFSU_SINGLE_CON& lfsuCon = lfsu.child(j);

				// compute con at integration point
				RF con = 0.0;
				for (size_type i=0; i<lfsuCon.size(); ++i)
					con += x(lfsuCon, i)*phi[i];

				// integration
				for (size_type i = 0; i < lfsuCon.size(); ++i)
				{
//				  debug_jochen << "[TLOP] con @" << lfsuCon.dofIndex(i)
//				      << (con * phi[i] * factor * volumeScale[i]) << std::endl;
					r.accumulate(lfsuCon, i, con * phi[i] * factor * volumeScale[i]);
				}
			}
		}
	}

#ifdef AX1_USE_JACOBIAN_METHODS_IN_OPERATOR
	template<typename EG_ORIG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume (const EG_ORIG& eg_orig, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                        M& mat) const
  {
    if(! useJacobianVolume)
    {
      Dune::PDELab::NumericalJacobianVolume<NernstPlanckTimeLocalOperator<PHYSICS,
        FiniteElementMap> >::jacobian_volume(eg_orig, lfsu, x, lfsv,mat);
      return;
    }

    //debug_jochen << "NernstPlanckTimeLocalOperator::alpha_volume @ "<< eg_orig.geometry().center() << std::endl;
    typedef typename Acme2CylGeometrySwitch::ElementGeometrySwitch<EG_ORIG>::type EG;
    const EG& eg(eg_orig);

    // select only concentration component
    //typedef typename LFSU::template Child<0>::Type PLFSU_CON;
    //const PLFSU_CON& plfsuCon = lfsu.template getChild<0>();

    typedef typename LFSU::template Child<0>::Type LFSU_SINGLE_CON;

    // Number of local power function space must equal number of ion species
    //assert(NUMBER_OF_SPECIES == physics.numOfSpecies());

    // domain and range field type
    typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType
        DF;
    typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType
        RF;
    typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType
        JacobianType;
    typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainType
        DomainType;
    typedef typename LFSU_SINGLE_CON::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
        RangeType;
    typedef typename LFSU::Traits::SizeType size_type;

    // dimensions
    const int dim = EG::Geometry::dimension;

    const LFSU_SINGLE_CON& lfsu0 = lfsu.child(0);

    // Assume equal FE order
    unsigned int feOrder = lfsu0.finiteElement().localBasis().order();
    const unsigned int intorder = 2 * feOrder + intorderadd;

    std::vector<double> volumeScale(lfsu0.size(), 1.0);
    // ============================================================================================
    // Determine Lagrange points for each basis (=test) function
    // => get node position the test functions belong to
    if(doVolumeScaling)
    {
      for (int k=0; k<dim; k++)
      {
        CoordinateEvaluation f(k);
        lfsu0.finiteElement().localInterpolation().interpolate(f,c);
        for (size_type i=0; i<lfsu0.size(); i++)
          local_position[i][k] = c[i];
      }
      for (size_type i=0; i<lfsu0.size(); i++)
      {
        local_volume[i] = physics.getNodeVolume(eg.geometry().global(local_position[i]));
        volumeScale[i] = 1./local_volume[i];
      }
    }
    // ============================================================================================


    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF, dim>& rule =
        Dune::QuadratureRules<DF, dim>::rule(gt, intorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF, dim>::const_iterator it = rule.begin(); it != rule.end(); ++it)
    {
      RF factor = it->weight() * eg.geometry().integrationElement(it->position());

      // evaluate basis functions on reference element
#if USE_CACHE
      const std::vector<RangeType>& phi = cache.evaluateFunction(it->position(),lfsu0.finiteElement().localBasis());
#else
      std::vector<RangeType> phi(lfsu0.size());
      lfsu0.finiteElement().localBasis().evaluateFunction(it->position(),phi);
#endif

      // *************** Concentration part ************************
      for (int k=0; k<NUMBER_OF_SPECIES; ++k)
      {
        // integration
        for (size_type j=0; j<lfsu.child(k).size(); j++)
        {
          for (size_type i=0; i<lfsu.child(k).size(); ++i)
          {
            mat.accumulate(lfsu.child(k),i,lfsu.child(k),j,( phi[j]*phi[i] )*factor*volumeScale[i]);
          }
        }
      }
    }
  }
#endif

private:
	PHYSICS& physics;
	const int intorderadd;
	typedef typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;
  Dune::PDELab::LocalBasisCache<LocalBasisType> cache;

  const int lfsuPotSize;

  const bool doVolumeScaling;

  mutable std::vector<double> volumeScale;
  mutable std::vector<Dune::FieldVector<double,dim> > local_position;
  mutable std::vector<double> local_volume;
  mutable std::vector<double> c;

  const bool useJacobianVolume;

};

#endif /* DUNE_AX1_ACME1MD_TOPERATOR_HH */
