/*
 * acme2_cyl_geometrycheck.hh
 *
 *  Created on: Oct 25, 2012
 *      Author: jpods
 */

#ifndef DUNE_ACME2_CYL_GEOMETRYCHECK_HH
#define DUNE_ACME2_CYL_GEOMETRYCHECK_HH


#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/utility/hierarchicsearch.hh>

#include <dune/ax1/common/constants.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_parametertree.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_geometrywrapper.hh>

class Acme2CylGeometryTools
{
  public:

    Acme2CylGeometryTools(const Acme2CylParameters& params_)
    : params(params_)
    {}

    template<typename GV, typename Real>
    void checkCylinderVolume(const GV& gv,
                             Real& sum,
                             unsigned qorder = 1)
    {
      typedef typename GV::template Codim<0>::
        template Partition<Dune::Interior_Partition>::Iterator EIterator;
      typedef typename GV::template Codim<0>::Geometry Geometry;
      typedef typename Geometry::ctype DF;
      static const int dimD = GV::dimension;
      typedef Dune::QuadratureRule<DF,dimD> QR;
      typedef Dune::QuadratureRules<DF,dimD> QRs;
      typedef typename QR::const_iterator QIterator;

      typedef Dune::FieldVector<Real,1> Range;

      sum = 0;
      Range val;
      const EIterator eend = gv.template end<0,
        Dune::Interior_Partition>();
      for(EIterator eit = gv.template begin<0,
          Dune::Interior_Partition>(); eit != eend; ++eit)
      {
        const Geometry& geoOrig = eit->geometry();
        typename Acme2CylGeometrySwitch::GeometrySwitch<Geometry>::type geo(geoOrig);

        Dune::GeometryType gt = geo.type();
        const QR& rule = QRs::rule(gt,qorder);
        const QIterator qend = rule.end();

        Range entitySum(0.0);
        for (QIterator qit=rule.begin(); qit != qend; ++qit)
        {
          // Test for volume: Integrate value of 1 over elements
          val = 1.0;

          debug_jochen << "GF Value @ " << geo.global(qit->position()) << ": " << val << std::endl;

          // accumulate error
          val *= qit->weight() * geo.integrationElement(qit->position());

          //debug_jochen << "quadrature weight " << qit->weight() << std::endl;
          //debug_jochen << "integration element " << geo.integrationElement(qit->position()) << std::endl;
          //debug_jochen << "element volume " << geo.volume() << std::endl;

          entitySum += val;
          sum += val;
        }

        double r1 = std::abs(geo.corner(0)[1]);
        double r2 = std::abs(geo.corner(geo.corners()-1)[1]);
        double h = std::abs(geo.corner(geo.corners()-1)[0] - geo.corner(0)[0]);
        debug_jochen << "= corner(0): " << geo.corner(0) << std::endl;
        debug_jochen << "= corner(" << (geo.corners()-1) << "): " << geo.corner(geo.corners()-1) << std::endl;

        double real_volume = getRealVolume(0, (r2 > r1 ? r1 : r2), (r2 > r1 ? r2 : r1), h, false);

        debug_jochen << "== Entity volume: " << geo.volume() << std::endl;
        debug_jochen << "== Integrated volume: " << entitySum[0] << std::endl;
        debug_jochen << "== Real volume: " << real_volume << std::endl;

        // Compare integral over entity with real volume
        if(std::abs(entitySum[0]-real_volume) > 1e-6)
          assert(std::abs(entitySum[0]/real_volume - 1.0) < 1e-6);
        // Compare geometry volume with real volume
        if(std::abs(geo.volume()-real_volume) > 1e-6)
          assert(std::abs(geo.volume()/real_volume - 1.0) < 1e-6);
        debug_jochen << "------------------------------------" << std::endl;
      }
      // Check complete cylinder volume
      double complete_volume = getRealVolume(0);
      sum = gv.comm().sum(sum);
      debug_jochen << "=== Full cylinder volume: " << sum << std::endl;
      debug_jochen << "=== Real cylinder volume: " << complete_volume << std::endl;
      assert(std::abs(sum/complete_volume - 1.0) < 1e-6);
    }

    template<typename GV, typename Real>
    void checkCylinderSurface(const GV& gv,
                              Real& sum,
                              unsigned qorder = 1) {
      typedef typename GV::template Codim<0>::
        template Partition<Dune::Interior_Partition>::Iterator EIterator;
      //typedef typename GV::template Codim<0>::Geometry Geometry;
      //typedef typename GV::template Codim<0>::
      //  template Partition<Dune::Interior_Partition>::Entity::LeafIntersectionIterator IIterator;
      typedef typename GV::IntersectionIterator IIterator;
      typedef typename IIterator::Intersection Intersection;
      typedef typename IIterator::Intersection Intersection;
      typedef typename Intersection::Geometry Geometry;
      typedef typename Intersection::Geometry::ctype DF;
      static const int dimD = GV::dimension - 1;
      typedef Dune::QuadratureRule<DF,dimD> QR;
      typedef Dune::QuadratureRules<DF,dimD> QRs;
      typedef typename QR::const_iterator QIterator;

      typedef Dune::FieldVector<Real,1> Range;
      sum = 0;
      Range val;
      const EIterator eend = gv.template end<0,
        Dune::Interior_Partition>();
      for(EIterator eit = gv.template begin<0,
          Dune::Interior_Partition>(); eit != eend; ++eit)
      {

        const IIterator iend = gv.iend(*eit);
        for(IIterator iit = gv.ibegin(*eit); iit != iend; ++iit)
        {
          const Geometry& geoOrig = iit->geometry();
          typename Acme2CylGeometrySwitch::GeometrySwitch<Geometry>::type geo(geoOrig);

          Dune::PDELab::IntersectionGeometry<Intersection> ig(*iit, -1);

          Dune::GeometryType gt = geo.type();
          const QR& rule = QRs::rule(gt,qorder);
          const QIterator qend = rule.end();

          Range entitySum(0.0);
          for (QIterator qit=rule.begin(); qit != qend; ++qit)
          {
            // Test for volume: Integrate over value of 1
            val = 1.0;

            debug_jochen << "GF Value @ " << geo.global(qit->position()) << ": " << val << std::endl;

            // accumulate error
            val *= qit->weight() * geo.integrationElement(qit->position());
            //debug_jochen << "quadrature weight " << qit->weight() << std::endl;
            //debug_jochen << "integration element " << geo.integrationElement(qit->position()) << std::endl;
            //debug_jochen << "element volume " << geo.volume() << std::endl;

            entitySum += val;
            sum += val;
          }

          bool isOrientedAxially = geo.isOrientedAxially();
          double r1 = std::abs(geo.corner(0)[1]);
          double r2 = std::abs(geo.corner(geo.corners()-1)[1]);
          double h = std::abs(geo.corner(geo.corners()-1)[0] - geo.corner(0)[0]);
          debug_jochen << "= corner(0): " << geo.corner(0) << std::endl;
          debug_jochen << "= corner(" << (geo.corners()-1) << "): " << geo.corner(geo.corners()-1) << std::endl;

          double real_volume = getRealVolume(1, (r2 > r1 ? r1 : r2), (r2 > r1 ? r2 : r1), h, isOrientedAxially);

          debug_jochen << "== Entity volume: " << geo.volume() << std::endl;
          debug_jochen << "== Integrated volume: " << entitySum[0] << std::endl;
          debug_jochen << "== Real volume: " << real_volume << std::endl;

          // Compare integral over entity with real volume
          if(std::abs(entitySum[0]-real_volume) > 1e-6)
            assert(std::abs(entitySum[0]/real_volume - 1.0) < 1e-6);
          // Compare geometry volume with real volume
          if(std::abs(geo.volume()-real_volume) > 1e-6)
            assert(std::abs(geo.volume()/real_volume - 1.0) < 1e-6);

          debug_jochen << "------------------------------------" << std::endl;
        }
      }
      // Check complete cylinder surface makes no sense, as this integration routine also takes all
      // the 'skeleton surfaces' into account!
    }


    template<typename GF>
    static void integrateGridFunctionOverCylinder(const GF& gf,
                               typename GF::Traits::RangeType& sum,
                               unsigned qorder = 1) {
      typedef typename GF::Traits::GridViewType GV;
      typedef typename GV::template Codim<0>::
        template Partition<Dune::Interior_Partition>::Iterator EIterator;
      typedef typename GV::template Codim<0>::Geometry GeometryOrig;
      typedef typename GF::Traits::RangeType Range;
      typedef typename GF::Traits::DomainFieldType DF;
      static const int dimD = GF::Traits::dimDomain;
      typedef Dune::QuadratureRule<DF,dimD> QR;
      typedef Dune::QuadratureRules<DF,dimD> QRs;
      typedef typename QR::const_iterator QIterator;

      sum = 0;
      Range val;
      const EIterator eend = gf.getGridView().template end<0,
        Dune::Interior_Partition>();
      for(EIterator eit = gf.getGridView().template begin<0,
          Dune::Interior_Partition>(); eit != eend; ++eit)
      {
        const GeometryOrig& geoOrig = eit->geometry();
        // Use switch to choose original or cylinder Geometry type
        typedef typename Acme2CylGeometrySwitch::GeometrySwitch<GeometryOrig>::type Geometry;
        const Geometry& geo(geoOrig);

        Dune::GeometryType gt = geo.type();
        const QR& rule = QRs::rule(gt,qorder);
        const QIterator qend = rule.end();

        for (QIterator qit=rule.begin(); qit != qend; ++qit)
        {
          // evaluate the given grid functions at integration point
          gf.evaluate(*eit,qit->position(),val);

          // accumulate error
          val *= qit->weight() * geo.integrationElement(qit->position());
          sum += val;
        }
      }
    }

    template<typename GF, typename PHYSICS>
    static void integrateGridFunctionOverCylinderSubdomain(const GF& gf, PHYSICS& physics,
                               int subdomain,
                               typename GF::Traits::RangeType& sum,
                               unsigned qorder = 1)
    {
      std::vector<int> subDomains(1, subdomain);
      integrateGridFunctionOverCylinderSubdomain(gf, physics, subDomains, sum, qorder);
    }

    //! Integrate a GridFunction only over a certain subdomain
    /**
     * Stolen and modified from
     * \code
     *   #include <dune/pdelab/common/functionutilities.hh>
     * \endcode
     */
    template<typename GF, typename PHYSICS>
    static void integrateGridFunctionOverCylinderSubdomain(const GF& gf, PHYSICS& physics,
                               const std::vector<int>& subdomain,
                               typename GF::Traits::RangeType& sum,
                               unsigned qorder = 1) {
      typedef typename GF::Traits::GridViewType GV;
      typedef typename GV::template Codim<0>::
        template Partition<Dune::Interior_Partition>::Iterator EIterator;
      typedef typename GV::template Codim<0>::Geometry Geometry;
      typedef typename GF::Traits::RangeType Range;
      typedef typename GF::Traits::DomainFieldType DF;
      static const int dimD = GF::Traits::dimDomain;
      typedef Dune::QuadratureRule<DF,dimD> QR;
      typedef Dune::QuadratureRules<DF,dimD> QRs;
      typedef typename QR::const_iterator QIterator;

      //debug_verb << " --- GROUP = " << elementGroup << std::endl;

      sum = 0;
      Range val;
      const EIterator eend = gf.getGridView().template end<0,
        Dune::Interior_Partition>();
      for(EIterator eit = gf.getGridView().template begin<0,
            Dune::Interior_Partition>(); eit != eend; ++eit)
      {

        bool matchingSubdomain = false;
        for(int s=0; s<subdomain.size(); ++s)
        {
          matchingSubdomain = matchingSubdomain || physics.getSubdomainIndex(*eit) == subdomain[s];
        }

        // Skip all elements on different subdomains
        if(! matchingSubdomain) continue;

        //debug_verb << "Visiting element #" << physics.getElementIndex(*eit)
        //    << ", subdomainIndex  = " << physics.getSubdomainIndex(*eit) << std::endl;

        const Geometry& geoOrig = eit->geometry();

        // This is either the above Geometry type or Acme2CylGeometryWrapper<Geometry>, depending on
        // whether cylinder coordinates are used or not
        const typename Acme2CylGeometrySwitch::GeometrySwitch<Geometry>::type& geo(geoOrig);

        Dune::GeometryType gt = geo.type();
        const QR& rule = QRs::rule(gt,qorder);
        const QIterator qend = rule.end();

        for (QIterator qit=rule.begin(); qit != qend; ++qit)
        {
          // evaluate the given grid functions at integration point
          gf.evaluate(*eit,qit->position(),val);

          //debug_jochen << "oo " << geo.center() << " : val = " << val << std::endl;

          // accumulate error
          val *= qit->weight() * geo.integrationElement(qit->position());
          sum += val;
        }
      }
    }


    template<typename GF>
    static void integrateIntersectionGridFunctionOverCylinder(const GF& gf,
                               typename GF::Traits::RangeType& sum,
                               unsigned qorder = 1,
                               bool onlyBoundary = true) {
      typedef typename GF::Traits::GridViewType GV;
      typedef typename GV::template Codim<0>::
        template Partition<Dune::Interior_Partition>::Iterator EIterator;
      //typedef typename GV::template Codim<0>::Geometry Geometry;
      //typedef typename GV::template Codim<0>::
      //  template Partition<Dune::Interior_Partition>::Entity::LeafIntersectionIterator IIterator;
      typedef typename GV::IntersectionIterator IIterator;
      typedef typename IIterator::Intersection Intersection;
      typedef typename IIterator::Intersection Intersection;
      typedef typename Intersection::Geometry GeometryOrig;
      typedef typename GF::Traits::RangeType Range;
      typedef typename GF::Traits::DomainFieldType DF;
      static const int dimD = GF::Traits::dimDomain;
      typedef Dune::QuadratureRule<DF,dimD> QR;
      typedef Dune::QuadratureRules<DF,dimD> QRs;
      typedef typename QR::const_iterator QIterator;

      sum = 0;
      Range val;
      const EIterator eend = gf.getGridView().template end<0,Dune::Interior_Partition>();
      for(EIterator eit = gf.getGridView().template begin<0,Dune::Interior_Partition>(); eit != eend; ++eit)
      {
        const IIterator iend = gf.getGridView().iend(*eit);
        for(IIterator iit = gf.getGridView().ibegin(*eit); iit != iend; ++iit)
        {
          if(onlyBoundary && !iit->boundary()) continue;

          const GeometryOrig& geoOrig = iit->geometry();
          // Use switch to choose original or cylinder Geometry type
          typedef typename Acme2CylGeometrySwitch::GeometrySwitch<GeometryOrig>::type Geometry;
          const Geometry& geo(geoOrig);

          // Update: Wrapping in IntersectionGeometry is not necessary anymore; boundary gridfunctions takes intersections
          // as first argument to evaluate() now!
//          typedef Dune::PDELab::IntersectionGeometry<Intersection> IG_ORIG;
//          IG_ORIG igOrig(*iit, -1);
//          // Use switch to choose original or cylinder IntersectionGeometry type
//          typedef typename Acme2CylGeometrySwitch::IntersectionGeometrySwitch<IG_ORIG>::type IG;
//          const IG& ig(igOrig);

          Dune::GeometryType gt = geo.type();
          const QR& rule = QRs::rule(gt,qorder);
          const QIterator qend = rule.end();

          for (QIterator qit=rule.begin(); qit != qend; ++qit)
          {
            // evaluate the given grid functions at integration point
            gf.evaluate(*iit,qit->position(),val);

            //debug_jochen << "@ x = " << geo.global(qit->position()) << ", y = " << val << std::endl;

            // accumulate error
            val *= qit->weight() * geo.integrationElement(qit->position());
            sum += val;
          }
        }
      }
    }

    template<typename GF, typename PHYSICS>
    static void integrateIntersectionGridFunctionOverCylinderSubdomain(const GF& gf,
                                PHYSICS& physics,
                                int subdomain,
                                typename GF::Traits::RangeType& sum,
                                unsigned qorder = 1,
                                bool onlyBoundary = true)
    {
      std::vector<int> subDomains(1, subdomain);
      integrateIntersectionGridFunctionOverCylinderSubdomain(gf, physics, subDomains, sum, qorder, onlyBoundary);
    }

    template<typename GF, typename PHYSICS>
    static void integrateIntersectionGridFunctionOverCylinderSubdomain(const GF& gf,
                               PHYSICS& physics,
                               const std::vector<int>& subdomain,
                               typename GF::Traits::RangeType& sum,
                               unsigned qorder = 1,
                               bool onlyBoundary = true) {
      typedef typename GF::Traits::GridViewType GV;
      typedef typename GV::template Codim<0>::
        template Partition<Dune::Interior_Partition>::Iterator EIterator;
      //typedef typename GV::template Codim<0>::Geometry Geometry;
      //typedef typename GV::template Codim<0>::
      //  template Partition<Dune::Interior_Partition>::Entity::LeafIntersectionIterator IIterator;
      typedef typename GV::IntersectionIterator IIterator;
      typedef typename IIterator::Intersection Intersection;
      typedef typename IIterator::Intersection Intersection;
      typedef typename Intersection::Geometry GeometryOrig;
      typedef typename GF::Traits::RangeType Range;
      typedef typename GF::Traits::DomainFieldType DF;
      static const int dimD = GF::Traits::dimDomain;
      typedef Dune::QuadratureRule<DF,dimD> QR;
      typedef Dune::QuadratureRules<DF,dimD> QRs;
      typedef typename QR::const_iterator QIterator;

      sum = 0;
      Range val;
      const EIterator eend = gf.getGridView().template end<0,Dune::Interior_Partition>();
      for(EIterator eit = gf.getGridView().template begin<0,Dune::Interior_Partition>(); eit != eend; ++eit)
      {
        bool matchingSubdomain = false;
        for(int s=0; s<subdomain.size(); ++s)
        {
          matchingSubdomain = matchingSubdomain || physics.getSubdomainIndex(*eit) == subdomain[s];
        }

        // Skip all elements on different subdomains
        if(! matchingSubdomain) continue;

        const IIterator iend = gf.getGridView().iend(*eit);
        for(IIterator iit = gf.getGridView().ibegin(*eit); iit != iend; ++iit)
        {
          if(onlyBoundary && !iit->boundary()) continue;

          const GeometryOrig& geoOrig = iit->geometry();
          // Use switch to choose original or cylinder Geometry type
          typedef typename Acme2CylGeometrySwitch::GeometrySwitch<GeometryOrig>::type Geometry;
          const Geometry& geo(geoOrig);

          // Update: Wrapping in IntersectionGeometry is not necessary anymore; boundary gridfunctions takes intersections
          // as first argument to evaluate() now!
//          typedef Dune::PDELab::IntersectionGeometry<Intersection> IG_ORIG;
//          IG_ORIG igOrig(*iit, -1);
//          // Use switch to choose original or cylinder IntersectionGeometry type
//          typedef typename Acme2CylGeometrySwitch::IntersectionGeometrySwitch<IG_ORIG>::type IG;
//          const IG& ig(igOrig);

          Dune::GeometryType gt = geo.type();
          const QR& rule = QRs::rule(gt,qorder);
          const QIterator qend = rule.end();

          for (QIterator qit=rule.begin(); qit != qend; ++qit)
          {
            // evaluate the given grid functions at integration point
            gf.evaluate(*iit,qit->position(),val);

            //debug_jochen << "@ x = " << geo.global(qit->position()) << ", y = " << val << std::endl;

            // accumulate error
            val *= qit->weight() * geo.integrationElement(qit->position());
            sum += val;
          }
        }
      }
    }




//    template<int codim>
//    void checkEntities(const Grid& grid)
//    {
//      typedef typename Grid::LeafGridView GV;
//      GV gv = grid.leafGridView();
//
//      typedef typename GV::template Codim<codim>::Entity Entity;
//      typedef typename GV::template Codim<codim>::Iterator EntityIterator;
//      typedef typename GV::template Codim<codim>::EntityPointer EntityPointer;
//      typedef Dune::SingleCodimSingleGeomTypeMapper<GV, codim> EntityMapper;
//
//      EntityMapper elementMapper(gv);
//
//      double volume = 0.0;
//      double area = 0.0;
//
//      // Tag elements by type
//      for(EntityIterator ep = gv.template begin<codim>(); ep != gv.template end<codim>(); ++ep)
//      {
//        // Geometry stuff
//        typedef typename EntityIterator::Entity::Geometry Geometry;
//        Geometry geo = ep->geometry();
//        //debug_jochen << Tools::getTypeName(geo) << std::endl;
//
//        typedef Acme2CylGeometryWrapper<Geometry::mydimension,Geometry::coorddimension,Grid> WrappedGeometry;
//        WrappedGeometry wrappedGeo(geo);
//        //Dune::Geometry<Geometry::mydimension,Geometry::coorddimension,Grid,WrappedGeometry> DuneGeo()
//
//        debug_jochen << "2D Volume: " << geo.volume() << std::endl;
//        debug_jochen << "Cylinder Volume: " << wrappedGeo.volume() << std::endl;
//
//        volume += wrappedGeo.volume();
//
//        typename EntityIterator::Entity::Geometry::GlobalCoordinate center = ep->geometry().center();
//
//        typename EntityIterator::Entity::Geometry::GlobalCoordinate diagonal
//          = ep->geometry().corner(ep->geometry().corners()-1);
//        diagonal -= ep->geometry().corner(0);
//
//        double h_y = std::abs(diagonal[1]);
//
//        typedef typename Entity::LeafIntersectionIterator EntityIntersectionIterator;
//        for(EntityIntersectionIterator iit = ep->ileafbegin(); iit != ep->ileafend(); ++iit)
//        {
//          typedef typename EntityIntersectionIterator::Intersection::Geometry EdgeGeometry;
//          EdgeGeometry edgeGeo = iit->geometry();
//          //debug_jochen << Tools::getTypeName(edgeGeo) << std::endl;
//
//          //typedef Acme2CylGeometryWrapper<1,2,Grid>
//          //  WrappedEdgeGeometry;
//          typedef Acme2CylGeometryWrapper<EdgeGeometry::mydimension,EdgeGeometry::coorddimension,Grid>
//            WrappedEdgeGeometry;
//          WrappedEdgeGeometry wrappedEdgeGeo(edgeGeo);
//          //Dune::Geometry<Geometry::mydimension,Geometry::coorddimension,Grid,WrappedGeometry> DuneGeo()
//
//
//          //debug_jochen << "wrapped geo mydim: " << EdgeGeometry::mydimension << std::endl;
//          //debug_jochen << "wrapped geo coorddim: " << EdgeGeometry::coorddimension << std::endl;
//
//          debug_jochen << "2D area: " << edgeGeo.volume() << std::endl;
//          debug_jochen << "Cylinder area: " << wrappedEdgeGeo.volume() << std::endl;
//
//          area += wrappedEdgeGeo.volume();
//        }
//      }
//
//      debug_jochen << "Total volume: " << volume << std::endl;
//      debug_jochen << "Real volume: " << getRealVolume(0) << std::endl;
//
//      debug_jochen << "Total area: " << area << std::endl;
//      debug_jochen << "Real area: " << getRealVolume(1) << std::endl;
//
//    }

  private:

    double getRealVolume(const int codim, double r1 = -1, double r2 = -1, double h = -1,
        bool isOrientedAxially = false)
    {

      double real_volume = 0.0;

      if(r1 < 0) r1 = 0;
      if(r2 < 0) r2 = params.yMax() - params.yMin();
      if(h < 0) h = params.xMax() - params.xMin();


      debug_jochen << "-Calculating real volume for this entity" << std::endl;
      debug_jochen << "-axial" << isOrientedAxially << std::endl;
      debug_jochen << "-r1 = " << r1 << std::endl;
      debug_jochen << "-r2 = " << r2 << std::endl;
      debug_jochen << "-h = " << h << std::endl;

      // Calculate volume
      if(codim == 0)
      {
        // Volume V = outer volume - inner volume
        real_volume = con_pi * r2 * r2 * h - con_pi * r1 * r1 * h;
      }

      // Calculate surface
      if(codim == 1)
      {
        if(isOrientedAxially)
        {
          assert(r1 == r2);
          // Side surface (A = 2*pi*r*h)
          real_volume = 2 * con_pi * r2 * h;
        } else {
          assert(r1 <= r2);
          // Base area (A = pi r2^2 - r1^2)
          real_volume = con_pi * r2 * r2 - con_pi * r1 * r1;
        }
      }

      return real_volume;

    }

  private:
    const Acme2CylParameters& params;
};



//! Simple as that: Return the volume of an element
template<typename GV, typename RF>
class Acme2CylVolumeGridFunction
  : public Dune::PDELab::GridFunctionBase<
          Dune::PDELab::GridFunctionTraits<GV,
                                           RF,
                                           1, // This should be enough for every channel (gbar, g and up to 3 gating particles)
                                           Dune::FieldVector<RF,1> >,
                                           Acme2CylVolumeGridFunction<GV,RF> >
{
  public:
    typedef Dune::PDELab::GridFunctionTraits<GV,   // grid view type
                                           RF, // image field type (double)
                                           1,                                        // number of components of image (k)
                                           Dune::FieldVector<RF,1> // image type (Dune::FieldVector<double, k>)
                                           > Traits;

    //! constructor
    Acme2CylVolumeGridFunction (const GV& gv_)
      : gridView(gv_)
    {
    }

    //! \copydoc GridFunctionBase::evaluate()
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                          typename Traits::RangeType& y) const
    {
      typename Acme2CylGeometrySwitch::GeometrySwitch<typename Traits::ElementType::Geometry>::type geo(e.geometry());
      y = geo.volume();
    }

    inline const typename Traits::GridViewType& getGridView () const
    {
      return gridView;
    }

  private:
    const typename Traits::GridViewType& gridView;
};

//! Return the plain 2D area of an element
template<typename GV, typename RF>
class Acme2CylPlain2DAreaGridFunction
  : public Dune::PDELab::GridFunctionBase<
          Dune::PDELab::GridFunctionTraits<GV,
                                           RF,
                                           1, // This should be enough for every channel (gbar, g and up to 3 gating particles)
                                           Dune::FieldVector<RF,1> >,
                                           Acme2CylPlain2DAreaGridFunction<GV,RF> >
{
  public:
    typedef Dune::PDELab::GridFunctionTraits<GV,   // grid view type
                                           RF, // image field type (double)
                                           1,                                        // number of components of image (k)
                                           Dune::FieldVector<RF,1> // image type (Dune::FieldVector<double, k>)
                                           > Traits;

    //! constructor
    Acme2CylPlain2DAreaGridFunction (const GV& gv_)
      : gridView(gv_)
    {
    }

    //! \copydoc GridFunctionBase::evaluate()
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                          typename Traits::RangeType& y) const
    {
      // Do not wrap geometry => get original area of 2D element, not cylinder volume
      y = e.geometry().volume();
    }

    inline const typename Traits::GridViewType& getGridView () const
    {
      return gridView;
    }

  private:
    const typename Traits::GridViewType& gridView;
};

//! Return the surface area of an element
template<typename GV, typename RF>
class Acme2CylSurfaceGridFunction
  : public Dune::PDELab::GridFunctionBase<
          Dune::PDELab::GridFunctionTraits<GV,
                                           RF,
                                           1, // This should be enough for every channel (gbar, g and up to 3 gating particles)
                                           Dune::FieldVector<RF,1> >,
                                           Acme2CylSurfaceGridFunction<GV,RF> >
{
  public:
    typedef Dune::PDELab::GridFunctionTraits<GV,   // grid view type
                                           RF, // image field type (double)
                                           1,                                        // number of components of image (k)
                                           Dune::FieldVector<RF,1> // image type (Dune::FieldVector<double, k>)
                                           > Traits;

    //! constructor
    Acme2CylSurfaceGridFunction (const GV& gv_)
      : gridView(gv_)
    {
    }

    //! \copydoc GridFunctionBase::evaluate()
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                          typename Traits::RangeType& y) const
    {
      y = 0.0;

      for(typename Traits::ElementType::LeafIntersectionIterator iit = gridView.ibegin(e);
          iit != gridView.iend(e); ++iit)
      {
        typedef typename Traits::ElementType::LeafIntersectionIterator::Intersection::Geometry IGEO_ORIG;
        typename Acme2CylGeometrySwitch::GeometrySwitch<IGEO_ORIG>::type geo(iit->geometry());

        y += geo.volume();
      }
    }

    inline const typename Traits::GridViewType& getGridView () const
    {
      return gridView;
    }

  private:
    const typename Traits::GridViewType& gridView;
};


//! Return the surface area of an element
template<typename GV, typename RF>
class Acme2CylSideSurfaceGridFunction
  : public Dune::PDELab::GridFunctionBase<
          Dune::PDELab::GridFunctionTraits<GV,
                                           RF,
                                           1, // This should be enough for every channel (gbar, g and up to 3 gating particles)
                                           Dune::FieldVector<RF,1> >,
                                           Acme2CylSideSurfaceGridFunction<GV,RF> >
{
  public:
    typedef Dune::PDELab::GridFunctionTraits<GV,   // grid view type
                                           RF, // image field type (double)
                                           1,                                        // number of components of image (k)
                                           Dune::FieldVector<RF,1> // image type (Dune::FieldVector<double, k>)
                                           > Traits;

    //! constructor
    Acme2CylSideSurfaceGridFunction (const GV& gv_)
      : gridView(gv_)
    {
    }

    //! \copydoc GridFunctionBase::evaluate()
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                          typename Traits::RangeType& y) const
    {
      y = 0.0;

      for(typename Traits::ElementType::LeafIntersectionIterator iit = gridView.ibegin(e);
          iit != gridView.iend(e); ++iit)
      {
        typedef typename Traits::ElementType::LeafIntersectionIterator::Intersection::Geometry IGEO_ORIG;

        // No use for Acme2CylGeometrySwitch here:
        // We need the cylinder geometry wrapper here in any case, as the isOrientedAxially() is needed later
        typedef Acme2CylGeometry<IGEO_ORIG> GEO;
        GEO geo(iit->geometry());

        // Only add surface area of axially oriented intersections
        if(geo.isOrientedAxially())
          y += geo.volume();
      }
    }

    inline const typename Traits::GridViewType& getGridView () const
    {
      return gridView;
    }

  private:
    const typename Traits::GridViewType& gridView;

};


//! Return the surface area of an element
template<typename GV, typename RF>
class Acme2CylBaseSurfaceGridFunction
  : public Dune::PDELab::GridFunctionBase<
          Dune::PDELab::GridFunctionTraits<GV,
                                           RF,
                                           1, // This should be enough for every channel (gbar, g and up to 3 gating particles)
                                           Dune::FieldVector<RF,1> >,
                                           Acme2CylBaseSurfaceGridFunction<GV,RF> >
{
  public:
    typedef Dune::PDELab::GridFunctionTraits<GV,   // grid view type
                                           RF, // image field type (double)
                                           1,                                        // number of components of image (k)
                                           Dune::FieldVector<RF,1> // image type (Dune::FieldVector<double, k>)
                                           > Traits;

    //! constructor
    Acme2CylBaseSurfaceGridFunction (const GV& gv_)
      : gridView(gv_)
    {
    }

    //! \copydoc GridFunctionBase::evaluate()
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                          typename Traits::RangeType& y) const
    {
      y = 0.0;

      for(typename Traits::ElementType::LeafIntersectionIterator iit = gridView.ibegin(e);
          iit != gridView.iend(e); ++iit)
      {
        typedef typename Traits::ElementType::LeafIntersectionIterator::Intersection::Geometry IGEO_ORIG;

        // No use for Acme2CylGeometrySwitch here:
        // We need the cylinder geometry wrapper here in any case, as the isOrientedAxially() is needed later
        typedef Acme2CylGeometry<IGEO_ORIG> GEO;
        GEO geo(iit->geometry());

        // Only add surface area of radially oriented intersections
        if(geo.isOrientedAxially())
          y += geo.volume();
      }
    }

    inline const typename Traits::GridViewType& getGridView () const
    {
      return gridView;
    }

  private:
    const typename Traits::GridViewType& gridView;

};




#endif /* DUNE_ACME2_CYL_GEOMETRYCHECK_HH */
