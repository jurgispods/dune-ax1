#ifndef DUNE_AX1_TOOLS_HH
#define DUNE_AX1_TOOLS_HH

#include <valarray>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cxxabi.h>
#include <iostream>
#include <fstream>

#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/scsgmapper.hh>

#include <dune/ax1/common/constants.hh>

class Tools
{
  
public:
  
  //! \brief get positions and solution out of the discrete grid function
  template<typename DGF, typename T>
  static void getSolutionVector(const DGF& udgf, int numSamplingPoints,
                                std::valarray<T>& pos,
                                std::valarray<T>& sol)
  {
    typedef typename DGF::Traits::GridViewType GV;
    typedef typename GV::template Codim<0>::Entity Entity;
    typedef typename GV::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GV::template Codim<0>::
      template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator ElementLeafIterator;
    
    const GV& gv = udgf.getGridView();
    int solSize = gv.size(0) * (numSamplingPoints + 2);
    
    pos.resize(solSize);
    sol.resize(solSize);
    
    typename DGF::Traits::DomainType x(0.);
    typename DGF::Traits::RangeType y;
    typename DGF::Traits::DomainType x_l;
    typename DGF::Traits::DomainType x_r;
    
    int i=0;
    for (ElementLeafIterator eit=gv.template begin<0,Dune::PartitionIteratorType::Interior_Partition>();
        eit!=gv.template end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++eit) {
      
      typename DGF::Traits::ElementType& e = *eit;
      
      x = typename DGF::Traits::DomainType(0.);
      x_l = e.geometry().corner(0);
      x_r = e.geometry().corner(1);
      
      for(int j=0; j<numSamplingPoints+2; ++j)
      {
        udgf.evaluate(e,x,y);
        pos[i] = (T) (x_l + (x_r-x_l)*x);
        sol[i] = (T) y[0];
        ++i;
        x += (1./(numSamplingPoints+1));
      }
      //udgf.evaluate(e,x,y);
    }
    
    // Handle last grid point
    //pos[i] = (T) (x_l + (x_r-x_l)*x);
    //sol[i] = (T) y[0];
  }

  // TODO Add a second function with multiple node output when the function is not continuous!!
  //! \brief get positions and solution out of the discrete grid function; insert extra entries into
  //! the solution vectors at subdomain interfaces (ionic milieu / membrane) when using different subdomains
  template<typename PHYSICS, typename DGF, typename T>
  static void getMultiGroupSolutionVector(PHYSICS& physics, const DGF& udgf, int numSamplingPoints,
                                std::valarray<T>& pos,
                                std::valarray<T>& sol)
  {
    typedef typename DGF::Traits::GridViewType GV;
    typedef typename GV::template Codim<0>::Entity Entity;
    typedef typename GV::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GV::template Codim<0>::
      template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator ElementLeafIterator;

    const GV& gv = udgf.getGridView();

    assert(physics.numberSubDomainInterfaces(gv) <= 2);

    // This is the number of subdomain interfaces!
    const int offset = physics.numberSubDomainInterfaces(gv);

    int solSize = gv.size(0) * (numSamplingPoints + 2); //+ gv.size(1) + offset;
    pos.resize(solSize);
    sol.resize(solSize);

    typename DGF::Traits::DomainType x(0.);
    typename DGF::Traits::RangeType y;
    typename DGF::Traits::DomainType x_l;
    typename DGF::Traits::DomainType x_r;

    int lastGroupIndex = -1;
    int subdomainIndex = -1;

    int i=0;
    for (ElementLeafIterator eit=gv.template begin<0,Dune::PartitionIteratorType::Interior_Partition>();
        eit!=gv.template end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++eit)
    {
      typename DGF::Traits::ElementType& e = *eit;

      subdomainIndex = physics.getSubdomainIndex(e);

      x = typename DGF::Traits::DomainType(0.);
      x_l = e.geometry().corner(0);
      x_r = e.geometry().corner(1);

      /*
      // Insert additional value at left corner of element
      // in case of a boundary/membrane interface
      if(subdomainIndex != lastGroupIndex)
      {
        pos[i] = (T) (x_l + (x_r-x_l)*x);
        udgf.evaluate(e,x,y);
        sol[i] = (T) y[0];
        ++i;
      }
      // Start on first subsampling point
      x += (1./(numSamplingPoints+1));
      */

      for(int j=0; j<numSamplingPoints+2; ++j)
      {
        pos[i] = (T) (x_l + (x_r-x_l)*x);
        udgf.evaluate(e,x,y);
        sol[i] = (T) y[0];
        ++i;
        x += (1./(numSamplingPoints+1));
      }
      lastGroupIndex = subdomainIndex;
    }
  }

  template<typename PHYSICS, typename DGF, typename T>
  static void getMembraneSolutionVector(PHYSICS& physics, const DGF& udgf, int numSamplingPoints,
                                std::valarray<T>& pos, // Extract domain type from DGF
                                std::valarray<typename DGF::Traits::RangeType>& sol)  // Extract range type from DGF
  {
    typedef typename DGF::Traits::GridViewType GV;
    typedef typename GV::template Codim<0>::Entity Entity;
    typedef typename GV::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GV::template Codim<0>::
      template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator ElementLeafIterator;

    typedef typename DGF::Traits::DomainType DT;
    typedef typename DGF::Traits::RangeType RT;

    const GV& gv = udgf.getGridView();

    int solSize = gv.size(0); // evaluate at element centers only
    pos.resize(solSize);
    sol.resize(solSize);

    typename DGF::Traits::DomainType x(0.);
    typename DGF::Traits::RangeType y;

    int i=0;
    for (ElementLeafIterator eit=gv.template begin<0,Dune::PartitionIteratorType::Interior_Partition>();
        eit!=gv.template end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++eit)
    {
      typename DGF::Traits::ElementType& e = *eit;

      x = typename DGF::Traits::DomainType(0.5);

      pos[i] = (T) e.geometry().center();
      udgf.evaluate(e,x,y);
      sol[i] = y;

      debug_verb << "+++ x = " << x << ", y = " << y << std::endl;
      ++i;
    }
  }

  template<typename PHYSICS, typename DGF, typename T>
  static void getMembraneSolutionVector(PHYSICS& physics, const DGF& udgf, int numSamplingPoints,
                                std::valarray<T>& pos, // Extract domain type from DGF
                                std::valarray<T>& sol)  // Extract range type from DGF
  {
    std::valarray<typename DGF::Traits::RangeType> yNest;
    Tools::getMembraneSolutionVector(physics, udgf, numSamplingPoints, pos, yNest);

    Tools::flattenVector(yNest,sol);
  }


  //! \brief get positions and solution out of the discrete grid function for DG elements
  template<typename DGF, typename T>
  static void getSolutionVectorDG(const DGF& udgf, int numSamplingPoints,
                                	std::valarray<T>& pos,
                                	std::valarray<T>& sol)
  {
    typedef typename DGF::Traits::GridViewType GV;
    typedef typename GV::template Codim<0>::Entity Entity;
    typedef typename GV::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GV::template Codim<0>::
      template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator ElementLeafIterator;

    const GV& gv = udgf.getGridView();
    int solSize = gv.size(0) * ( numSamplingPoints + 2 );

    pos.resize(solSize);
    sol.resize(solSize);

    typename DGF::Traits::DomainType x(0.);
    typename DGF::Traits::RangeType y;
    typename DGF::Traits::DomainType x_l;
    typename DGF::Traits::DomainType x_r;

    int i=0;
    for (ElementLeafIterator eit=gv.template begin<0,Dune::PartitionIteratorType::Interior_Partition>();
        eit!=gv.template end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++eit) {

      typename DGF::Traits::ElementType& e = *eit;

      x = typename DGF::Traits::DomainType(0.);
      x_l = e.geometry().corner(0);
      x_r = e.geometry().corner(1);


      for(int j=0; j<numSamplingPoints+2; ++j)
      {
        udgf.evaluate(e,x,y);
        pos[i] = (T) (x_l + (x_r-x_l)*x);
        sol[i] = (T) y[0];
        ++i;
        x += (1./(numSamplingPoints+1));
      }
    }
  }

  //! \brief get positions and values of the solution of the discrete grid function for elements centers
  template<typename DGF, typename T>
  static void getSolutionVectorCenter(const DGF& udgf,
                                			std::valarray<T>& pos,
                                			std::valarray<T>& sol)
  {
    typedef typename DGF::Traits::GridViewType GV;
    typedef typename GV::template Codim<0>::Entity Entity;
    typedef typename GV::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GV::template Codim<0>::
      template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator ElementLeafIterator;

    const GV& gv = udgf.getGridView();
    int solSize = gv.size(0);

    pos.resize(solSize);
    sol.resize(solSize);

    typename DGF::Traits::DomainType x(0.);

    typename DGF::Traits::DomainType x_l;
    typename DGF::Traits::DomainType x_r;

    typename DGF::Traits::RangeType y;
    typename DGF::Traits::RangeType y_l;
    typename DGF::Traits::RangeType y_r;

    int i=0;
    for (ElementLeafIterator eit=gv.template begin<0,Dune::PartitionIteratorType::Interior_Partition>();
        eit!=gv.template end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++eit)
    {
      typename DGF::Traits::ElementType& e = *eit;

      x = typename DGF::Traits::DomainType(0.);

      x = e.geometry().center();

      udgf.evaluate(e,e.geometry().local(x),y);

      pos[i] = (T) x;
      sol[i] = (T) y;
      ++i;
    }
  }

  //! \brief compute order of convergence for given error(dof)-curve through central differences
  template<typename T>
  static void convergenceOrder(const std::valarray<T>& dof,
                               const std::valarray<T>& err,
                                     std::valarray<T>& mdof,
                                     std::valarray<T>& ord)
  {
    for (int i=0; i<mdof.size(); ++i)
    {
      mdof[i] = 0.5 * ( dof[i] + dof[i+1] ); // midpoint between two dof values
      
      ord[i] = - ( log( err[i+1] ) - log( err[i] ) ) / ( log( dof[i+1] ) - log( dof[i] ) );
    }
  }
  
  //! \brief compute L2 norm of a function
  template<typename T>
  static T normL2(const std::valarray<T>& x, const std::valarray<T>& y)
  {
    T bernd = 0.0;
    for (int i=0; i<x.size()-1; ++i)
      bernd += ( x[i+1] - x[i] ) * 0.5 * ( y[i+1] * y[i+1] + y[i] * y[i] );
    return sqrt( bernd );
  }
  
  //! \brief locally refine 1d grid coordinates
  template<typename T>
  static void localRefineVector( std::vector<T>& v, const int level, T& upperBound )
  {
    T length = std::abs(upperBound - v[0]);
    int nPoints = std::pow(2, level)-1;

    int i = 0;
    while(v[i] < upperBound)
    {
      T delX = (v[i+1]-v[i])/(nPoints+1);

      T currX = v[i];
      for (int j=0; j<nPoints; ++j)
      {
        currX = currX + delX;
        v.push_back( currX );
      }
     i++;
    }

    sort( v.begin(), v.end() );
  }

  //! \brief globally refine 1d grid coordinates
  template<typename T>
  static void globalRefineVector( std::vector<T>& v, const int level, bool refineMembrane = false)
  {
    T length = std::abs(v[v.size()-1] - v[0]);
    int nPoints = std::pow(2, level)-1;
    T delX = length/(nPoints+1);

    T currX = v[0];
    for (int j=0; j<nPoints; ++j)
    {
      currX = currX + delX;
      v.push_back( currX );
    }

    // Refine also inside membrane
    if(refineMembrane and v[0] > 0)
    {
      nPoints = std::floor(v[0]/delX);
      v.push_back(0.0);
      currX = 0.0;
      delX = v[0]/(nPoints+1);
      debug_verb << nPoints << " -- " << delX << std::endl;
      for (int j=0; j<nPoints; ++j)
      {
        currX = currX + delX;
        v.push_back( currX );
      }
    }

    sort( v.begin(), v.end() );
  }

  //! \brief globally refine 1d grid coordinates nonuniformly (finer near membrane)
  template<typename T>
  static void globalRefineVectorLogarithmic( std::vector<T>& v, const int level )
  {
    T length = std::abs(v[v.size()-1] - v[0]);
    T exp = std::log(1+length);
    int nPoints = std::pow(2, level)-1;
    T delExp = exp/(nPoints+1);

    T currExp = v[0];
    for (int j=0; j<nPoints; ++j)
    {
      T currExp = currExp + delExp;
      v.push_back( v[0] -1 + std::exp(currExp));
    }
    sort( v.begin(), v.end() );
  }
  
  //! \brief compute exact (for linear finite elements) potential values at point p
  template<typename T>
  static DUNE_DEPRECATED T potentialExact(const std::valarray<T>& x,
                          const std::valarray<T>& r,
                          const std::valarray<T>& permittivity,
                          const T& p)
  {
    T bernd=0.0;
    for (int i=0; i<x.size()-1; ++i)
    {
      T s=(r[i+1]-r[i])/(x[i+1]-x[i]);

      T zwischenBernd=-0.5/permittivity[i]*(
        (s*std::pow(x[i+1],3)/3.0+(r[i]-s*(p+x[i]))*std::pow(x[i+1],2)/2.0+(s*x[i]-r[i])*p*x[i+1])
       -(s*std::pow(x[i  ],3)/3.0+(r[i]-s*(p+x[i]))*std::pow(x[i  ],2)/2.0+(s*x[i]-r[i])*p*x[i  ]));

      //if ( p<=x[i]   ) bernd += zwischenBernd;
      if ( p>=x[i+1] ) zwischenBernd = -zwischenBernd;

      /*
      debug_verb << "slope[" << i << "] = " << s << ", ";
      debug_verb << "source[" << i << "] = " << r[i] << ", ";
      debug_verb << "|x - x'| = " << (p-x[i]) << ", ";
      debug_verb << "zwischenbernd[" << i << "] = " << zwischenBernd << std::endl;
      */

      bernd += zwischenBernd;
    }
    return bernd;
  }
  
  //! \brief Multipol expansion of a given charge distribution in 1D at point p
  template<typename T>
  static T multipoleExpansion1D(const std::valarray<T>& x,
                                const std::valarray<T>& source,
                                const std::valarray<T>& permittivity,
                                const T& p,
                                const unsigned int order)
  {
    T bernd=0.0;
    for (int j=0; j<order+1; ++j)
    {
      T moment=0.0;
      for (int i=0; i<x.size()-1; ++i)
      {
        T slope=(source[i+1]-source[i])/(x[i+1]-x[i]);
        
        moment+=(((source[i]-slope*x[i])/(1+j)+slope/(2+j)*x[i+1])*std::pow(x[i+1],1+j)-
                 ((source[i]-slope*x[i])/(1+j)+slope/(2+j)*x[i  ])*std::pow(x[i  ],1+j))
                /permittivity[i];
        
        //moment+=0.5*(std::pow(x[i+1],j)*source[i+1]+std::pow(x[i],j)*source[i])*(x[i+1]-x[i]);
      }
      bernd+=moment*std::pow(p,j)/std::abs(std::pow(p,2*j+1));
    }
    return bernd * 0.25 / con_pi;
  }
  
  //! \brief compute second derivative of a valarray
  template<typename T>
  static void secondDerivative(const std::valarray<T>& x,
                               const std::valarray<T>& f,
                                     std::valarray<T>& d2f_dx2)
  {
    d2f_dx2 = 0.0;
    
    for (int i=1; i<x.size()-1; ++i)
    d2f_dx2[i] = -(( f[i+1] - f[i] ) / ( x[i+1] - x[i] ) - ( f[i] - f[i-1] ) / ( x[i] - x[i-1] ))
    * 2.0 / ( x[i+1] - x[i-1] );
  }
  
  //! \brief  compute integral of f(x)
  template<typename T>
  static T integral(const std::valarray<T>& x, const std::valarray<T>& f)
  {
    T bernd = 0.0;
    for (int i=0; i<x.size()-1; ++i)
    {
      bernd += 0.5 * ( f[i+1] + f[i] ) * ( x[i+1] - x[i] );
    }
    return bernd;
  }
  

  /** \brief Converts a variable of type T to std::string using stringstreams.
   *
   *  \tparam    T Type of variable to be put into a string.
   *  \param[in] t Value to be converted.
   *
   *  \returns   std::string containing the value converted to a string.
   */
  template <class T>
  static std::string toString(const T& t)
  {
    std::stringstream myStringStream;
    myStringStream << t;
    return myStringStream.str();
  }


  /** \brief Converts a variable from std::string to type T using stringstreams.
   *
   *  \tparam    T Type of variable to be read from the string.
   *  \param[in] str String to be converted.
   *
   *  \returns   T containing the value extracted from string.
   */
  template <class T>
  static T fromString(const std::string& str)
  {
    std::stringstream myStringStream(str);
    T t;
    myStringStream >> t;
    return t;
  }


  /** \brief Converts a variable of any float compliant type T to std::string using
   *         stringstreams with prescribed precision.
   *
   *  \tparam    T    Type of variable to be put into a string.
   *  \param[in] t    Value to be converted.
   *  \param[in] prec Precision used for conversion.
   *
   *  \returns   std::string containing the value converted to a string.
   */
  template <class T>
  static std::string floatToString(const T& t, int prec)
  {
    // write formatted number to String
    std::stringstream myStringStream;
    myStringStream << std::fixed << std::setprecision(prec) << t;
    return myStringStream.str();
  }

  template<class U>
  static void expCoefficientVector(U& u)
  {
    for(std::size_t i = 0; i<u.N(); ++i)
    {
      u[i] = std::exp(u[i]);
    }
  }

  // explicit time step computations ###############################################

  //! \brief compute maximum timestep for explicit advection
  template<typename T>
  static T explicitTimeStepAdvection(const T & dx, const T & velocity)
  {
    return dx / std::abs(velocity);
  }

  //! \brief compute maximum timestep for explicit diffusion
  template<typename T>
  static T explicitTimeStepDiffusion(const T & dx, const T & diffusionCoefficient)
  {
    return 0.25 * dx * dx / diffusionCoefficient;
  }

  //! \brief compute explicit time step
	template<typename T>
	static T explicitTimeStep(const T & dx, const T & velocity, const T & diffusionCoefficient)
	{
		if (std::abs(velocity) < 1e-12)
		{
			return Tools::explicitTimeStepDiffusion(dx, diffusionCoefficient);
		}
		else
		{
			T dtAdve = Tools::explicitTimeStepAdvection(dx, velocity);
			T dtDiff = Tools::explicitTimeStepDiffusion(dx, diffusionCoefficient);
			return std::min(dtAdve,dtDiff);
		}
	}

	/**
   * do for every element
   *   v = - D * z * potGrad;
   *   dx = ... // cell width
   *   dtMin = explicitTimeStep(dx, v, D)
   * dt = min(dtMin) // minimum over all cells
   *
   * @param leckererDuftenderHackbratenMitSpeckKruste Nomen est omen!
   * @return
   */
  template<typename GV, typename PHYSICS, typename DGF_GRAD_POT>
  static double getTimeStep(GV& gv, PHYSICS& physics, DGF_GRAD_POT& leckererDuftenderHackbratenMitSpeckKruste)
  {
    typedef double Real;
    const int dim = GV::dimension;
    const int intorder = 2;

    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef typename Dune::SingleCodimSingleGeomTypeMapper<GV, 0> ElementMapper;
    typedef typename DGF_GRAD_POT::Traits::DomainType DT;
    typedef typename DGF_GRAD_POT::Traits::DomainFieldType DF;
    typedef typename DGF_GRAD_POT::Traits::RangeType RT;
    typedef typename DGF_GRAD_POT::Traits::RangeFieldType RF;

    ElementMapper elementMapper(gv);

    RF dt = 17.0e45; // aha!
    RF dtAdv = dt;
    RF dtDiff = dt;

    Real dtAdvLocal = dt;
    Real dtDiffLocal = dt;
    for(int j=0; j<NUMBER_OF_SPECIES; ++j)
    {
      int z = physics.getValence(j);

      for (ElementIterator eit=gv.template begin<0>(); eit!=gv.template end<0>(); ++eit)
      {
        int elemIndex = elementMapper.map(*eit);
        Real D = physics.getDiffCoeff(j, elemIndex);

        DT cellWidth(eit->geometry().corner(0));
        cellWidth -= eit->geometry().corner(1);

        // TODO This works in 1D only!
        // Think about how to determine the smallest cell width for more dimensions!
        Real dx = cellWidth.infinity_norm();

        RT potGrad(0.0);
        leckererDuftenderHackbratenMitSpeckKruste.evaluate(*eit, eit->geometry().center(), potGrad);

        Real v = -D * z * potGrad;
        //Real dtLocal = Tools::explicitTimeStep(dx, v, D);
        dtAdvLocal = Tools::explicitTimeStepAdvection(dx, v);
        dtDiffLocal = Tools::explicitTimeStepDiffusion(dx, D);


        //debug_verb << "D=" << D << ", potGrad=" << potGrad << ", dx=" << dx << std::endl;
        //debug_verb << "dt for element " << elemIndex << " = " << dt << std::endl;

        dtAdv = std::min(dtAdv, dtAdvLocal);
        dtDiff = std::min(dtDiff, dtDiffLocal);
      }
    }
    debug_verb << "dt [adv] = " << dtAdv << std::endl;
    debug_verb << "dt [diff] = " << dtDiff << std::endl;
    dt = std::min(dtAdv, dtDiff);
    return dt;
  }

  // demangle !!!!!!!!!!! yeah! ####################################################

  template<typename T>
  static const std::string getTypeName(const T& type)
  {
    //char* demangled = __cxxabiv1::__cxa_demangle(typeid(type).name(),0,0,0);
    char* demangled = abi::__cxa_demangle(typeid(type).name(),0,0,0);
    std::string typeName(demangled);
    delete demangled;
    return typeName;
  }

  template<typename T>
  static const std::string getTypeName(const T* type)
  {
    char* demangled = abi::__cxa_demangle(typeid(type).name(),0,0,0);
    std::string typeName(demangled);
    delete demangled;
    return typeName;
  }

  template<typename GFS, typename U, typename U_CHILD>
  static void compositeToChildCoefficientVector(GFS& gfs, U& uFull, U_CHILD& uChild, int child)
  {
    typedef typename U_CHILD::Backend CHILD_BE;
    typedef typename U::Backend BE;

    int n = uChild.flatsize();
    for(int i=0; i<n; ++i)
    {
      CHILD_BE::access(uChild, i) = BE::access(uFull, gfs.subMap(child, i));
    }
  }

  template<typename GFS, typename U, typename U_CHILD>
  static void childToCompositeCoefficientVector(GFS& gfs, U_CHILD& uChild, int child, U& uFull)
  {
    typedef typename U_CHILD::Backend CHILD_BE;
    typedef typename U::Backend BE;

    int n = uChild.flatsize();
    for(int i=0; i<n; ++i)
    {
      BE::access(uFull, gfs.subMap(child, i)) = CHILD_BE::access(uChild, i);
    }
  }

  // TODO This method assumes the diffusion coefficient to be isotropic, this should be revised!
  template<typename GV, typename PARAMS, typename DGF_GRAD_POT>
  static void pecletNumber(const GV& gv, PARAMS& params, DGF_GRAD_POT& dgfGradPot)
  {
    //typedef typename GF::Traits::GridViewType GV;
    typedef typename GV::template Codim<0>::Entity Entity;
    typedef typename GV::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GV::template Codim<0>::
      template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator ElementLeafIterator;

    //const GV& gv = initalChargeDensity.getGridView();
    int solSize = gv.size(0);

    typename PARAMS::Traits::DomainType xMax(0.0);
    double pecletMax = 0;
    int ionMax = 0;

    for (ElementLeafIterator eit=gv.template begin<0,Dune::PartitionIteratorType::Interior_Partition>();
         eit!=gv.template end<0,Dune::PartitionIteratorType::Interior_Partition>(); ++eit)
    {

      const Entity& e = *eit;
      typename PARAMS::Traits::DomainType x(0.5);
      typename PARAMS::Traits::RangeType P(0.0);
      // Field type for scalar diffusion coefficient
      typename PARAMS::Traits::RangeFieldType diffusionCoefficient(0.0);

      typename PARAMS::Traits::DomainType cellWidth = e.geometry().corner(e.geometry().corners()-1);
      cellWidth -= e.geometry().corner(0);

      //debug_jochen << "-- x = " << e.geometry().global(x) << std::endl;
      typename DGF_GRAD_POT::Traits::RangeType gradPot(0.0);
      dgfGradPot.evaluate(e, x, gradPot);

      for(int i=0; i<NUMBER_OF_SPECIES; ++i)
      {
        // convection velocity
        P = params[i]->b(e, x, gradPot);
        diffusionCoefficient = params[i]->A(e, x)[0][0]; // TODO Only valid in isotropic case!

        // This is a vector with a Peclet number component in each direction
        for(int j=0; j<P.size(); j++)
        {
          P[j] *= (cellWidth[j] / diffusionCoefficient);
        }

        // Take the maximum of the "Peclet number vector"
        double cellPecletNumber = P.infinity_norm();

        if(cellPecletNumber > pecletMax)
        {
          pecletMax = cellPecletNumber;
          xMax = e.geometry().center();
          ionMax = i;
        }

        //debug_verb << ION_NAMES[i] << " Peclet number : "
          //  << cellPecletNumber << std::endl;
      }
    }
    debug_info << " == MAX Peclet number on this processor [" << ION_NAMES[ionMax] << "] @ x = " << xMax << " : "
       << pecletMax << std::endl;

    pecletMax = gv.comm().max(pecletMax);
    debug_info << " == MAX Peclet number (global) : " << pecletMax << std::endl;

  }

  //! Integrate a GridFunction only over a certain subdomain
  /**
   * ATTENTION! When using cylinder coordinates, use the corresponding method from
   * Acme2CylinderGeometryTools instead!
   *
   * Stolen and modified from
   * \code
   *  #include <dune/pdelab/common/functionutilities.hh>
   * \endcode
   *
   *
   */
  template<typename GF, typename PHYSICS>
  static void integrateGridFunctionOverSubdomain(const GF& gf, PHYSICS& physics,
                             int subdomain,
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

    assert(! USE_CYLINDER_COORDINATES);

    //debug_verb << " --- GROUP = " << elementGroup << std::endl;

    sum = 0;
    Range val;
    const EIterator eend = gf.getGridView().template end<0,Dune::Interior_Partition>();
    for(EIterator eit = gf.getGridView().template begin<0,Dune::Interior_Partition>(); eit != eend; ++eit)
    {
      // Skip all elements on a different subdomain
      if(physics.getSubdomainIndex(*eit) != subdomain) continue;

      //debug_verb << "Visiting element #" << physics.getElementIndex(*eit)
      //    << ", subdomainIndex  = " << physics.getSubdomainIndex(*eit) << std::endl;

      const Geometry& geo = eit->geometry();

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

  //! \brief Calculcated reversal potential for each ion species using the Nernst equation
  template<typename T>
  static T calcReversalPotential(T temp, T z, T conRatio)
  {
    T revPotential = - std::log(conRatio) * (temp * con_k) / (z * con_e);
    return revPotential;
  }

  //! \brief Calculate transmembrane flux from potential jump and channel conductance
  template<typename PHYSICS, typename RangeFieldType>
  static void calcTransMembraneFlux(PHYSICS& physics, const RangeFieldType& conductance,
      const RangeFieldType& conCytosol, const RangeFieldType& conExtracellular, const RangeFieldType& potJump,
      const int ionSpecies, RangeFieldType& flux)
  {
    typename PHYSICS::FieldType temp = physics.getElectrolyte().getTemperature();
    typename PHYSICS::FieldType stdCon = physics.getElectrolyte().getStdCon();
    int valence = physics.getValence(ionSpecies);
    double d_memb = physics.getParams().dMemb();

    const RangeFieldType conUp = std::max(conCytosol, conExtracellular);
    const RangeFieldType conDown = std::min(conCytosol, conExtracellular);

    const RangeFieldType conJump = conCytosol - conExtracellular;
    debug_jochen << "  [n] = " << conJump << ", n_up = " << conUp <<  ", [pot] = " << potJump << std::endl;
    flux = 0.0;

    // Direction decides in which direction the diffusion flux flows
    // (depends on which side of the membrane conUp and conDown are located)
    double direction = 1.;
    if(conJump < 0) direction = -1.;

    double diffTerm = direction * (conUp * std::log(conUp/conDown));
    if(std::abs(conJump) < 1e-8) diffTerm = 0.0; // Check for equal (esp. zero!) concentrations on both sides

    // drift term already has the right orientation (by signs of valence and -potJump)
    double driftTerm = (valence * conUp * -potJump);

    // Raw flux (del n + z * n_up * del phi) from Nernst-Planck equation
    flux = diffTerm;   // positive conJump => negative flux
    flux += driftTerm; // positive potJump => negative flux

    // Multiply with diffusion coefficient 10*g*(k*T)/(e^2 * z^2 * stdCon)
    double D = (10 * con_k * temp) / (valence * valence * con_e * con_e * stdCon) * conductance;
    flux *= D;

    debug_jochen << "  D = " << D << std::endl;
    debug_jochen << "  diffusion term = " << diffTerm
      << ", drift term = " << driftTerm
      << " (" << std::abs(diffTerm)/std::abs(diffTerm + driftTerm) << ")" << std::endl;
  }

  //! \brief Calculate transmembrane flux from potential jump and channel conductance
    template<typename PHYSICS, typename RangeFieldType>
    static void calcTransMembraneFluxOLD(PHYSICS& physics, const RangeFieldType& conductance,
        const RangeFieldType& conCytosol, const RangeFieldType& conExtracellular, const RangeFieldType& potJump,
        const int ionSpecies, RangeFieldType& flux)
    {
      typename PHYSICS::FieldType temp = physics.getElectrolyte().getTemperature();
      typename PHYSICS::FieldType stdCon = physics.getElectrolyte().getStdCon();
      int valence = physics.getValence(ionSpecies);
      double d_memb = physics.getParams().dMemb();

      const RangeFieldType conUp = std::max(conCytosol, conExtracellular);
      const RangeFieldType conDown = std::min(conCytosol, conExtracellular);
      const RangeFieldType conJump = conCytosol - conExtracellular;

      // Direction decides in which direction the diffusion flux flows
      // (depends on which side of the membrane conUp and conDown are located)
      double direction = 1.;
      if(conJump < 0) direction = -1.;

      debug_jochen << "[n] = " << conJump << ", n_up = " << conUp <<  ", [pot] = " << potJump << std::endl;

      flux = 0.0;
      // Raw flux (del n + z * n_up * del phi) from Nernst-Planck equation
      double diffTerm = conJump; // positiv conJump => negative flux
      double driftTerm = (valence * conUp * -potJump); // positive potJump => negativ flux
      // Multiply with diffusion coefficient 10*g*(k*T)/(e^2 * z^2 * stdCon)
      double D = (10 * con_k * temp) / (valence * valence * con_e * con_e * stdCon) * conductance;

      flux += diffTerm;
      flux += driftTerm;
      flux *= D;

      debug_jochen << "D = " << D << std::endl;
      debug_jochen << "diffusion term = " << diffTerm
        << ", drift term = " << driftTerm
        << " (" << std::abs(diffTerm)/std::abs(diffTerm + driftTerm) << ")" << std::endl;
    }

  //! \brief Calculate Na concentration to be injected into the cell to
  //! correspond to a current injection of I_inj [pA] over a time interval
  //! of dt [ms]
  template<typename PHYSICS, typename Real>
  static Real calculateNaInjectionFlux(const PHYSICS& physics, const Real& I_inj, const Real& dt)
  {
    int valence = physics.getValence(Na);
    Real naInjection = 0.0;
    naInjection = 1e-15 * I_inj * dt / (con_mol * valence * con_e);
    return naInjection;
  }

  /**
   * Flattens a nested valarray;
   * Note: Only a depth of 1 is used to dig into the nested vector for now
   *
   * @param v A nested vector such as std::valarray<Dune::FieldVector<double,3> >
   * @return A flat std::valarray such as std::valarray<double>
   */
  template<typename T, typename V>
  static void flattenVector(std::valarray<V>& v, std::valarray<T>& result)
  {
    result.resize(v.size() * v[0].size());
    int i=0;
    for(int j=0; j<v.size(); j++)
    {
      for(int k=0; k<v[0].size(); k++)
      {
        result[i] = v[j][k];
        i++;
      }
    }
  }

  /**
   * Flattens a nested valarray;
   * Note: Only a depth of 1 is used to dig into the nested vector for now
   *
   * @param v A nested vector such as std::valarray<Dune::FieldVector<double,3> >
   * @return A flat std::vector such as std::vector<double>
   */
  template<typename T, typename V>
  static void flattenVector(std::valarray<V>& v, std::vector<T>& result)
  {
    result.resize(v.size() * v[0].size());
    int i=0;
    for(int j=0; j<v.size(); j++)
    {
      for(int k=0; k<v[0].size(); k++)
      {
        result[i] = v[j][k];
        i++;
      }
    }
  }

  template<typename T, typename V>
  static void flattenVector(std::vector<V>& v, std::vector<T>& result)
  {
    result.resize(v.size() * v[0].size());
    int i=0;
    for(int j=0; j<v.size(); j++)
    {
      for(int k=0; k<v[0].size(); k++)
      {
        result[i] = v[j][k];
        i++;
      }
    }
  }

  template<typename T, typename V>
  static void flattenVector(V& v, std::vector<T>& result)
  {
    result.resize(v.size() * v[0].size());
    int i=0;
    for(int j=0; j<v.size(); j++)
    {
      for(int k=0; k<v[0].size(); k++)
      {
        result[i] = v[j][k];
        i++;
      }
    }
  }

  template<typename T, typename V>
  static void extractComponentFromVector(V& v, std::vector<T>& result, const int component = 0)
  {
    result.resize(v.size());
    for(int j=0; j<v.size(); j++)
    {
      result[j] = v[j][component];
    }
  }


  template<typename T, int size>
  static bool lessThan(const Dune::FieldVector<T,size>& v1, const Dune::FieldVector<T,size>& v2)
  {
    bool result = true;
    for(int i=0; i<size; i++)
    {
     result = result && v1[i] < v2[i];
    }
    return result;
  }

  template<typename T, int size>
  static bool lessOrEqualThan (const Dune::FieldVector<T,size>& v1, const Dune::FieldVector<T,size>& v2,
      T eps = 0.0)
  {
    bool result = true;
    for(int i=0; i<size; i++)
    {
      result = result && v1[i]-eps <= v2[i];
    }

    //debug_jochen << "Comparing " << v1 << " <= " << v2 << ": " << result << std::endl;

    return result;
  }

  template<typename T, int size>
  static bool greaterThan(const Dune::FieldVector<T,size>& v1, const Dune::FieldVector<T,size>& v2)
  {
    bool result = true;
    for(int i=0; i<size; i++)
    {
      result = result && v1[i] > v2[i];
    }
    return result;
  }

  template<typename T, int size>
  static bool greaterOrEqualThan(const Dune::FieldVector<T,size>& v1, const Dune::FieldVector<T,size>& v2,
      T eps = 0.0)
  {
    bool result = true;
    for(int i=0; i<size; i++)
    {
     result = result && v1[i]+eps >= v2[i];
    }

    //debug_jochen << "Comparing " << v1 << " >= " << v2 << ": " << result << std::endl;

    return result;
  }

  static std::vector<std::string> &splitString(const std::string &s, char delim, std::vector<std::string> &elems) {
      std::stringstream ss(s);
      std::string item;
      while (std::getline(ss, item, delim)) {
          elems.push_back(item);
      }
      return elems;
  }

  static std::vector<std::string> splitString(const std::string &s, char delim) {
      std::vector<std::string> elems;
      Tools::splitString(s, delim, elems);
      return elems;
  }

  /**
   * This function searches a give file for lines beginning with a certain prefix. If the flag
   * findValue is set, it additionally searches for a value following after prefix. nSkip is the
   * number of lines to skip at the beginning of the file, and tol specifies the tolerance when
   * comparing tokens with value (which is, although, templated, implicitly assumed to be of numeric
   * type when using in conjunction with tol!)
   */
  template<typename T>
  static int findLine(std::ifstream& file, T& value, std::string prefix = std::string(""),
      bool findValue = false, int nSkip = 0, T tol = 1e-4)
  {
    // Make sure we are at the beginning of the file
    file.clear();
    file.seekg(0, std::ios::beg);

    std::string line("");
    int nLine = 0;
    T find = value;

    for(int i=0; i<nSkip; i++)
    {
      std::getline(file, line);
      nLine++;
    }

    while(file.good())
    {
      std::getline(file, line);
      nLine++;

      if(line.find(prefix) == std::string::npos)
        continue;

      std::stringstream line_str(line.substr(prefix.size()));
      // Try to extract value
      line_str >> value;

      // If requested type could not be extracted, try extraction as long until
      // either success or input sequence was completely parsed.
      while(!line_str.eof() && line_str.bad())
      {
        line_str >> value;
      }

      if(! findValue)
        return nLine;

//      debug_jochen << "find=" << find << ", value=" << value << " in line " << nLine << std::endl;
//      debug_jochen << "  diff=" << std::abs(find-value) << nLine << std::endl;

      if(findValue && std::abs(find-value)<tol)
        return nLine;

    }
    return -1;
  }

  template<typename T>
  static int findValue(std::ifstream& file, T& value, const int nCol = 1, std::string prefix = std::string(""),
      bool findValue = false, int nSkip = 0, T tol = 1e-4)
  {
    // Make sure we are at the beginning of the file
    file.clear();
    file.seekg(0, std::ios::beg);

    std::string line("");
    int nLine = 0;
    T find = value;

    for(int i=0; i<nSkip; i++)
    {
      std::getline(file, line);
      nLine++;
    }

    while(file.good())
    {
      std::getline(file, line);
      nLine++;

      if(line.find(prefix) == std::string::npos)
        continue;

      int nFound = 0;
      std::stringstream line_str(line.substr(prefix.size()));
      // Try to extract value
      line_str >> value;

      nFound += (! line_str.bad());

      // If requested type could not be extracted, try extraction as long until
      // either success or input sequence was completely parsed.
      while(!line_str.eof() && nFound < nCol)
      {
        line_str >> value;

        nFound += (! line_str.bad());
      }

      if(! findValue)
        return nLine;

      if(findValue && std::abs(find-value)<tol)
        return nLine;

    }
    return -1;
  }

  //! Custom DOF permutation; deprecated; use PHYSICS::setupDOFPermutation instead
  template <typename PHYSICS>
  static void setupDOFPermutation(const PHYSICS& physics, const int n_con, const int n_pot,
      std::vector<std::size_t>& permute, std::vector<std::size_t>& inv_permute)
  {
    int nVerticesCytosol = physics.nVerticesCytosol();
    int nVerticesMembrane = physics.nVerticesMembrane();
    int nVerticesExtracellular = physics.nVerticesExtracellular();

    int nVerticesThisProcess = physics.nVertices();

    debug_info << "Calculated #vertices on this process, cytosol: " << nVerticesCytosol
        << ", membrane (inner): " << nVerticesMembrane << ", extracellular: " << nVerticesExtracellular
        << ", TOTAL: " << nVerticesThisProcess << std::endl;
    assert(nVerticesCytosol + nVerticesMembrane + nVerticesExtracellular == nVerticesThisProcess);

    int nDOFsCytosol = nVerticesCytosol * (NUMBER_OF_SPECIES+1);
    int nDOFsMembrane = nVerticesMembrane;
    int nDOFsExtracellular = nVerticesExtracellular * (NUMBER_OF_SPECIES+1);
    int nDOFsThisProcess = nDOFsCytosol + nDOFsMembrane + nDOFsExtracellular;

    debug_info << "Calculated #DOFs on this process, cytosol: " << nDOFsCytosol
        << ", membrane (inner): " << nDOFsMembrane << ", extracellular: " << nDOFsExtracellular
        << ", TOTAL: " << nDOFsThisProcess << std::endl;

    debug_info << "Setting up DOF permutation vectors..." << std::endl;

    const int n_total = n_con + n_pot;
    debug_jochen << "n_con = " << n_con  << std::endl;
    debug_jochen << "n_pot = " << n_pot  << std::endl;
    debug_jochen << "n_total = " << n_total  << std::endl;

    assert(n_con + n_pot == n_total);

    permute.resize(n_total);
    inv_permute.resize(n_total);

    int conIndex = 0;
    int potIndex = n_con;
    int count = 0;
    //debug_jochen << "New index mapping: " << std::endl;

    for (int i = 0; i < n_total; i++)
    {
      if (count < NUMBER_OF_SPECIES && conIndex < n_con
          && (i<nDOFsCytosol || i>= nDOFsCytosol+nDOFsMembrane))
      {
        permute[i] = conIndex;
        inv_permute[conIndex] = i;
        conIndex++;
        count++;
      } else
      {
        permute[i] = potIndex;
        inv_permute[potIndex] = i;
        potIndex++;
        count = 0; // reset counter
      }
      //debug_jochen << i << " -> " << permute[i] << std::endl;
    }

    //debug_jochen << "conIndex = " << conIndex << std::endl;
    //debug_jochen << "potIndex = " << potIndex << std::endl;
  }

  template<typename AmgParams>
  static void printAmgParams(const AmgParams& amgParams)
  {
    debug_info << "AMG parameters:" << std::endl;
    debug_info << "-----------------------------------------------" << std::endl;
    debug_info << "debugLevel: " << amgParams.debugLevel() << std::endl;
    debug_info << "getNoPreSmoothSteps: " << amgParams.getNoPreSmoothSteps() << std::endl;
    debug_info << "getNoPostSmoothSteps: " << amgParams.getNoPostSmoothSteps() << std::endl;
    debug_info << "getGamma: " << amgParams.getGamma() << std::endl;
    debug_info << "getAdditive: " << amgParams.getAdditive() << std::endl;
    debug_info << "maxLevel: " << amgParams.maxLevel() << std::endl;
    debug_info << "coarsenTarget: " << amgParams.coarsenTarget() << std::endl;
    debug_info << "minCoarsenRate: " << amgParams.minCoarsenRate() << std::endl;
    debug_info << "accumulate: " << amgParams.accumulate() << std::endl;
    debug_info << "getProlongationDampingFactor: " << amgParams.getProlongationDampingFactor() << std::endl;
    debug_info << "maxDistance: " << amgParams.maxDistance() << std::endl;
    debug_info << "skipIsolated: " << amgParams.skipIsolated() << std::endl;
    debug_info << "minAggregateSize: " << amgParams.minAggregateSize() << std::endl;
    debug_info << "maxAggregateSize: " << amgParams.maxAggregateSize() << std::endl;
    debug_info << "maxConnectivity: " << amgParams.maxConnectivity() << std::endl;
    debug_info << "beta: " << amgParams.beta() << std::endl;
    debug_info << "alpha: " << amgParams.alpha() << std::endl;
    debug_info << "-----------------------------------------------" << std::endl;
  }

  static void replaceAll(std::string& str, const std::string& from, const std::string& to)
  {
    if(from.empty())
      return;
    std::size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos)
    {
      str.replace(start_pos, from.length(), to);
      start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}


};




#endif // DUNE_AX1_TOOLS_HH
