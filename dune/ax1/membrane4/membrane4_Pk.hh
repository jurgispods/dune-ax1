
#ifndef MEMBRANE1_PK_HH
#define MEMBRANE1_PK_HH

#include <valarray>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/io/file/gnuplot.hh>
#include <dune/pdelab/finiteelementmap/pk1dbasis.hh>
#include <dune/pdelab/newton/newton.hh>

#include <dune/ax1/membrane4/membrane4_operator.hh>

template<class GV>
void membrane4_Pk (const GV& gv, Physics<double>& physics,
                   std::valarray<double>& x, std::valarray<double>& f)
{
  
	// <<<1>>> Choose domain and range field type
	typedef typename GV::Grid::ctype Coord;
	typedef double Real;
	
	// <<<2>>> Make grid function space
	typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM;
	FEM fem(2);
  
	//typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ConformingDirichletConstraints CONSTRAINTS; // constraints class
	typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
	typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CONSTRAINTS,VBE> GFS;
	GFS gfs(gv,fem);
	BCTypeParam bctype; // boundary condition type
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints( bctype, gfs, cc ); // assemble constraints
  std::cout << "constrained dofs=" << cc.size() 
    << " of " << gfs.globalSize() << std::endl;
	
	// <<<3>>> Make grid operator space
	typedef Membrane4LocalOperator<BCTypeParam> LOP;                // operator including boundary
	LOP lop( bctype, physics, 4); //#######
	typedef VBE::MatrixBackend MBE;
	typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
	GO go(gfs,cc,gfs,cc,lop);
  
	// <<<4>>> Make FE function extending Dirichlet boundary conditions
  typedef typename GO::Traits::Domain U;
  U u(gfs,0.0);
  typedef BCExtension<GV,Real> G;                               // boundary value + extension
  G g(gv, physics);
  Dune::PDELab::interpolate(g,gfs,u);                           // interpolate coefficient vector
  
	// <<<5>>> Select a linear solver backend
	//typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
	//LS ls(5000,true);
	typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
  LS ls;
  
	// <<<6>>> solve nonlinear problem
	typedef typename GO::Traits::Domain U;
	//typedef typename Dune::PDELab::BackendVectorSelector<GFS,int>::Type U;

	//U u(gfs,0.0); // initial value
	Dune::PDELab::Newton<GO,LS,U> newton(go,u,ls);                         // <= NEW
	newton.setReassembleThreshold(0.0);
	newton.setVerbosityLevel(0); // (2)
	newton.setReduction(1e-14);
	newton.setMinLinearReduction(1e-4);
	newton.setMaxIterations(25);
	newton.setLineSearchMaxIterations(10);
  try {
    newton.apply();
  }
  catch (Dune::Exception e) {
    std::cout << e << std::endl;
    std::cout << "blub" << std::endl;
  }
	
	// <<<7>>> solution vector
	typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
	DGF udgf(gfs,u);
  
  /*
	Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
	vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
	vtkwriter.write("membrane3_Pk",Dune::VTKOptions::binaryappended);       // <= NEW
  */
  
	Tools::getSolutionVector(udgf, 8, x, f);
  
}

#endif
