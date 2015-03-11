
#ifndef MEMBRANE1_PK_HH
#define MEMBRANE1_PK_HH

#include <dune/common/exceptions.hh>
#include <dune/grid/io/file/gnuplot.hh>
#include <dune/pdelab/finiteelementmap/pk1dbasis.hh>
#include <dune/pdelab/newton/newton.hh>

#include <dune/ax1/common/electrolyte.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/membrane1/membrane1_operator.hh>


template<class GV>
void membrane1_Pk (const GV& gv)
{
  
  // electrolyte definition #################
  
  Ion<double> sodium(1.0, "na");
  Ion<double> chloride(-1.0, "cl");
  
  Solvent<double> water(1.0);
  
  Electrolyte<double> elec1(water, 300.0, 1.0e8, 1);
  Electrolyte<double> elec2(water, 300.0, 1.0e8, 1);
  
  elec1.addIon(sodium);
  elec2.addIon(chloride);
  
  //#########################################
  
	// <<<1>>> Choose domain and range field type
	typedef typename GV::Grid::ctype Coord;
	typedef double Real;
	
	// <<<2>>> Make grid function space
	typedef Dune::PDELab::Pk1dLocalFiniteElementMap<Coord,Real> FEM;

	FEM fem(2);
	typedef Dune::PDELab::NoConstraints CON;
	typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
	typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
	GFS gfs(gv,fem);
	typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
	
	// <<<3>>> Make grid operator space
	typedef Membrane1LocalOperator LOP;                                     // <= NEW
	LOP lop(elec1, elec2, 4); //########
	typedef VBE::MatrixBackend MBE;
	typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
	GO go(gfs,gfs,lop);
	
	// <<<4>>> Select a linear solver backend
	typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
	LS ls(5000,true);
	
	// <<<5>>> solve nonlinear problem
	typedef typename GO::Traits::Domain U;
	//typedef typename Dune::PDELab::BackendVectorSelector<GFS,int>::Type U;

	U u(gfs,0.0); // initial value
	Dune::PDELab::Newton<GO,LS,U> newton(go,u,ls);                         // <= NEW
	newton.setReassembleThreshold(0.0);
	newton.setVerbosityLevel(2);
	newton.setReduction(1e-10);
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
	
	
	// <<<6>>> graphical output
	typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
	DGF udgf(gfs,u);
  
	Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
	vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
	vtkwriter.write("membrane1_Pk",Dune::VTK::OutputType::appendedraw);       // <= NEW

  std::valarray<typename DGF::Traits::RangeType> x;
  std::valarray<typename DGF::Traits::RangeType> f;
  //Tools::getSolutionVector(udgf, 5, x, f);
  //Tools::gnuplotSolution("membrane1_Pk.dat", x, f);
}

#endif
