template<class GV>
void example02_Q1 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem;
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
  typedef Example02LocalOperator<BCTypeParam> LOP;                        // operator including boundary
  LOP lop( bctype );
  typedef VBE::MatrixBackend MBE;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop);

  // <<<4>>> Make FE function extending Dirichlet boundary conditions
  typedef typename GO::Traits::Domain U;
  U u(gfs,0.0);
  typedef BCExtension<GV,Real> G;                               // boundary value + extension
  G g(gv);
  Dune::PDELab::interpolate(g,gfs,u);                           // interpolate coefficient vector

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);

  // <<<6>>> assemble and solve linear problem
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
  SLP slp(go,u,ls,1e-10);
  slp.apply();

  // <<<7>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  DGF udgf(gfs,u);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::DataMode::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
  vtkwriter.write("example02_Q1",Dune::VTK::OutputType::appendedraw);
}
