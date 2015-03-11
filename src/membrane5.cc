// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file

    \brief Solve elliptic problem in unconstrained spaces with conforming finite elements
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include <dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>

#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#include <dune/grid/onedgrid.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_ALUGRID
#include<dune/grid/alugrid.hh>
#include<dune/grid/io/file/dgfparser/dgfalu.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif

#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/p12dfem.hh>
#include<dune/pdelab/finiteelementmap/pk2dfem.hh>
#include<dune/pdelab/finiteelementmap/pk3dfem.hh>
#include<dune/pdelab/finiteelementmap/q12dfem.hh>
#include<dune/pdelab/finiteelementmap/q22dfem.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/finiteelementmap/p1fem.hh>
#include<dune/pdelab/finiteelementmap/rannacher_turek2dfem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/constraints.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>

#include<dune/ax1/common/tools.hh>
#include<dune/ax1/common/output.hh>
#include<dune/ax1/common/membrane_physics.hh>

#include<dune/ax1/membrane4/membrane4_bctype.hh>
#include<dune/ax1/membrane4/membrane4_bcextension.hh>
#include<dune/ax1/membrane4/membrane4_Pk.hh>

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
	std::cout<< "This is Membrane5. Dick." << std::endl;
  try
  {
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
    {
      std::cout<< "This is a sequential program." << std::endl;
    } else {
      if(helper.rank()==0)
        std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
    }
    if (argc<3)
    {
      if(helper.rank()==0)
        std::cout << "usage: ./membrane3 <min_level> <max_level>" << std::endl;
      return 1;
    }
    
    // level loop #####################################################
    
    const int minLevel = atoi(argv[1]);
    const int maxLevel = atoi(argv[2]);
    
    std::valarray<double> dof(maxLevel-minLevel), err(maxLevel-minLevel);
    
    /*
    Output::gnuplotInitialize("error.dat");
    Output::gnuplotInitialize("error_plot.dat");
    Output::gnuplotInitialize("rel_error_plot.dat");
    */
    
    for ( int level=minLevel; level<maxLevel; level++ )
    {
      //sscanf(argv[1],"%d",&level);
      
      std::cout << "" << std::endl;
      std::cout << "##### convergence loop, level: " << level << "#####" << std::endl;
      
      const int dim = 1;
      
      // sequential version
      if (helper.size()==1)
      {
        
        // electrolyte definition ###########################
        
        Ion<double> na_in(12.0,"na_in");
        Ion<double> k_in(155.0,"k_in");
        Ion<double> ca_in(167.0265,"ca_in");

        Ion<double> na_ex(145.0,"na_ex");
        Ion<double> k_ex(4.0,"k_ex");
        Ion<double> ca_ex(149.0,"ca_ex");
        
        Solvent<double> water(80.0);
        
        //Electrolyte<double> elec1(water, 300.0, 1.0e8, 1);
        //Electrolyte<double> elec2(water, 300.0, 1.0e8, 1);
        Electrolyte<double> elec1(water, 279.45, con_mol, 1e-9);
        Electrolyte<double> elec2(water, 279.45, con_mol, 1e-9);
        
        elec1.addIon(na_in);
        elec1.addIon(k_in);
        elec1.addIon(ca_in);
        
        elec2.addIon(na_ex);
        elec2.addIon(k_ex);
        elec2.addIon(ca_ex);
        
        double d = 10.0;
        Membrane<double> memb( d );
        
        double dl = elec1.getDebyeLength();
        
        double dPhi = 4.0;
        
        Physics<double> physics( elec1, memb, elec2 );
        
        std::cout << "thickness of membrane: " << d << std::endl;
        
        // grid generation ###################################
        
        /*
         Dune::FieldVector<double,dim> L(1.0);
         Dune::FieldVector<int,dim> N(1);
         Dune::FieldVector<bool,dim> periodic(false);
         int overlap=0;
         Dune::YaspGrid<dim> grid(L,N,periodic,overlap);
         */
        
        //std::vector<double> coords = {-5., -4., -3., -2., -1., 0., 1., 2., 3., 4., 5.};
        double x_l = -50.0;
        double x_r = 50.0;
        
        std::vector<double> coords = { x_l, -0.5*d };
        
        Tools::globalRefineVector( coords, level );
        
        //int currentSize = coords.size();
        
        /*
        std::vector<double> coordsMemb = { - 0.5 * d, 0.5 * d };
        Tools::globalRefineVector( coordsMemb, level );
        for (int i=1; i<coordsMemb.size()-1; i++) coords.push_back( coordsMemb[i] );
        */
        std::vector<double> coords2 = { 0.5*d, x_r };
        Tools::globalRefineVector( coords2, level );
        
        for (int i=0; i<coords2.size(); i++) coords.push_back( coords2[i] ); // mirror grid
        
        //for (int i=0; i<coords.size(); i++) std::cout << coords[i] << std::endl;
        
        //exit(0);
        
        Dune::OneDGrid grid(coords);
        
        //grid.globalRefine(level);
        //typedef Dune::YaspGrid<dim>::LeafGridView GV;
        typedef Dune::OneDGrid::LeafGridView GV;
        const GV& gv=grid.leafView();
        
        dof[level-minLevel] = gv.size(1);
        
        // analytical solutions ##############################
        
        double boundl = 0.0;
        double boundr = 0.231;
        
        physics.setBoundaries( boundl, boundr );
        std::cout << "outer boundaries set to: " << boundl << " and " << boundr << std::endl;
        
        //std::cout << "#bla# " << physics.potential2( 0.5 * d, dl,  1,  d, dPhi ) << std::endl;
        
        physics.setMembranePos( -0.5*d, 0.5*d );
        
        // computation #######################################
        
        std::valarray<double> position, solution;
        
        //############################################
        membrane4_Pk(gv, physics, position, solution);
        //############################################
        
        int half = ( solution.size() + 0 ) / 2;
        
        std::cout << " slopes: " << 0.5 * dPhi / d << " " << 
        ( solution[half+2] - solution[half+1] ) / ( position[half+2] - position[half+1] ) << std::endl;
        
        std::valarray<double> error(solution), concentrationl(half), concentrationr(half),
                                                                positionl(half), positionr(half);
        
        std::vector<std::valarray<double> > concL, concR;
        
        for ( int i=0; i<physics.elecl.numOfSpecies(); i++ ) concL.push_back( concentrationl );
        for ( int i=0; i<physics.elecr.numOfSpecies(); i++ ) concR.push_back( concentrationr );
        
        for ( int i=0; i<position.size(); i++ )
        {
          if ( position[i] <= physics.getMembranePosL() )
          {
            //error[i] = solution[i] - physics.potential2( position[i], dl, -1, -d, -4.0/3.0*dPhi );
            for (int j=0; j<physics.elecl.numOfSpecies(); j++ )
              concL[j][i] = physics.elecl.getConcentration (j, solution[i]);
            positionl[i] = position[i];
          }
          else if ( position[i] >= physics.getMembranePosR() )
          {
            //error[i] = solution[i] - physics.potential2( position[i], dl,  2,  d,  2.0/3.0*dPhi );
            for (int j=0; j<physics.elecr.numOfSpecies(); j++ )
              concR[j][i-half] = physics.elecr.getConcentration (j, solution[i]);
            positionr[i-half] = position[i];
          }
          else
            error[i] = 0.0;
        }
        
        //err[level-minLevel] = (abs(error)).max();
        
        // graphical output ##################################
        
        //Output::gnuplotAppend("error.dat", (double) dof[level-minLevel], err[level-minLevel] );
        
        Output::gnuplotArray("solution.dat", position, solution);
        
        /*
        std::stringstream infoStream;
        infoStream << "dof: " << dof[level-minLevel];
        Output::gnuplotAppendArray("error_plot.dat", position, error, infoStream.str());
        error /= (abs(error)).max();
        Output::gnuplotAppendArray("rel_error_plot.dat", position, error, infoStream.str());
        */
        
        Output::gnuplotMultiArray("concentrationl.dat", positionl, concL);
        Output::gnuplotMultiArray("concentrationr.dat", positionr, concR);
      }
    } // level loop
    
    /*
    std::valarray<double> mdof(maxLevel-minLevel-1), ord(maxLevel-minLevel-1);
    
    Tools::convergenceOrder( dof, err, mdof, ord );
    
    Output::gnuplotArray("order.dat", mdof, ord );
     */
  }
  
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
	return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
	return 1;
  }
}
