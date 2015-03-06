#ifndef GEOGRIDFACTORY_HH
#define GEOGRIDFACTORY_HH

#include"gridvector.hh"

template<int dim>
class GeoGridFactory
{};

template<>
class GeoGridFactory<3>
{
  public:
  typedef Dune::UGGrid<3> Grid;
  typedef Dune::UGGrid<3>::LeafGridView GV;

  //constructor
  GeoGridFactory(GridVector<Grid> x_,GridVector<Grid> y_,GridVector<Grid> z_) : x(x_), y(y_), z(z_) {}

  GeoGridFactory(int holes,int level)
  {
    if (holes==1)
    {
      x.equidistant_n_xend(4,0.5);
      x.setEntityChange();
      x.geometric_n_h0_xend(15,0.125,3000);
 
      z.start(-4500);
      z.equidistant_n_xend(1,-4000);
      z.setEntityChange();
      z.equidistant_n_xend(1,-3800.5);
      z.setEntityChange();
      z.equidistant_n_xend(1,-3799.5);
      z.setEntityChange();
      z.equidistant_n_xend(1,-3300.5);
      z.setEntityChange();
      z.equidistant_n_xend(1,-3299.5);
      z.setEntityChange();
      z.equidistant_n_xend(1,-3100);
      z.setEntityChange();
      z.equidistant_n_xend(1,-3000);

      y=x;
    }
    if (holes==2)
    {
      x.start(-3000.0);
      x.geometric_n_hend_xend(5*(level+1),0.014,-1000.14);
      x.setEntityChange();
      x.equidistant_n_xend(level+1,-999.86);
      x.setEntityChange();
      x.geometric_n_h0_xend(5*(level+1),0.014,0);
      x.geometric_n_hend_xend(5*(level+1),0.014,999.86);
      x.setEntityChange();
      x.equidistant_n_xend(level+1,1000.14);
      x.setEntityChange();
      x.geometric_n_h0_xend(5*(level+1),0.014,3000);

      y.start(-2000.0);
      y.geometric_n_hend_xend(5*(level+1),0.014,-0.014);
      y.setEntityChange();
      y.equidistant_n_xend(level+1,0.014);
      y.setEntityChange();
      y.geometric_n_h0_xend(10,0.014,2000.0);
 
      z.start(-3900);
      z.equidistant_n_xend(2*(level+1),-3400);
      z.setEntityChange();
      z.equidistant_n_xend(level+1,-3399.86);
      z.setEntityChange();
      z.equidistant_n_xend(3*(level+1),-3101);
      z.setEntityChange();
      z.equidistant_n_xend(level+1,-3100);
      z.setEntityChange();
      z.equidistant_n_xend(level+1,-3000);
    }
  }  

  //make grid
  Grid* makeGrid()
  {
    //get GridFactory instance
    Dune::GridFactory<Dune::UGGrid<3> > factory;

    //insert vertices
    GridVector<Grid>::iterator itxend = x.end();
    GridVector<Grid>::iterator ityend = y.end();
    GridVector<Grid>::iterator itzend = z.end();
    for (GridVector<Grid>::iterator itz = z.begin(); itz!=itzend; ++itz)
      for (GridVector<Grid>::iterator ity = y.begin(); ity!=ityend; ++ity)
        for (GridVector<Grid>::iterator itx = x.begin(); itx!=itxend; ++itx)
        {
          //set vertex position
          Dune::FieldVector<Grid::ctype,3> pos;
          pos[0] = *itx;
          pos[1] = *ity;
          pos[2] = *itz;
          factory.insertVertex(pos);
        }

    //get geometry type
    const Dune::GeometryType gt(Dune::GeometryType::cube,3);
  
    //insert elements
    for (unsigned int i=0; i<x.size()-1; i++)
      for (unsigned int j=0; j<y.size()-1; j++)
        for (unsigned int k=0; k<z.size()-1; k++)
        {
          //set node positions
          std::vector<unsigned int> pos;
          pos.push_back(k*x.size()*y.size()+j*x.size()+i);
          pos.push_back(k*x.size()*y.size()+j*x.size()+i+1);
          pos.push_back(k*x.size()*y.size()+(j+1)*x.size()+i);
          pos.push_back(k*x.size()*y.size()+(j+1)*x.size()+i+1);
          pos.push_back((k+1)*x.size()*y.size()+j*x.size()+i);
          pos.push_back((k+1)*x.size()*y.size()+j*x.size()+i+1);
          pos.push_back((k+1)*x.size()*y.size()+(j+1)*x.size()+i);
          pos.push_back((k+1)*x.size()*y.size()+(j+1)*x.size()+i+1);
          factory.insertElement(gt,pos);
        }

    //make grid
    return factory.createGrid();
  }

  GridVector<Grid> x,y,z;
};

template<>
class GeoGridFactory<2>
{
  public:
  typedef Dune::UGGrid<2> Grid;
  typedef Dune::UGGrid<2>::LeafGridView GV;

  //constructor
  GeoGridFactory(GridVector<Grid> x_, GridVector<Grid> y_) : x(x_), y(y_) {}

  GeoGridFactory(int holes,int level)
  {
    if (holes==1)
    {
      x.equidistant_n_xend(4,0.5);
      x.setEntityChange();
      x.geometric_n_h0_xend(15,0.125,3000);
 
      y.start(-4500);
      y.equidistant_n_xend(1,-4000);
      y.setEntityChange();
      y.equidistant_n_xend(1,-3800.5);
      y.setEntityChange();
      y.equidistant_n_xend(1,-3799.5);
      y.setEntityChange();
      y.equidistant_n_xend(1,-3300.5);
      y.setEntityChange();
      y.equidistant_n_xend(1,-3299.5);
      y.setEntityChange();
      y.equidistant_n_xend(1,-3100);
      y.setEntityChange();
      y.equidistant_n_xend(1,-3000);
    }
    if (holes==2)
    {
      x.start(-3000.0);
      x.geometric_n_hend_xend(5*(level+1),0.014,-1000.14);
      x.setEntityChange();
      x.equidistant_n_xend(level+1,-999.86);
      x.setEntityChange();
      x.geometric_n_h0_xend(5*(level+1),0.014,0);
      x.geometric_n_hend_xend(5*(level+1),0.014,999.86);
      x.setEntityChange();
      x.equidistant_n_xend(level+1,1000.14);
      x.setEntityChange();
      x.geometric_n_h0_xend(5*(level+1),0.014,3000);

      y.start(-3900);
      y.equidistant_n_xend(2*(level+1),-3400);
      y.setEntityChange();
      y.equidistant_n_xend(3*(level+1),-3101);
      y.setEntityChange();
      y.equidistant_n_xend(level+1,-3100);
      y.setEntityChange();
      y.equidistant_n_xend(level+1,-3000);
    }
  }  
  //make grid
  Grid* makeGrid()
  {
    //get GridFactory instance
    Dune::GridFactory<Grid> factory;

    //insert vertices
    GridVector<Grid>::iterator itxend = x.end();
    GridVector<Grid>::iterator ityend = y.end();
    for (GridVector<Grid>::iterator ity = y.begin(); ity!=ityend; ++ity)
      for (GridVector<Grid>::iterator itx = x.begin(); itx!=itxend; ++itx)
      {
        //set vertex position
        Dune::FieldVector<Grid::ctype,2> pos;
        pos[0] = *itx;
        pos[1] = *ity;
        factory.insertVertex(pos);
      }
  
    //get geometry type
    const Dune::GeometryType gt(Dune::GeometryType::cube,2);

    //insert elements
    for (unsigned int i=0; i<x.size()-1;i++)
      for (unsigned int j=0; j<y.size()-1;j++)
      {
        //set node positions
        std::vector<unsigned int> pos;
        pos.push_back(i+j*x.size());
        pos.push_back(i+1+j*x.size());
        pos.push_back(i+(j+1)*x.size());
        pos.push_back(i+1+(j+1)*x.size());
        factory.insertElement(gt,pos);  
      } 

    //make grid
    return factory.createGrid();
  }

  GridVector<Grid> x,y;  
};

#endif
