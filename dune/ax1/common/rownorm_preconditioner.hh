#ifndef _ROWNORM_PRECONDITIONER_HH
#define _ROWNORM_PRECONDITIONER_HH

#include <dune/common/float_cmp.hh>

#include <dune/ax1/common/constants.hh>

/*
 this class changes the matrix and the right side
 */
template<class GV, int BlockSize = 1>
class RowNormPreconditioner
{

  public:
    // Grid
    typedef typename GV::Grid Grid;
    typedef typename GV::IndexSet IndexSet;
    typedef GV GridView;
    typedef typename Grid::ctype ctype;
    enum
    {
      dim = Grid::dimension
    /* !<dimension of the grid we are using */};

    //  typedef UDGF UDGFunction;

    // Maximal number of shapefunctions ( on cube ) for velocity and
    // pressure respectively
    // enum{ n_sfs_v_max = UDGFunction::BaseType::n_sfs_v_max};
    // enum{ n_sfs_p_max = UDGFunction::BaseType::n_sfs_p_max};

    typedef Dune::FieldVector<ctype, BlockSize> LocalVectorBlock;
    typedef Dune::FieldMatrix<ctype, BlockSize, BlockSize> LocalMatrixBlock;

    typedef Dune::BCRSMatrix<LocalMatrixBlock> Matrix;
    typedef Dune::BlockVector<LocalVectorBlock> SolutionVector;

    typedef std::vector<size_t> MatrixIndexMap;

    typedef Dune::FieldVector<ctype, dim - 1> FacePoint;
    typedef Dune::FieldVector<ctype, dim> Point;

  private:

    //const GridView & gridview;

    // Mapper for grid entities
    template<int dim>
    struct MapperLayout
    {
        bool contains(Dune::GeometryType gt)
        {
          if (gt.dim() == dim) return true;
          return false;
        }
    };
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, MapperLayout> Mapper;
    typedef typename GridView::template Codim<0>::Iterator Iterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

  public:

    //! Perform equilibration of matrix rows
    RowNormPreconditioner(Matrix &A, SolutionVector &B)
    {

      typedef typename Matrix::RowIterator RowIterator;
      typedef typename Matrix::ColIterator ColIterator;
      typedef typename Matrix::size_type size_type;

      {

        RowIterator rit = A.begin();
        const RowIterator erit = A.end();
        for (; rit != erit; ++rit)
        {
          const size_type r_index = rit.index();

          LocalMatrixBlock & diag = (*rit)[r_index];
          LocalVectorBlock vnorm;

          //	ctype vnorm(0),pnorm(0);

          // Get equilibration coefficients (= row norms)
          getBlockNorm(diag, vnorm);

          //	std::cout << vnorm << " " << pnorm << std::endl;

          // Equilibrate row and rhs
          rowNormalizeBlock<true> (diag, vnorm);
          rowNormalizeRhs(B[r_index], vnorm);

          // Iterate off diagonal elements
          ColIterator cit = rit->begin();
          const ColIterator ecit = rit->end();

          for (; cit != ecit; ++cit)
          {

            const size_type c_index = cit.index();

            // Skipt diagonal
            if (c_index == r_index) continue;

            rowNormalizeBlock<false> (*cit, vnorm);

          }// cit
        }// rit
      }

    }

    /*
     This functions are used to scale the both rows in the jacobian matrix. The scale factor is
     the infinity norm in the every row => if the biggest value is on the diagonal,
     the diagonal in jacobian matrix consists of 1
     */
    void getBlockNorm(LocalMatrixBlock &m, LocalVectorBlock & norm)
    {
      typename LocalMatrixBlock::RowIterator fit = m.begin();
      typename LocalMatrixBlock::RowIterator efit = m.end();

      for (; fit != efit; ++fit)
        norm[fit.index()] = fit->infinity_norm();
    }

    template<bool diagonal>
    void rowNormalizeBlock(LocalMatrixBlock & m, const LocalVectorBlock & norm)
    {

      typename LocalMatrixBlock::Iterator fit = m.begin();
      typename LocalMatrixBlock::Iterator efit = m.end();

      for (; fit != efit; ++fit)
      {
        const ctype n = norm[fit.index()];
        if (n)
          *fit /= n;
        else if (diagonal) (*fit)[fit.index()] = 1.0;
      }
    }

    void rowNormalizeRhs(LocalVectorBlock & m, const LocalVectorBlock & norm)
    {

      typename LocalVectorBlock::Iterator fit = m.begin();
      typename LocalVectorBlock::Iterator efit = m.end();

      for (; fit != efit; ++fit)
      {
        const ctype n = norm[fit.index()];
        if (n) *fit /= n;
      }
    }

    /*
    //This functions are used to scale the second row in the jacobian matrix. The first row
    //stays unchanged. The diagonal elements in the submatrix are the same.
    void getBlockNorm(LocalMatrixBlock &m, LocalVectorBlock & norm)
    {
      typename LocalMatrixBlock::RowIterator fit = m.begin();
      typename LocalMatrixBlock::RowIterator efit = m.end();

      for (; fit != efit; ++fit)
      {
        norm[fit.index()] = fit->infinity_norm();

      }
    }

    template<bool diagonal>
    void rowNormalizeBlock(LocalMatrixBlock & m, const LocalVectorBlock & norm)
    {

      typename LocalMatrixBlock::Iterator fit = m.begin();
      typename LocalMatrixBlock::Iterator efit = m.end();

      for (; fit != efit; ++fit)
      {
        if (fit.index() == 0) continue;
        const ctype n = norm[fit.index()] / norm[0];
        if (n)
          *fit /= n;
        else if (diagonal) (*fit)[fit.index()] = 1.0;
      }
    }

    void rowNormalizeRhs(LocalVectorBlock & m, const LocalVectorBlock & norm)
    {

      typename LocalVectorBlock::Iterator fit = m.begin();
      typename LocalVectorBlock::Iterator efit = m.end();

      for (; fit != efit; ++fit)
      {
        if (fit.index() == 0) continue;
        const ctype n = norm[fit.index()] / norm[0];
        if (n) *fit /= n;
      }
    }
    */

    void getBlockNorm(LocalMatrixBlock &m, ctype & vnorm, ctype & pnorm)
    {
      typename LocalMatrixBlock::RowIterator fit = m.begin();
      typename LocalMatrixBlock::RowIterator efit = m.end();

      ctype *norm = &vnorm;

      for (; fit != efit; ++fit)
      {
        if (fit.index() == dim * 0) norm = &pnorm;

        const ctype & current = (*fit).two_norm();
        if (current > *norm)
        {
          *norm = current;
        }
      }
    }

    template<bool diagonal>
    void rowNormalizeBlock(LocalMatrixBlock & m, ctype & vn, ctype & pn)
    {

      typename LocalMatrixBlock::RowIterator fit = m.begin();
      typename LocalMatrixBlock::RowIterator efit = m.end();

      ctype * norm = &vn;

      for (; fit != efit; ++fit)
      {
        if (fit.index() == dim * 0) norm = &pn;

        if (*norm)
          (*fit) /= *norm;
        else if (diagonal) (*fit)[fit.index()] = 1.0;
      }

    }

    void rowNormalizeRhs(LocalVectorBlock & r, ctype & vn, ctype & pn)
    {

      typename LocalVectorBlock::Iterator fit = r.begin();
      typename LocalVectorBlock::Iterator efit = r.end();

      ctype * norm = &vn;

      for (; fit != efit; ++fit)
      {
        if (fit.index() == dim * 0) norm = &pn;

        if (*norm) (*fit) /= *norm;
      }

    }

};

#endif /* _ROWNORM_PRECONDITIONER_HH */
