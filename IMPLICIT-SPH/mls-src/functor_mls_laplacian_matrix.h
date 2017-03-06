#pragma once
#ifndef __FUNCTOR_MLS_LAPLACIAN_MATRIX_H__
#define __FUNCTOR_MLS_LAPLACIAN_MATRIX_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace MLS {

    template<class PairIsph>
    class FunctorOuterLaplacianMatrix : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterLaplacianMatrix(PairIsph *isph, 
                                  double alpha, 
                                  double *material = NULL)
        : FunctorOuter<PairIsph>(isph),
          _alpha(alpha),
          _material(material) 
      { }

      void enterFor();
      void operator()(const int ii);
      void exitFor();

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      Epetra_IntSerialDenseVector _idx;
      Epetra_SerialDenseVector _val, _q;

      double _alpha, *_material;
    };

    template<class PairIsph> inline void
    FunctorOuterLaplacianMatrix<PairIsph>::enterFor() {
      if (_filter == NULL)
        _pair->error->all(FLERR, "FunctorLaplacianMatrix:: Filter is required to consider dirichlet-like boundary");

      const int n = _pair->tags_in_cut->MaxNumIndices();
      _idx.Size(n);
      _val.Size(n);

      const int ndof = _pair->basis->ndof(ISPH_MLS_MAX_ORDER) + 2;
      _q.Size(ndof);
    }

    template<class PairIsph> inline void
    FunctorOuterLaplacianMatrix<PairIsph>::operator()(const int ii) {
      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype);

      // extract pre-allocated row
      double *val = _val.Values();
      int *idx = _idx.Values();

      if (!_filter->yes(ikind)) {
        idx[0] = _tag[i];

        // when we do not compute for i, leave the diagonal entree zero
        // and later need to change diagonal
        val[0] = 0.0; 
        _pair->A.crs->ReplaceGlobalValues(_tag[i], 1, val, idx);
        return;
      }

      // order corresponding to the particle
      const int np = _pair->np;
      const int ndof = _pair->basis->ndof(np) + (_pair->isCompactPoisson() & 
                                                 ikind == PairIsph::Boundary);

      // extract pointer to pre-allocated basis
      double *qq = _q.Values();
      memset(qq, 0, sizeof(double)*ndof);

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];
      int cnt = 0;

      // index to compute derivatives
      const int p[3][3] = { {2,0,0},    // dx2
                            {0,2,0},    // dy2
                            {0,0,2} };  // dz2

      // grab a mass matrix and initialize diag value
      double *M = _pair->M[i];
      double aii = 0.0;

      // compute a row in the laplacian matrix
      for (int jj=0;jj<jnum;++jj) {
        // j th particle; always put neighmask two bits are reserved for LAMMPS
        const int j = (jlist[jj] & NEIGHMASK);
        const int jtype = _type[j];
        const int jkind = _pair->getParticleKind(jtype);

        // compute distance between partile i and j
        double rsq = 0.0, rij[3] = {};
        for (int k=0;k<_dim;++k) {
          rij[k] = _x[i][k] - _x[j][k];
          rsq += (rij[k]*rij[k]);
        }

        const double r   = sqrt(rsq) + ISPH_EPSILON;
        const double rth = LOOKUP(_pair->h, itype, jtype);

        if (r < rth) {
          if (_filter->yes(ikind, jkind)) {
            // compute basis val
            _pair->basis->val(np, rij, rth, qq);

            const double wij = _pair->kernel->val(r, rth);

            // compute laplacian
            double aij = 0.0;

            for (int k=0;k<_dim;++k) {
              // compute derivatives
              int idx = 0; double dq = 0.0;
              _pair->basis->dval(np, rth,   p[k][0], p[k][1], p[k][2],   idx, dq);

              // grab a row corresponding to the basis->dval
              double *row = &VIEW2(M, ndof, idx, 0);

              // compute laplacian value for (i,j)
              aij += dq*util.dotVectors(ndof, row, ndof, qq, 1);
            }

            // multiply kernel value
            aij *= wij;
            aii -= aij;

            val[cnt] = aij;
            idx[cnt] = _tag[j];

            ++cnt;
          }
        }
      }

      // self contribution
      if (_pair->is_interpolation_enabled) {
        val[cnt] = aii;
        idx[cnt] = _tag[i];

        ++cnt;
      } else {
        const double rij[3] = {};
        const double r   = 0.0;
        const double rth = LOOKUP(_pair->h, itype, itype);

        // compute basis val
        _pair->basis->val(np, rij, rth, qq);

        const double wij = _pair->kernel->val(r, rth);

        // compute laplacian
        aii = 0.0;
        for (int k=0;k<_dim;++k) {
          // compute derivatives
          int idx = 0; double dq = 0.0;
          _pair->basis->dval(np, rth,   p[k][0], p[k][1], p[k][2],   idx, dq);

          // grab a row corresponding to the basis->dval
          double *row = &VIEW2(M, ndof, idx, 0);

          // compute laplacian value for (i,j)
          aii += dq*util.dotVectors(ndof, row, ndof, qq, 1);
        }

        // multiply kernel value
        aii *= wij; 

        val[cnt] = aii;
        idx[cnt] = _tag[i];

        ++cnt;
      }

      // scale the row by alpha and material
      _val.Scale(_alpha*(_material == NULL ? 1.0 : _material[i]));
      _pair->A.crs->SumIntoGlobalValues(_tag[i], cnt, val, idx);
    }

    template<class PairIsph> inline
    void FunctorOuterLaplacianMatrix<PairIsph>::exitFor() {

      // finalize connectivity construction
      if (_pair->A.is_filled == 0) {
        _pair->A.crs->FillComplete();
        _pair->A.crs->OptimizeStorage();
        _pair->A.is_filled = 1;
      }
    }

  }
}
#endif

