#pragma once
#ifndef __FUNCTOR_MLS_GRADIENT_OPERATOR_H__
#define __FUNCTOR_MLS_GRADIENT_OPERATOR_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"
#include "functor_mls_helper.h"

namespace LAMMPS_NS {

  namespace MLS {

    template<class PairIsph>
    class FunctorOuterGradientOperator : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterGradientOperator(PairIsph *isph, double alpha = 1.0) 
        : FunctorOuter<PairIsph>(isph),
          _alpha(alpha),
          _cnt(0),
          _lda(0),
          _is_entered(false) { }

      void enterFor();
      void operator()(const int ii);

      int getSize() const;
      int getNeighbor(const int jj);
      double getValue(const int jj, const int idx);

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      double _alpha;

      int _cnt, _lda;
      Epetra_IntSerialDenseVector _idx;
      Epetra_SerialDenseVector _val, _q;

      bool _is_entered;
    };


    template<class PairIsph>
    inline int
    FunctorOuterGradientOperator<PairIsph>::getSize() const {
      return _cnt;
    }
    template<class PairIsph>
    inline int
    FunctorOuterGradientOperator<PairIsph>::getNeighbor(const int jj) {
      return _idx[jj];
    }
    template<class PairIsph>
    inline double
    FunctorOuterGradientOperator<PairIsph>::getValue(const int jj, const int k) {
      return _val[jj + k*_lda];
    }

    template<class PairIsph>
    inline void
    FunctorOuterGradientOperator<PairIsph>::enterFor() {
      if (!_is_entered) {
        
        _lda = _pair->tags_in_cut->MaxNumIndices();

        _idx.Size(_lda);
        _val.Size(_lda*3);

        _q.Size(_pair->basis->ndof(ISPH_MLS_MAX_ORDER) + 2);

        _is_entered = true;
      }
    }

    template<class PairIsph> inline void
    FunctorOuterGradientOperator<PairIsph>::operator()(const int ii) {
      if (!_is_entered)
        enterFor();

      const int i = _ilist[ii];
      const int itype = _type[i];

      // extract pointers
      double *val = _val.Values();
      int *idx = _idx.Values();

      // clean up the work arrays
      memset(idx, 0, sizeof(int)*_lda);
      memset(val, 0, sizeof(double)*_lda*3);

      _cnt = 0;

      if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
        return;

      // basis dof
      const int np   = _pair->np;
      const int ndof = _pair->basis->ndof(np) + _pair->isCompactPoisson();
      
      // extract pointer to pre-allocated basis
      double *qq = _q.Values();

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // index to compute derivatives
      const int p[3][3] = { {1,0,0},    // dx
                            {0,1,0},    // dy
                            {0,0,1} };  // dz
      
      // mass matrix
      double *M = _pair->M[i];
      
      // self contribution
      if (_pair->is_interpolation_enabled) {
        idx[_cnt++] = i;
        // the value is initialized as zero and later sum -aij
      } else {
        const double rij[3] = {};
        const double r   = 0.0;
        const double rth = LOOKUP(_pair->h, itype, itype);

        // compute basis and kernel values
        _pair->basis->val(np, rij, rth, qq);
        const double wij = _pair->kernel->val(r,rth);

        // compute laplacian                                                                                  
        for (int k=0;k<_dim;++k) {
          // compute derivatives                                                                              
          int idx = 0; double dq = 0.0;
          _pair->basis->dval(np, rth,   p[k][0], p[k][1], p[k][2],   idx, dq);

          // grab a row corresponding to the basis->dval                                                      
          double *row = &VIEW2(M, ndof, idx, 0);

          // compute operator for (i,j)
          const double ijtmp = dq*util.dotVectors(ndof, row, ndof, qq, 1);
          VIEW2(val, _lda, _cnt, k) += ijtmp*wij;
        }
        idx[_cnt++] = i;
      }

      // compute operator
      for (int jj=0;jj<jnum;++jj) {
        // j th particle; always put neighmask two bits are reserved for LAMMPS
        const int j = (jlist[jj] & NEIGHMASK);
        const int jtype = _type[j];

        // compute distance between partile i and j
        double rsq = 0.0, rij[3] = {};
        for (int k=0;k<_dim;++k) {
          rij[k] = _x[i][k] - _x[j][k];
          rsq += (rij[k]*rij[k]);
        }

        const double r   = sqrt(rsq) + ISPH_EPSILON;
        const double rth = LOOKUP(_pair->h, itype, jtype);

        if (r < rth) {
          if (_filter != NULL && _filter->yes(_pair->getParticleKind(itype),
                                              _pair->getParticleKind(jtype))) {
            
            // compute basis and kernel value
            _pair->basis->val(np, rij, rth, qq);
            const double wij = _pair->kernel->val(r,rth);

            // compute operator
            for (int k=0;k<_dim;++k) {
              // compute derivatives
              int idx = 0; double dq = 0.0;
              _pair->basis->dval(np, rth,   p[k][0], p[k][1], p[k][2],   idx, dq);

              // grab a row corresponding to the basis->dval
              double *row = &VIEW2(M, ndof, idx, 0);

              // compute operator for (i,j)
              const double ijtmp = dq*util.dotVectors(ndof, row, ndof, qq, 1)*wij;
              VIEW2(val, _lda, _cnt, k) += ijtmp;

              // substract aij from aii
              if (_pair->is_interpolation_enabled)
                VIEW2(val, _lda, 0, k) -= ijtmp;
            }
            _idx[_cnt++] = j;
          }
        }
      }
      
      _val.Scale(_alpha);
    }
  }
}
#endif

