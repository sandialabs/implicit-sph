#pragma once
#ifndef __FUNCTOR_MLS_HELPER_H__
#define __FUNCTOR_MLS_HELPER_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace MLS {

    // this compute temporary vector for evluating \sum P(x_j) W_ij f_j
    template<class PairIsph>
    class FunctorOuterHelper : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterHelper(PairIsph *isph, double *f)
        : FunctorOuter<PairIsph>(isph),
          _f(f),
          _f3(NULL),
          _is_entered(false) { }
      FunctorOuterHelper(PairIsph *isph, double **f3) 
        : FunctorOuter<PairIsph>(isph),
          _f(NULL),
          _f3(f3),
          _is_entered(false) { }
      FunctorOuterHelper(PairIsph *isph, double *f, double **f3) 
        : FunctorOuter<PairIsph>(isph),
          _f(f),
          _f3(f3),
          _is_entered(false) { }

      virtual void enterFor();
      virtual void operator()(const int ii);

      double* getScalarHelper();
      double* getVectorHelper(const int k);

      bool hasScalarField() const;  
      bool hasVectorField() const; 

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      double *_f,**_f3;

      int _lda;
      Epetra_SerialDenseVector _q, _q_at_i, _q3_at_i;

      bool _is_entered;
    };

    template<class PairIsph> inline double*
    FunctorOuterHelper<PairIsph>::getScalarHelper() {
      return &_q_at_i[0];
    }
    template<class PairIsph> inline double*
    FunctorOuterHelper<PairIsph>::getVectorHelper(const int k) {
      return &_q3_at_i[k*_lda];
    }
    template<class PairIsph> inline bool
    FunctorOuterHelper<PairIsph>::hasScalarField() const {
      return (_f != NULL);
    }
    template<class PairIsph> inline bool
    FunctorOuterHelper<PairIsph>::hasVectorField() const {
      return (_f3 != NULL);
    }

    template<class PairIsph> inline void
    FunctorOuterHelper<PairIsph>::enterFor() {
      if (!_is_entered) {
        // ndof of basis functions
        _lda = _pair->basis->ndof(ISPH_MLS_MAX_ORDER);
        
        // workspace for j
        _q.Size(_lda);
        
        // scalar field workspace for i
        if (_f != NULL) 
          _q_at_i.Size(_lda);
        
        // vector field workspace for i
        if (_f3 != NULL) 
          _q3_at_i.Size(_lda*_dim);
        
        _is_entered = true;
      }
    }

    template<class PairIsph> inline void
    FunctorOuterHelper<PairIsph>::operator()(const int ii) {
      if (!_is_entered)
        enterFor();

      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype);

      // basis dof
      const int np   = _pair->np;
      const int ndof = _pair->basis->ndof(np);

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // pointers to q
      double *q = _q.Values();
      double *q_at_i = _q_at_i.Values();
      double *q3_at_i = _q3_at_i.Values();

      // f at i: offset value when the interpolation property is set
      double f_at_i = 0.0;
      double f3_at_i[3] = {};

      // initialization for early return
      if (_f != NULL) {
        memset(q_at_i, 0, sizeof(double)*_lda);
        if (_pair->is_interpolation_enabled)
          f_at_i = _f[i];
      }
      if (_f3 != NULL) {
        memset(q3_at_i, 0, sizeof(double)*_lda*_dim);
        if (_pair->is_interpolation_enabled)
          memcpy(f3_at_i, _f3[i], sizeof(double)*_dim);
      }
      
      if (_filter != NULL && !_filter->yes(ikind))
        return;

      for (int jj=0;jj<jnum;++jj) {
        const int j = (jlist[jj] & NEIGHMASK);
        const int jtype = _type[j];
        const int jkind = _pair->getParticleKind(jtype);

        if (_filter != NULL && !_filter->yes(ikind, jkind))
          continue;
      
        // compute distance between partile i and j
        double rsq = 0.0, rij[3] = {};
        for (int k=0;k<_dim;++k) {
          rij[k] = _x[i][k] - _x[j][k];
          rsq += (rij[k]*rij[k]);
        }

        const double r   = sqrt(rsq) + ISPH_EPSILON;
        const double rth = LOOKUP(_pair->h, itype, jtype); 
        if (r < rth) {
          // compute kernel
          const double w = _pair->kernel->val(r, rth);

          // ** compute P W f_j
          {
            // compute basis
            memset(q, 0, sizeof(double)*ndof);
            _pair->basis->val(np, rij, rth, q); 
            
            // scalar field 
            if (_f != NULL) 
              for (int k=0;k<ndof;++k)
                q_at_i[k] += q[k]*w*(_f[j] - f_at_i);
            
            // vector field; component-wise
            if (_f3 != NULL) 
              for (int k2=0;k2<_dim;++k2) 
                for (int k1=0;k1<ndof;++k1) 
                  VIEW2(q3_at_i, _lda, k1, k2) += q[k1]*w*(_f3[j][k2] - f3_at_i[k2]);
          }
        }
      }

      // apply self contribution if the interpolation property is not set
      if (!_pair->is_interpolation_enabled) {
        const double r   = 0.0;
        const double rij[3] = {};
        const double rth = LOOKUP(_pair->h, itype, itype); 

        // compute kernel
        const double w = _pair->kernel->val(r, rth);

        // compute basis
        _pair->basis->val(np, rij, rth, q); 

        // scalar field 
        if (_f != NULL) 
          for (int k=0;k<ndof;++k)
            q_at_i[k] += q[k]*w*_f[i];
        
        // vector field; component-wise
        if (_f3 != NULL) 
          for (int k2=0;k2<_dim;++k2) 
            for (int k1=0;k1<ndof;++k1) 
              VIEW2(q3_at_i, _lda, k1, k2) += q[k1]*w*_f3[i][k2];
      }

    }
  }
}
#endif

