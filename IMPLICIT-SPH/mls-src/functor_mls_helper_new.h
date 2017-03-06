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
      FunctorOuterHelper(PairIsph *isph, double *f=NULL)
        : FunctorOuter<PairIsph>(isph),_f(f),_f3(NULL),_is_entered(false) { }
      FunctorOuterHelper(PairIsph *isph, double **f3) 
        : FunctorOuter<PairIsph>(isph),_f(NULL),_f3(f3),_is_entered(false) { }
      FunctorOuterHelper(PairIsph *isph, double *f, double **f3) 
        : FunctorOuter<PairIsph>(isph),_f(f),_f3(f3),_is_entered(false) { }

      virtual void enterFor();
      virtual void operator()(const int ii);

      double* getScalarHelper();
      double* getVectorHelper(const int k);
      double* getScalarDevHelper(const int k);

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      double *_f,**_f3;
      Epetra_SerialDenseVector _q, _q_at_i, _q3_at_i[3];
      Epetra_SerialDenseMatrix _dq_df_at_i;

      bool _is_entered;
    };

    template<class PairIsph> inline double*
    FunctorOuterHelper<PairIsph>::getScalarHelper() {
      return _q_at_i.Values();
    }

    //TODO: change orrible name
    template<class PairIsph> inline double*
    FunctorOuterHelper<PairIsph>::getScalarDevHelper(const int jj) {
      return _dq_df_at_i[jj].Values();
    }

    template<class PairIsph> inline double*
    FunctorOuterHelper<PairIsph>::getVectorHelper(const int k) {
      return _q3_at_i[k].Values();
    }

    template<class PairIsph> inline void
    FunctorOuterHelper<PairIsph>::enterFor() {
      if (!_is_entered) {
        int ndof = _pair->basis->ndof(ISPH_MLS_MAX_ORDER);
        
        // workspace for j
        _q.Size(ndof);
        
        // scalar field workspace for i
        if (_f != NULL) 
          _q_at_i.Size(ndof);
        
        // vector field workspace for i
        if (_f3 != NULL) 
          for (int k=0;k<_dim;++k) 
            _q3_at_i[k].Size(ndof);
        
        _is_entered = true;
      }
    }

    template<class PairIsph> inline void
    FunctorOuterHelper<PairIsph>::operator()(const int ii) {
      if (!_is_entered)
        enterFor();

      int i = _ilist[ii];
      int itype = _type[i];

      // basis dof
      int np   = _pair->np;
      int ndof = _pair->basis->ndof(np);

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // pointers to q
      double *q = _q.Values();
      double *q_at_i = _q_at_i.Values();
      double *q3_at_i[3];
      for (int k=0;k<_dim;++k)
        q3_at_i[k] = _q3_at_i[k].Values();

      // initialization for early return
      if (_f != NULL)
        memset(q_at_i, 0, sizeof(double)*ndof);

      if (_f3 != NULL) 
        for (int k=0;k<_dim;++k)
          memset(q3_at_i[k], 0, sizeof(double)*ndof);
      
      if (_f != NULL && _f3 != NULL)
        if (_dq_df_at_i.ColDim() < (jnum+1))
          _dq_df_at_i.Shape(3,jnum+1);

      if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
        return;

      for (int jj=0;jj<jnum;++jj) {
        int j = (jlist[jj] & NEIGHMASK);
        int jtype = _type[j];

        if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype), 
                                             _pair->getParticleKind(jtype)))
          continue;
      
        // compute distance between partile i and j
        double rsq = 0.0, rij[3] = {};
        for (int k=0;k<_dim;++k) {
          rij[k] = _x[i][k] - _x[j][k];
          rsq += (rij[k]*rij[k]);
        }

        double r   = sqrt(rsq) + ISPH_EPSILON;
        double rth = LOOKUP(_pair->h, itype, jtype); 
        if (r < rth) {

          // compute kernel
          double w = _pair->kernel->val(r, rth);

          // compute basis
          _pair->basis->val(np, rij, rth, q); 

          // scalar field 
          if (_f != NULL) 
            for (int k=0;k<ndof;++k)
              q_at_i[k] += q[k]*w*_f[j];

          // vector field; component-wise
          if (_f3 != NULL) 
            for (int k2=0;k2<_dim;++k2) 
              for (int k1=0;k1<ndof;++k1) 
                q3_at_i[k2][k1] += q[k1]*w*_f3[j][k2];

          // compute dq_df
          if (_f != NULL && _f3 != NULL)
            for (int k=0;k<ndof;++k)
              _dq_df_at_i[jj][k] = q[k]*w;
        }
      }

      // self contribution
      {
        double r   = 0.0;
        double rij[3] = {};
        double rth = LOOKUP(_pair->h, itype, itype); 

        // compute kernel
        double w = _pair->kernel->val(r, rth);

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
              q3_at_i[k2][k1] += q[k1]*w*_f3[i][k2];

        // compute dq_df
        if (_f != NULL && _f3 != NULL)
        for (int k=0;k<ndof;++k)
          _dq_df_at_i[jnum][k] = q[k]*w;

      }
      
    }
  }
}
#endif

