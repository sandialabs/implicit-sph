#pragma once
#ifndef __FUNCTOR_MLS_GRADIENT_H__
#define __FUNCTOR_MLS_GRADIENT_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"
#include "functor_mls_helper.h"

namespace LAMMPS_NS {

  namespace MLS {

    template<class PairIsph>
    class FunctorOuterGradient : public FunctorOuterHelper<PairIsph> {
    public:
      FunctorOuterGradient(PairIsph *isph, 
                           double *f, 
                           double alpha = 1.0, 
                           double **grad = NULL)
        : FunctorOuterHelper<PairIsph>(isph, f),
          _alpha(alpha),
          _grad(grad) { }
      FunctorOuterGradient(PairIsph *isph, 
                           double **f3, 
                           double alpha = 1.0) 
        : FunctorOuterHelper<PairIsph>(isph, f3),
          _alpha(alpha),
          _grad(NULL) { }
      FunctorOuterGradient(PairIsph *isph, 
                           double *f, 
                           double **f3, 
                           double alpha = 1.0) 
        : FunctorOuterHelper<PairIsph>(isph, f, f3),
          _alpha(alpha),
          _grad(NULL) { }

      void operator()(const int ii);

      double* getScalarGradient();
      double* getVectorGradient(const int k = 0);

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      double _alpha,**_grad;
      double _f_grad_at_i[3],_f3_grad_at_i[3][3];
    };

    template<class PairIsph> inline double*
    FunctorOuterGradient<PairIsph>::getScalarGradient() {
      return _f_grad_at_i;
    }
    template<class PairIsph> inline double*
    FunctorOuterGradient<PairIsph>::getVectorGradient(const int k) {
      return _f3_grad_at_i[k];
    }

    template<class PairIsph> inline void
    FunctorOuterGradient<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];

      // initialize with zero
      memset(_f3_grad_at_i, 0, sizeof(double)*9);
      memset(_f_grad_at_i, 0, sizeof(double)*3);

      if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
        return;
      
      // compute helper; q_at_i = sum P(x_j) W_ij f_j
      FunctorOuterHelper<PairIsph>::operator()(ii); 

      // basis arguments
      int np   = _pair->np;
      int ndof = _pair->basis->ndof(np);
      double rth = LOOKUP(_pair->h, itype, itype);

      const int p[3][3] = { {1,0,0},    // dx
                            {0,1,0},    // dy
                            {0,0,1} };  // dz

      // mass matrix
      double *M = _pair->M[i];

      // compute gradient
      for (int k2=0;k2<_dim;++k2) {
        int idx = 0;
        double dq = 0.0, *qq, *row = NULL;

        // d/d{x,y,z}, pick a row correspond to the combination of p
        _pair->basis->dval(np, rth,   p[k2][0], p[k2][1], p[k2][2],   idx, dq); 
        row = &VIEW2(M, ndof, idx, 0);
        
        // scalar field gradient with scaling of alpha 
        qq = FunctorOuterHelper<PairIsph>::getScalarHelper();
        if (qq != NULL)
          _f_grad_at_i[k2] = dq*util.dotVectors(ndof, row, ndof, qq, 1)*_alpha;

        // vector field gradient with scaling of alpha;
        // check the first component of the vector field
        qq = FunctorOuterHelper<PairIsph>::getVectorHelper(0); 
        if (qq != NULL)
          for (int k1=0;k1<_dim;++k1) {
            qq = FunctorOuterHelper<PairIsph>::getVectorHelper(k1);
            _f3_grad_at_i[k1][k2] = dq*util.dotVectors(ndof, row, ndof, qq, 1)*_alpha;
          }
      }

      // store gradient at i if grad is not null and computed
      if (_grad != NULL && 
          FunctorOuterHelper<PairIsph>::getScalarHelper() != NULL) 
        memcpy(_grad[i], _f_grad_at_i, sizeof(double)*_dim);
    }

  }
}
#endif

