#pragma once
#ifndef __FUNCTOR_MLS_GRADIENT_COMPACT_POISSON_H__
#define __FUNCTOR_MLS_GRADIENT_COMPACT_POISSON_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"
#include "functor_mls_helper_compact_poisson.h"

namespace LAMMPS_NS {

  namespace MLS {

    template<class PairIsph>
    class FunctorOuterGradientCompactPoisson : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterGradientCompactPoisson(PairIsph *isph, 
                                         double *u, 
                                         double alpha, 
                                         double *f,
                                         double **normal, double *g,
                                         double **u_grad = NULL)
        :  FunctorOuter<PairIsph>(isph),
           _u_helper(isph, u),
           _p_helper(isph, f, normal, g),
           _alpha(alpha),
           _u_grad(u_grad) { }

      void operator()(const int ii);

      double* getScalarGradient();

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      FunctorOuterHelper<PairIsph> _u_helper;
      FunctorOuterHelperCompactPoisson<PairIsph> _p_helper;

      double _alpha;
      double _u_grad_at_i[3], **_u_grad;
    };

    template<class PairIsph> inline double*
    FunctorOuterGradientCompactPoisson<PairIsph>::getScalarGradient() {
      return _u_grad_at_i;
    }

    template<class PairIsph> inline void
    FunctorOuterGradientCompactPoisson<PairIsph>::operator()(const int ii) {
      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype); 

      // initialize with zero
      memset(_u_grad_at_i, 0, sizeof(double)*3);

      if (_filter != NULL && !_filter->yes(ikind))
        return;
      
      // compute helper; q_at_i = sum P(x_j) W_ij u_j
      if (_u_helper.hasScalarField() || _u_helper.hasVectorField())
        _u_helper(ii);

      if (_p_helper.hasScalarField() || _p_helper.hasVectorField())
        _p_helper(ii);

      // basis arguments
      const int np   = _pair->np;
      const int ndof = _pair->basis->ndof(np);
      const int ndofpl = ndof + (ikind == PairIsph::Boundary);

      const double rth = LOOKUP(_pair->h, itype, itype);

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
        row = &VIEW2(M, ndofpl, idx, 0);

        // scalar field gradient with scaling of alpha 
        {
          if (_u_helper.hasScalarField()) {
            qq = _u_helper.getScalarHelper();
            _u_grad_at_i[k2] += dq*util.dotVectors(ndofpl, row, ndofpl, qq, 1)*_alpha;
          }
          
          // penalty 
          if (_u_helper.hasScalarField()) {
            qq = _p_helper.getScalarHelper();
            _u_grad_at_i[k2] += dq*util.dotVectors(ndofpl, row, ndofpl, qq, 1)*_alpha;
          }
        }

      }

      // store gradient at i if grad is not null and computed
      if (_u_grad != NULL)
        memcpy(_u_grad[i], _u_grad_at_i, sizeof(double)*_dim);
    }

  }
}
#endif

