#pragma once
#ifndef __FUNCTOR_MLS_DIVERGENCE_H__
#define __FUNCTOR_MLS_DIVERGENCE_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"
#include "functor_mls_helper.h"

namespace LAMMPS_NS {

  namespace MLS {

    template<class PairIsph>
    class FunctorOuterDivergence : public FunctorOuterHelper<PairIsph> {
    public:
      FunctorOuterDivergence(PairIsph *isph, double **f,  
                             double alpha, 
                             double *div = NULL)
        : FunctorOuterHelper<PairIsph>(isph, f),_alpha(alpha),_div(div) { }

      void operator()(const int ii);

      double getDivergence();

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      double _alpha,*_div;
      double _div_at_i;
    };

    template<class PairIsph> inline double
    FunctorOuterDivergence<PairIsph>::getDivergence() {
      return _div_at_i;
    }

    template<class PairIsph> inline void
    FunctorOuterDivergence<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];

      // initialize
      _div_at_i = 0.0;

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

      // compute divergence
      for (int k=0;k<_dim;++k) {
        int idx = 0;
        double dq = 0.0, *qq = NULL, *row = NULL;

        // d/d{x,y,z}, pick a row correspond to the combination of p1,p2,p3
        _pair->basis->dval(np, rth,   p[k][0], p[k][1], p[k][2],   idx, dq); 
        row = &VIEW2(M, ndof, idx, 0);
        
        // divergence
        qq = FunctorOuterHelper<PairIsph>::getVectorHelper(k);
        _div_at_i += dq*util.dotVectors(ndof, row, ndof, qq, 1);
      }

      // scale the divergence value
      _div_at_i *= _alpha;

      // store divergence at i if grad is not null and computed
      if (_div != NULL)
        _div[i] = _div_at_i;
    }

  }
}
#endif

