#pragma once
#ifndef __FUNCTOR_MLS_CURL_H__
#define __FUNCTOR_MLS_CURL_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"
#include "functor_mls_helper.h"

namespace LAMMPS_NS {

  namespace MLS {

    template<class PairIsph>
    class FunctorOuterCurl : public FunctorOuterHelper<PairIsph> {
    public:
      FunctorOuterCurl(PairIsph *isph, 
                       double **f,  
                       double alpha = 1.0,
                       double **curl = NULL)
        : FunctorOuterHelper<PairIsph>(isph, f),
          _alpha(alpha),
          _curl(curl) { }

      void operator()(const int ii);

      double* getCurl();

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      double _alpha,**_curl;
      double _curl_at_i[3];
    };

    template<class PairIsph> inline double*
    FunctorOuterCurl<PairIsph>::getCurl() {
      return _curl_at_i;
    }

    template<class PairIsph> inline void
    FunctorOuterCurl<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];

      // initialize
      memset(_curl_at_i, 0, sizeof(double)*3);

      if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
        return;
      
      // compute helper; q_at_i = sum P(x_j) W_ij f_j
      FunctorOuterHelper<PairIsph>::operator()(ii); 

      // basis arguments
      int np   = _pair->np;
      int ndof = _pair->basis->ndof(np);
      double rth = LOOKUP(_pair->h, itype, itype);

      const int p[3][2][4] = { { {2, /**/ 0,1,0},      // Fz,y
                                 {1, /**/ 0,0,1} },    // Fy,z
                               { {0, /**/ 0,0,1},      // Fx,z
                                 {2, /**/ 1,0,0} },    // Fz,x
                               { {1, /**/ 1,0,0},      // Fy,x
                                 {0, /**/ 0,1,0} } };  // Fx,y
      
      // mass matrix
      double *M = _pair->M[i];

      // compute laplacian
      for (int k=0;k<3;++k) {
        int idx = 0;
        double tmp[2];

        for (int l=0;l<2;++l) {
          // d/d{x,y,z}, pick a row correspond to the combination of p
          if (_dim == 2 && (p[k][l][0] == 2 || p[k][l][3] > 0)) {
            tmp[l] = 0.0;
          } else {
            double dq = 0.0, *qq = NULL, *row = NULL;
            _pair->basis->dval(np, rth,   p[k][l][1], p[k][l][2], p[k][l][3],   idx, dq); 
            row = &VIEW2(M, ndof, idx, 0);

            qq = FunctorOuterHelper<PairIsph>::getVectorHelper(p[k][l][0]);
            tmp[l] = dq*util.dotVectors(ndof, row, ndof, qq, 1);
          }
        }
        _curl_at_i[k] = (tmp[0] - tmp[1])*_alpha;
      }

      if (_curl != NULL)
        memcpy(_curl[i], _curl_at_i, sizeof(double)*3);

    }

  }
}
#endif

