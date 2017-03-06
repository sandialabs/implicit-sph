#pragma once
#ifndef __FUNCTOR_MLS_CURLCURL_H__
#define __FUNCTOR_MLS_CURLCURL_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"
#include "functor_mls_helper.h"

namespace LAMMPS_NS {

  namespace MLS {

    template<class PairIsph>
    class FunctorOuterCurlCurl : public FunctorOuterHelper<PairIsph> {
    public:
      FunctorOuterCurlCurl(PairIsph *isph, 
                           double **f,  
                           double alpha = 1.0,
                           double **curlcurl = NULL)
        : FunctorOuterHelper<PairIsph>(isph, f),
          _alpha(alpha),
          _curlcurl(curlcurl) { }

      void operator()(const int ii);

      double* getCurlCurl();

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      double _alpha,**_curlcurl;
      double _curlcurl_at_i[3];
    };

    template<class PairIsph> inline double*
    FunctorOuterCurlCurl<PairIsph>::getCurlCurl() {
      return _curlcurl_at_i;
    }

    template<class PairIsph> inline void
    FunctorOuterCurlCurl<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];

      // initialize
      memset(_curlcurl_at_i, 0, sizeof(double)*3);

      if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
        return;
      
      // compute helper; q_at_i = sum P(x_j) W_ij f_j
      FunctorOuterHelper<PairIsph>::operator()(ii); 

      // basis arguments
      int np   = _pair->np;
      int ndof = _pair->basis->ndof(np);
      double rth = LOOKUP(_pair->h, itype, itype);

      const int p[3][4][4] = { { {1, /**/ 1,1,0},      // Fy,xy  --
                                 {0, /**/ 0,2,0},      // Fx,yy  --
                                 {0, /**/ 0,0,2},      // Fx,zz
                                 {2, /**/ 1,0,1} },    // Fz,xz
                               { {2, /**/ 0,1,1},      // Fz,yz
                                 {1, /**/ 0,0,2},      // Fy,zz
                                 {1, /**/ 2,0,0},      // Fy,xx  --
                                 {0, /**/ 1,1,0} },    // Fx,xy  --
                               { {0, /**/ 1,0,1},      // Fx,xz
                                 {2, /**/ 2,0,0},      // Fz,xx
                                 {2, /**/ 0,2,0},      // Fz,yy
                                 {1, /**/ 0,1,1} } };  // Fy,yz
      
      // mass matrix
      double *M = _pair->M[i];

      // compute laplacian
      for (int k=0;k<_dim;++k) {
        int idx = 0;
        double tmp[4];

        for (int l=0;l<4;++l) {
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
        _curlcurl_at_i[k] = ((tmp[0] - tmp[1]) - (tmp[2] - tmp[3]))*_alpha;
      }

      if (_curlcurl != NULL)
        memcpy(_curlcurl[i], _curlcurl_at_i, sizeof(double)*_dim);

    }

  }
}
#endif

