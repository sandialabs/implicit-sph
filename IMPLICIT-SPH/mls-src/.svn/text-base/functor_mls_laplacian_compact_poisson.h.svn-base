#pragma once
#ifndef __FUNCTOR_MLS_LAPLACIAN_COMPACT_POISSON_H__
#define __FUNCTOR_MLS_LAPLACIAN_COMPACT_POISSON_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"
#include "functor_mls_helper_compact_poisson.h"

namespace LAMMPS_NS {

  namespace MLS {

    template<class PairIsph>
    class FunctorOuterLaplacianCompactPoisson : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterLaplacianCompactPoisson(PairIsph *isph, 
                                          double *u, 
                                          double alpha, 
                                          double *f,
                                          double **normal, double *g,
                                          double *material = NULL,
                                          double *u_laplace = NULL)
        :  FunctorOuter<PairIsph>(isph),
           _u_helper(isph, u),
           _p_helper(isph, f, normal, g),
           _alpha(alpha),
           _material(material),
           _u_laplace(u_laplace)
      { }
      
      void operator()(const int ii);
      
      double getScalarLaplacian();

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      FunctorOuterHelper<PairIsph> _u_helper;
      FunctorOuterHelperCompactPoisson<PairIsph> _p_helper;

      double _alpha, *_material;
      double _u_laplace_at_i, *_u_laplace;
    };
    
    template<class PairIsph> inline double
    FunctorOuterLaplacianCompactPoisson<PairIsph>::getScalarLaplacian() {
      return _u_laplace_at_i;
    }

    template<class PairIsph> inline void
    FunctorOuterLaplacianCompactPoisson<PairIsph>::operator()(const int ii) {
      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype); 

      // initialize with zero
      _u_laplace_at_i = 0.0;

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

      const int p[3][3] = { {2,0,0},    // dx
                            {0,2,0},    // dy
                            {0,0,2} };  // dz

      // mass matrix
      double *M = _pair->M[i];
      
      // compute laplacian
      for (int k2=0;k2<_dim;++k2) {
        int idx = 0;
        double dq = 0.0, *qq, *row = NULL;
        
        // d/d{x,y,z}, pick a row correspond to the combination of p
        _pair->basis->dval(np, rth,   p[k2][0], p[k2][1], p[k2][2],   idx, dq); 
        row = &VIEW2(M, ndofpl, idx, 0);
        
        {
          // scalar field laplacian 
          if (_u_helper.hasScalarField()) {
            qq = _u_helper.getScalarHelper();
            _u_laplace_at_i += dq*util.dotVectors(ndofpl, row, ndofpl, qq, 1);
          }
          
          // penalty 
          if (_p_helper.hasScalarField()) {
            qq = _p_helper.getScalarHelper();
            _u_laplace_at_i += dq*util.dotVectors(ndofpl, row, ndofpl, qq, 1);
          }
        }

      }

      if (_material != NULL) {
        _u_laplace_at_i *= _material[i];
        //_u_laplace_at_i += util.dotVectors(_dim, _u_grad.getScalarGradient(), _material_grad.getScalarGradient());
      }

      _u_laplace_at_i *= _alpha;

      // store laplacian at i if _u_laplace is not null
      if (_u_laplace != NULL)
        _u_laplace[i] = _u_laplace_at_i;
    }

  }
}
#endif

