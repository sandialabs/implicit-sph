#pragma once
#ifndef __FUNCTOR_MLS_LAPLACIAN_H__
#define __FUNCTOR_MLS_LAPLACIAN_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"
#include "functor_mls_helper.h"

namespace LAMMPS_NS {

  namespace MLS {

    template<class PairIsph>
    class FunctorOuterLaplacian : public FunctorOuterHelper<PairIsph> {
    public:
      FunctorOuterLaplacian(PairIsph *isph, 
                            double *f,  
                            double alpha, double* scal = NULL)
        : FunctorOuterHelper<PairIsph>(isph, f),
          _alpha(alpha),
          _scal(scal),
          _grad_f(NULL),
          _grad_scal(NULL)
          { }

      void enterFor();

      void operator()(const int ii);

      void exitFor();

      double getScalarLaplacian();

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      double _alpha;
      double* _scal;
      double _laplace_at_i;
      FunctorOuterGradient<PairIsph> *_grad_f, *_grad_scal;
    };

    template<class PairIsph> void
    FunctorOuterLaplacian<PairIsph>::enterFor() {
      FunctorOuterHelper<PairIsph>::enterFor();
      if(_scal != NULL) {
        _grad_f = new FunctorOuterGradient<PairIsph>(_pair, this->_f, 1.0);
        _grad_scal = new FunctorOuterGradient<PairIsph>(_pair, _scal, 1.0);
        _grad_f->setFilter(_filter);
        _grad_scal->setFilter(_filter);
      }
    }


    template<class PairIsph> inline double
    FunctorOuterLaplacian<PairIsph>::getScalarLaplacian() {
      return _laplace_at_i;
    }

    template<class PairIsph> inline void
    FunctorOuterLaplacian<PairIsph>::operator()(const int ii) {

      int i = _ilist[ii];
      int itype = _type[i];

      // initialize
      _laplace_at_i = 0.0;

      if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
        return;
      
      // compute helper; q_at_i = sum P(x_j) W_ij f_j
      FunctorOuterHelper<PairIsph>::operator()(ii); 

      // basis arguments
      int np   = _pair->np;
      int ndof = _pair->basis->ndof(np);
      double rth = LOOKUP(_pair->h, itype, itype);

      const int p[3][3] = { {2,0,0},    // dx2
                            {0,2,0},    // dy2
                            {0,0,2} };  // dz2
      
      // mass matrix
      double *M = _pair->M[i];

      // compute laplacian
      for (int k=0;k<_dim;++k) {
        int idx = 0;
        double dq = 0.0, *qq = NULL, *row = NULL;

        // d/d{x,y,z}, pick a row correspond to the combination of p
        _pair->basis->dval(np, rth,   p[k][0], p[k][1], p[k][2],   idx, dq); 
        row = &VIEW2(M, ndof, idx, 0);
        
        // laplacian
        qq = FunctorOuterHelper<PairIsph>::getScalarHelper();
        _laplace_at_i += dq*util.dotVectors(ndof, row, ndof, qq, 1);
      }

      //add grad(scal) . grad(f)
      if(_scal != NULL) {
        (*_grad_f)(ii);
        (*_grad_scal)(ii);
        _laplace_at_i *= _scal[i];
        _laplace_at_i += util.dotVectors(_dim, _grad_f->getScalarGradient(), _grad_scal->getScalarGradient());
      }

      _laplace_at_i *= _alpha;
    }

    template<class PairIsph> void
    FunctorOuterLaplacian<PairIsph>::exitFor() {
      if(_scal != NULL) {
        delete _grad_f;
        delete _grad_scal;
      }
      FunctorOuterHelper<PairIsph>::exitFor();
    }

  }
}
#endif

