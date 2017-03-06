#pragma once
#ifndef __FUNCTOR_UNCORRECTED_GRADIENT_H__
#define __FUNCTOR_UNCORRECTED_GRADIENT_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {
  
    template<class PairIsph>
    class FunctorOuterGradient : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterGradient(PairIsph *isph, double alpha = 1.0)
        : FunctorOuter<PairIsph>(isph),_mirror(NULL),_f(NULL),_f3(NULL),_alpha(alpha),_grad(NULL) { }
      FunctorOuterGradient(PairIsph *isph, double *f, double alpha = 1.0, 
                           double **grad = NULL)
        : FunctorOuter<PairIsph>(isph),_mirror(NULL),_f(f),_f3(NULL),_alpha(alpha),_grad(grad) { }
      FunctorOuterGradient(PairIsph *isph, double **f3, double alpha = 1.0)
        : FunctorOuter<PairIsph>(isph),_mirror(NULL),_f(NULL),_f3(f3),_alpha(alpha),_grad(NULL) { }
      FunctorOuterGradient(PairIsph *isph, double *f, double **f3, double alpha = 1.0) 
        : FunctorOuter<PairIsph>(isph),_mirror(NULL),_f(f),_f3(f3),_alpha(alpha),_grad(NULL) { }
      virtual~FunctorOuterGradient() { 
        if (_mirror != NULL)
          delete _mirror; 
      }

      void operator()(const int ii);
      double* getScalarGradient();
      double* getVectorGradient(const int k);
      double* getGradientOperator(const int k);

      virtual void createMirrorFunction(); 
    protected:
      MirrorFunction *_mirror; 

      double *_f,**_f3,_alpha,**_grad;
      double _f3_grad_at_i[3][3], _f_grad_at_i[3];
      Epetra_SerialDenseVector _grad_operator[3];

      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline double*
    FunctorOuterGradient<PairIsph>::getScalarGradient() {
      return _f_grad_at_i;
    }
    template<class PairIsph> inline double*
    FunctorOuterGradient<PairIsph>::getVectorGradient(const int k) {
      return _f3_grad_at_i[k];
    }
    template<class PairIsph> inline double*
    FunctorOuterGradient<PairIsph>::getGradientOperator(const int k) {
      return _grad_operator[k].Values();
    }
    template<class PairIsph> inline void
    FunctorOuterGradient<PairIsph>::createMirrorFunction() {
      _mirror = new MirrorNothing(_dim, _type, _pair->h);
    }

    template<class PairIsph> inline void
    FunctorOuterGradient<PairIsph>::operator()(const int ii) {
      if (_mirror == NULL)
        createMirrorFunction();

      int i = _ilist[ii];
      int itype = _type[i];

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // gradient is computed on the given storage or member variable
      memset(_f3_grad_at_i, 0, sizeof(double)*9);
      memset(_f_grad_at_i, 0, sizeof(double)*3);

      for(int k=0; k< _dim; ++k)
        _grad_operator[k].Size(jnum+1);

      if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
        return;

      // mirroring for particle i
      _mirror->setFreeParticle(i);

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

        double cutsq = LOOKUP(_pair->cutsq, itype, jtype); 
        if (rsq < cutsq) {
          _mirror->setMirrorParticle(j);
          double coeff = 1.0;
          if(_pair->getParticleKind(jtype) == PairIsph::Solid)
            coeff = _mirror->computeMirrorCoefficient(sqrt(cutsq));

          double r = sqrt(rsq) + ISPH_EPSILON;
          double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));
        
          // 1st order consistent corrected gradient
          // \grad_{1} psi_{ij}^{k} = (psi_{j} - psi_{i}) G_i^{kq} dwdr (r_{ij}^{q} / r) V_{j}
          //                        = (psi_{j} - psi_{i}) ( G_i^{kq} r_{ij}^{q} ) dwdr / r V_{j} 
          for (int k2=0;k2<_dim;++k2) {   // k

            double vjtmp = dwdr/r*_vfrac[j]*coeff;          
            if (_f != NULL) {
              _f_grad_at_i[k2] += rij[k2]*(_f[j] - _f[i])*vjtmp;
            } 
            if ( _f3 != NULL) {
              for (int k1=0;k1<_dim;++k1)
                _f3_grad_at_i[k1][k2] += rij[k2]*(_f3[j][k1] - _f3[i][k1])*vjtmp;
            }
            _grad_operator[k2][jj+1] = rij[k2]*vjtmp*_alpha;
            _grad_operator[k2][0] -= rij[k2]*vjtmp*_alpha;
          }
        }
      }
    
      // scale with alpha
      if (_f != NULL)
        for (int k=0;k<_dim;++k)
          _f_grad_at_i[k] *= _alpha;

      if (_f3 != NULL)
        for (int k2=0;k2<_dim;++k2)
          for (int k1=0;k1<_dim;++k1)
            _f3_grad_at_i[k2][k1] *= _alpha;
          
      // store gradient at i if grad is not null
      if (_grad != NULL && _f != NULL) 
        memcpy(_grad[i], _f_grad_at_i, sizeof(double)*_dim);
    }

  }
}

#endif

