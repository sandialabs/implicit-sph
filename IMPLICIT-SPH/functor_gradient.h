#pragma once
#ifndef __FUNCTOR_GRADIENT_H__
#define __FUNCTOR_GRADIENT_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {
  
    template<class PairIsph, bool AntiSymmetric = false>
    class FunctorOuterGradient : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterGradient(PairIsph *isph, 
                           double *f, 
                           double alpha = 1.0, 
                           double **grad = NULL)
        : FunctorOuter<PairIsph>(isph),
          _mirror(NULL),
          _f(f),
          _f3(NULL),
          _alpha(alpha),
          _grad(grad) { }
      FunctorOuterGradient(PairIsph *isph, 
                           double **f3, 
                           double alpha = 1.0)
        : FunctorOuter<PairIsph>(isph),
          _mirror(NULL),
          _f(NULL),
          _f3(f3),
          _alpha(alpha),
          _grad(NULL) { }
      FunctorOuterGradient(PairIsph *isph, 
                           double *f, 
                           double **f3, 
                           double alpha = 1.0) 
        : FunctorOuter<PairIsph>(isph),
          _mirror(NULL),
          _f(f),
          _f3(f3),
          _alpha(alpha),
          _grad(NULL) { }

      virtual~FunctorOuterGradient() { 
        if (_mirror != NULL)
          delete _mirror; 
      }

      void operator()(const int ii);
      double* getScalarGradient();
      double* getVectorGradient(const int k = 0);
      
      virtual void createMirrorFunction(); 
    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      MirrorFunction *_mirror; 

      double *_f,**_f3,_alpha,**_grad;
      double _f3_grad_at_i[3][3], _f_grad_at_i[3];
    };

    template<class PairIsph, bool AntiSymmetric> inline double*
    FunctorOuterGradient<PairIsph,AntiSymmetric>::getScalarGradient() {
      return _f_grad_at_i;
    }
    template<class PairIsph, bool AntiSymmetric> inline double*
    FunctorOuterGradient<PairIsph,AntiSymmetric>::getVectorGradient(const int k) {
      return _f3_grad_at_i[k];
    }
    template<class PairIsph, bool AntiSymmetric> inline void
    FunctorOuterGradient<PairIsph,AntiSymmetric>::createMirrorFunction() {
      _mirror = new MirrorNothing(_dim, _type, _pair->h);
    }

    template<class PairIsph, bool AntiSymmetric> inline void
    FunctorOuterGradient<PairIsph,AntiSymmetric>::operator()(const int ii) {
      if (_mirror == NULL)
        createMirrorFunction();

      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype);

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // gradient is computed on the given storage or member variable
      if (_f3 != NULL)
        memset(&_f3_grad_at_i[0][0], 0, sizeof(double)*9);

      if (_f != NULL)
        memset(&_f_grad_at_i[0], 0, sizeof(double)*3);

      if (_filter != NULL && !_filter->yes(ikind))
        return;

      // mirroring for particle i
      _mirror->setFreeParticle(i); 

      // initialize the correction matrix
      const double *G = (AntiSymmetric ? _pair->Gi : _pair->Gc[i]);

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

        const double cutsq = LOOKUP(_pair->cutsq, itype, jtype); 
        if (rsq < cutsq) {
          _mirror->setMirrorParticle(j);
          double coeff = 1.0;
          if (!(ikind & PairIsph::Solid) && (jkind & PairIsph::Solid))
            coeff = _mirror->computeMirrorCoefficient(sqrt(cutsq));

          const double r = sqrt(rsq) + ISPH_EPSILON;
          const double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));
        
          // 1st order consistent corrected gradient
          // \grad_{1} psi_{ij}^{k} = (psi_{j} - psi_{i}) G_i^{kq} dwdr (r_{ij}^{q} / r) V_{j}
          //                        = (psi_{j} - psi_{i}) ( G_i^{kq} r_{ij}^{q} ) dwdr / r V_{j} 
          const double vfrac = (AntiSymmetric ? sqrt(_vfrac[i]*_vfrac[j]) : _vfrac[j]);
          const double vjtmp = dwdr/r*vfrac*coeff;          
          for (int k2=0;k2<_dim;++k2) {   // k
            double gitmp = 0.0;
            for (int k1=0;k1<_dim;++k1)   // q
              gitmp += VIEW2(G, _dim, k1, k2)*rij[k1];
            
            const double ijtmp = gitmp*vjtmp;

            if (_f != NULL) 
              _f_grad_at_i[k2] += ijtmp*(sphOperator<AntiSymmetric>(_f[i], _f[j]));

            if ( _f3 != NULL) 
              for (int k1=0;k1<_dim;++k1)
                _f3_grad_at_i[k1][k2] += ijtmp*(sphOperator<AntiSymmetric>(_f3[i][k1], _f3[j][k1])); 
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

    template<class PairIsph> using FunctorOuterGradientSymmetric     = FunctorOuterGradient<PairIsph,false>;
    template<class PairIsph> using FunctorOuterGradientAntiSymmetric = FunctorOuterGradient<PairIsph,true>;
  }
}

#endif

