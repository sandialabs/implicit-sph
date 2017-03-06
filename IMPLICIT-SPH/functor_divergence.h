#pragma once
#ifndef __FUNCTOR_DIVERGENCE_H__
#define __FUNCTOR_DIVERGENCE_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {


    template<class PairIsph, bool AntiSymmetric>
    class FunctorOuterDivergence : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterDivergence(PairIsph *isph, double **f, 
                             double alpha = 1.0, 
                             double *div = NULL)
        : FunctorOuter<PairIsph>(isph),
          _mirror(NULL),
          _f(f),
          _alpha(alpha),
          _div(div) { }
      virtual~FunctorOuterDivergence() {
        if (_mirror != NULL)
          delete _mirror;
      }

      void operator()(const int ii);
      double getDivergence();

      virtual void createMirrorFunction();
    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      MirrorFunction *_mirror;

      double **_f, _alpha, *_div;
      double _div_at_i;
    };

    template<class PairIsph, bool AntiSymmetric> inline double
    FunctorOuterDivergence<PairIsph,AntiSymmetric>::getDivergence() {
      return _div_at_i;
    }
    template<class PairIsph, bool AntiSymmetric> inline void
    FunctorOuterDivergence<PairIsph,AntiSymmetric>::createMirrorFunction() {
      _mirror = new MirrorNothing(_dim, _type, _pair->h);
    }

    template<class PairIsph, bool AntiSymmetric> inline void
    FunctorOuterDivergence<PairIsph,AntiSymmetric>::operator()(const int ii) {
      if (_mirror == NULL)
        createMirrorFunction();

      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype);

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // divergence is computed on the given storage or member variable
      _div_at_i = 0.0;

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
        
          // 1st order consistent corrected divergence
          // \grad_{1} psi_{ij}^{k} = (psi_{j} - psi_{i}) G_i^{kq} dwdr (r_{ij}^{q} / r) V_{j}
          //                        = (psi_{j} - psi_{i}) ( G_i^{kq} r_{ij}^{q} ) dwdr / r V_{j} 
          const double vfrac = (AntiSymmetric ? sqrt(_vfrac[i]*_vfrac[j]) : _vfrac[j]);
          const double vjtmp = dwdr/r*vfrac*coeff;
          for (int k2=0;k2<_dim;++k2) {   // k
            double gitmp = 0.0;
            for (int k1=0;k1<_dim;++k1)   // q
              gitmp += VIEW2(G, _dim, k1, k2)*rij[k1];

            _div_at_i += gitmp*(sphOperator<AntiSymmetric>(_f[i][k2], _f[j][k2]))*vjtmp; 
          }
        }
      }
    
      // scale with alpha
      _div_at_i *= _alpha;

      if (_div != NULL) 
        _div[i] = _div_at_i;
    }

    template<class PairIsph> using FunctorOuterDivergenceSymmetric     = FunctorOuterDivergence<PairIsph,false>;
    template<class PairIsph> using FunctorOuterDivergenceAntiSymmetric = FunctorOuterDivergence<PairIsph,true>;
  }
}

#endif

