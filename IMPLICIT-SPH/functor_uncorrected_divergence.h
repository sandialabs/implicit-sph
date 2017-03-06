#pragma once
#ifndef __FUNCTOR_UNCORRECTED_DIVERGENCE_H__
#define __FUNCTOR_UNCORRECTED_DIVERGENCE_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {

    template<class PairIsph>
    class FunctorOuterDivergence : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterDivergence(PairIsph *isph, double **f, double alpha, 
                             double *div)
        : FunctorOuter<PairIsph>(isph),_f(f),_alpha(alpha),_div(div) { }
      void operator()(const int ii);
      double getDivergence();
    protected:
      double **_f, _alpha, *_div;
      double _div_at_i;

      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline double
    FunctorOuterDivergence<PairIsph>::getDivergence() {
      return _div_at_i;
    }
    // Reference version
    template<class PairIsph> inline void
    FunctorOuterDivergence<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // divergence is computed on the given storage or member variable
      _div_at_i = 0.0;

      if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
        return;

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
      
        if (rsq < LOOKUP(_pair->cutsq, itype, jtype)) {
          double r = sqrt(rsq) + ISPH_EPSILON;
          double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));
        
          // 1st order consistent divergence
          for (int k=0;k<_dim;++k) {   // k
            _div_at_i += rij[k]*(_f[j][k] - _f[i][k])*dwdr/r*_vfrac[j];
          }
        }
      }
    
      // scale with alpha
      _div_at_i *= _alpha;

      if (_div != NULL) 
        _div[i] = _div_at_i;
    }

  }
}

#endif

