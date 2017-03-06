#pragma once
#ifndef __FUNCTOR_GRADIENT_CORRECTION_H__
#define __FUNCTOR_GRADIENT_CORRECTION_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {

    template<class PairIsph>
    class FunctorOuterGradientCorrection : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterGradientCorrection(PairIsph *isph) : FunctorOuter<PairIsph>(isph) { }
      void operator()(const int ii);
    protected:
      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline void 
    FunctorOuterGradientCorrection<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // initialize the correction matrix
      memset(_pair->Gc[i], 0, sizeof(double)*_dimsq);

      if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
        return;
    
      // local correction matrix
      double G[9] = {};

      for (int jj=0;jj<jnum;++jj) {
        int j = (jlist[jj] & NEIGHMASK);
        int jtype = _type[j];

        if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype),
                                             _pair->getParticleKind(jtype)))
          continue;

        // compute distance between partile i and j
        double rsq = 0.0, rij[3];
        for (int k=0;k<_dim;++k) {
          rij[k] = _x[i][k] - _x[j][k];
          rsq += (rij[k]*rij[k]);
        }

        if (rsq < LOOKUP(_pair->cutsq, itype, jtype)) {
          double r = sqrt(rsq) + ISPH_EPSILON;
          double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));
        
          // 3x3 rank 1 update should be vectorized
          // (G^{mn})^{-1} = r^{m} dwdr r^{n} / r V
          //               = r^{m} r^{n} dwdr/r V ... symmetry
          for (int k2=0;k2<_dim;++k2)   // n
            for (int k1=0;k1<_dim;++k1) // m
              VIEW2(G, _dim, k1, k2) -= rij[k1]*rij[k2]*dwdr/r*_vfrac[j];
        }
      }
      // invert G and store it to Gc
      double ipiv[3];
      util.invertDenseMatrix(_dim, G, _pair->Gc[i], 3, &ipiv[0]);
    }

  }
}

#endif

