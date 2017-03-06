#pragma once
#ifndef __FUNCTOR_PHASE_DIVERGENCE_H__
#define __FUNCTOR_PHASE_DIVERGENCE_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {

    template<class PairIsph>
    class FunctorOuterPhaseDivergence : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterPhaseDivergence(PairIsph *isph, double **normal, double *mag,  
                                  double *div = NULL)
        : FunctorOuter<PairIsph>(isph),
          _normal(normal), _mag(mag), 
          _div(div) { }

      double getDivergence();

      void operator()(const int ii);      

    protected:
      double **_normal, *_mag, *_div;
      double _div_at_i;
      
      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline double
    FunctorOuterPhaseDivergence<PairIsph>::getDivergence() {
      return _div_at_i;
    }

    template<class PairIsph> inline void
    FunctorOuterPhaseDivergence<PairIsph>::operator()(const int ii) {
      const int i = _ilist[ii];
      const int itype = _type[i];

      const int ikind  = _pair->getParticleKind(itype);
      const int iphase = _pair->getParticlePhase(itype);      

      // neighbor around i
      const int *jlist = _firstneigh[i];
      const int jnum = _numneigh[i];

      // divergence is computed on the given storage or member variable
      _div_at_i = 0.0;

      if (_filter->yes(ikind) && _mag[i] > ISPH_EPSILON) {

        // initialize the correction matrix
        double *G = _pair->Gc[i];
        
        for (int jj=0;jj<jnum;++jj) {
          const int j = (jlist[jj] & NEIGHMASK);
          const int jtype = _type[j];
          
          const int jkind  = _pair->getParticleKind(jtype);
          const int jphase = _pair->getParticlePhase(jtype);
          
          if (_filter->yes(ikind, jkind)) {
            if (_pair->st.csf.color->yes(iphase, jphase)) {
              if (_mag[j] > ISPH_EPSILON) {
                // compute distance between partile i and j
                double rsq = 0.0, rij[3] = {};
                for (int k=0;k<_dim;++k) {
                  rij[k] = _x[i][k] - _x[j][k];
                  rsq += (rij[k]*rij[k]);
                }
                
                if (rsq < LOOKUP(_pair->cutsq, itype, jtype)) {
                  const double r = sqrt(rsq) + ISPH_EPSILON;
                  const double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));
                  
                  // 1st order consistent corrected divergence
                  for (int k2=0;k2<_dim;++k2) {   // k
                    double gitmp = 0.0;
                    for (int k1=0;k1<_dim;++k1)   // q
                      gitmp += VIEW2(G, _dim, k1, k2)*rij[k1];
                    
                    double sign = _pair->st.csf.color->sign(iphase, jphase);
                    _div_at_i += gitmp*(sign*_normal[j][k2] - _normal[i][k2])*dwdr/r*_vfrac[j]; 
                  }
                }
              }
            }
          }
        }
      }
      if (_div != NULL) 
        _div[i] = _div_at_i;
    }
    
  }
}

#endif

