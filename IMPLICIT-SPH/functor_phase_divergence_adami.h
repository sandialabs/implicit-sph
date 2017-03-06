#pragma once
#ifndef __FUNCTOR_PHASE_DIVERGENCE_ADAMI_H__
#define __FUNCTOR_PHASE_DIVERGENCE_ADAMI_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"
#include "color.h"

namespace LAMMPS_NS {

  namespace Corrected {

    template<class PairIsph>
    class FunctorOuterPhaseDivergenceAdami : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterPhaseDivergenceAdami(PairIsph *isph,
                                       double **normal, double *mag,  
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
    FunctorOuterPhaseDivergenceAdami<PairIsph>::getDivergence() {
      return _div_at_i;
    }

    template<class PairIsph> inline void
    FunctorOuterPhaseDivergenceAdami<PairIsph>::operator()(const int ii) {
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

        double denominator = 0.0;
        double numerator = 0.0;
        
        for (int jj=0;jj<jnum;++jj) {
          const int j = (jlist[jj] & NEIGHMASK);
          const int jtype = _type[j];
          
          const int jkind  = _pair->getParticleKind(jtype);
          const int jphase = _pair->getParticlePhase(jtype);
          
          if (_filter->yes(ikind, jkind)) {
            if (_pair->st.csf.color->yes(iphase, jphase)) {
              if (_mag[j] > ISPH_EPSILON) {
                double rsq = 0.0, rij[3] = {};
                for (int k=0;k<_dim;++k) {
                  rij[k] = _x[i][k] - _x[j][k];
                  rsq += (rij[k]*rij[k]);
                }
                
                if (rsq < LOOKUP(_pair->cutsq, itype, jtype)) {
                  const double r = sqrt(rsq) + ISPH_EPSILON;
                  const double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));
                  
                  double nij[3], sign = _pair->st.csf.color->sign(iphase, jphase);
                  for (int k=0;k<_dim;++k)
                    nij[k] = (_normal[i][k] - sign*_normal[j][k]);
                  
                  const double vjtmp = dwdr*_vfrac[j];
                  denominator += r*vjtmp;
                  numerator   += util.dotVectors(_dim, nij, rij)/r*vjtmp;
                }
              }
            }
          }
        }
        
        if (denominator != 0.0)
          _div_at_i = _dim*numerator/denominator;
      }

      if (_div != NULL) 
        _div[i] = _div_at_i;
    }
    
  }
}

#endif

