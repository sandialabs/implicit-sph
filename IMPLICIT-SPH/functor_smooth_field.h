#pragma once
#ifndef __FUNCTOR_SMOOTH_FIELD_H__
#define __FUNCTOR_SMOOTH_FIELD_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {
  
    template<class PairIsph>
    class FunctorOuterSmoothField : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterSmoothField(PairIsph *isph, 
                              double **f3, double **sf3)
        : FunctorOuter<PairIsph>(isph),
          _f(NULL),
          _sf(NULL),
          _f3(f3),
          _sf3(sf3) { }

      FunctorOuterSmoothField(PairIsph *isph, 
                              double *f, double *sf,
                              double **f3 = NULL, double **sf3 = NULL)
        : FunctorOuter<PairIsph>(isph),
          _f(f),
          _sf(sf),
          _f3(f3),
          _sf3(sf3) { }

      void operator()(const int ii);

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      double *_f,*_sf,**_f3,**_sf3;
    };

    template<class PairIsph> inline void
    FunctorOuterSmoothField<PairIsph>::operator()(const int ii) {
      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype);

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      if (_filter != NULL && !_filter->yes(ikind))
        return;

      double sf_at_i, sf3_at_i[3];

      // self contribution
      {
        const double w = _pair->kernel->val(0.0, LOOKUP(_pair->h, itype, itype));
        const double tmp = w*_vfrac[i];
        
        if (_f != NULL)
          sf_at_i = _f[i]*tmp;
        
        if (_f3 != NULL)
          for (int k1=0;k1<_dim;++k1)
            sf3_at_i[k1] = _f3[i][k1]*tmp;
      }

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
          const double r = sqrt(rsq) + ISPH_EPSILON;
          const double w = _pair->kernel->val(r, LOOKUP(_pair->h, itype, jtype));
          const double tmp = w*_vfrac[j];
          
          if (_f != NULL)
            sf_at_i += _f[j]*tmp;

          if (_f3 != NULL)
            for (int k1=0;k1<_dim;++k1)
              sf3_at_i[k1] += _f3[j][k1]*tmp;
        }
      }
      
      if (_sf != NULL) 
        _sf[i] = sf_at_i;
      
      if (_sf3 != NULL)
        memcpy(_sf3[i], sf3_at_i, sizeof(double)*_dim);        
    }

  }
}

#endif

