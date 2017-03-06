#pragma once
#ifndef __FUNCTOR_PAIRWISE_FORCE_H__
#define __FUNCTOR_PAIRWISE_FORCE_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {
  
  template<class PairIsph>
  class FunctorOuterPairwiseForce : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterPairwiseForce(PairIsph *isph, double **f)
      : FunctorOuter<PairIsph>(isph),
        _f(f) { }

    void operator()(const int ii);
    
  protected:
    double **_f, _f_at_i[3];
  public:
    double _f_sum[3];

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  template<class PairIsph>
  inline void
  FunctorOuterPairwiseForce<PairIsph>::operator()(const int ii) {
    const int i = _ilist[ii];
    const int itype = _type[i];

    const int ikind  = _pair->getParticleKind(itype);
    const int iphase = _pair->getParticlePhase(itype);

    // neighbor around i
    const int *jlist = _firstneigh[i];
    const int jnum = _numneigh[i];

    // _f_at_i
    memset(_f_at_i, 0, sizeof(double)*_dim);
    
    if (_filter != NULL && _filter->yes(ikind)) {
      for (int jj=0;jj<jnum;++jj) {
        const int j = (jlist[jj] & NEIGHMASK);
        const int jtype = _type[j];
        
        const int jkind  = _pair->getParticleKind(jtype);
        const int jphase = _pair->getParticlePhase(jtype);
        
        if (_filter != NULL && _filter->yes(ikind, jkind)) {
          // compute distance between partile i and j
          double rsq = 0.0, rij[3] = { };
          for (int k=0;k<_dim;++k) {
            rij[k] = _x[i][k] - _x[j][k];
            rsq += (rij[k] * rij[k]);
          }
          
          const double cutsq = LOOKUP(_pair->cutsq, itype, jtype);
          if (rsq < cutsq) { 
            const double s = LOOKUP(_pair->st.pf.s, iphase, jphase);
            const double r = sqrt(rsq) + ISPH_EPSILON;
            const double c = sqrt(cutsq);

            const double f = _pair->st.pf.force->val(s, r, c);
            
            for (int k=0;k<_dim;++k)
              _f_at_i[k] += -f*rij[k]/r;
          }
        }
      }

      // accumulate force
      for (int k=0;k<_dim;++k) {
        _f[i][k] += _f_at_i[k];
        _f_sum[k] += _f_at_i[k];
      }
    }

    
  }

}

#endif

