#pragma once
#ifndef __FUNCTOR_COMPUTE_SHIFT_H__
#define __FUNCTOR_COMPUTE_SHIFT_H__


#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {
  
  template<class PairIsph>
  class FunctorOuterComputeShift : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterComputeShift(PairIsph *isph, 
                             double alpha,
                             double shiftcut,
                             double nonfluidweight,
                             double **dr)
      : FunctorOuter<PairIsph>(isph),
        _alpha(alpha),
        _shiftcutsq(shiftcut*shiftcut),
        _nonfluidweight(nonfluidweight),
        _dr(dr){ }

    void enterFor();
    void operator()(const int ii);

  protected:
    double _alpha,_shiftcutsq,_nonfluidweight,**_dr;
    FilterBinary _filter_op;

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  template<class PairIsph> inline void
  FunctorOuterComputeShift<PairIsph>::enterFor() {
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorComputeShift:: Filter is not allowed - internal filter behave same as one in FunctorApplyShift");
    
    _filter_op.setPairYes(PairIsph::Fluid, PairIsph::All); 
  }

  template<class PairIsph> inline void
  FunctorOuterComputeShift<PairIsph>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];
    int ikind = _pair->getParticleKind(itype);

    // neighbor around i
    int *jlist = _firstneigh[i];
    int jnum = _numneigh[i];

    memset(&_dr[i][0], 0, sizeof(double)*_dim);

    // this should be applied to Fluid particles only
    if (!_filter_op.yes(ikind)) 
      return;

    // initialize ri
    int cnt = 0;
    double ri = 0.0;
    for (int jj=0;jj<jnum;++jj) {
      // j th particle; always put neighmask two bits are reserved for LAMMPS
      int j = (jlist[jj] & NEIGHMASK);
      int jtype = _type[j];
      int jkind = _pair->getParticleKind(jtype);

      if (!_filter_op.yes(ikind, jkind))
        continue;
      
      // compute distance between partile i and j
      double rsq = 0.0;
      for (int k=0;k<_dim;++k) 
        rsq += pow(_x[i][k] - _x[j][k], 2);

      // integrate
      double rth = min(LOOKUP(_pair->cutsq, itype, jtype), _shiftcutsq);
      if (rsq < rth) {
        ++cnt;
        ri += sqrt(rsq);
      }
    }

    if (cnt)
      ri /= (double)cnt;
    
    for (int jj=0;jj<jnum;++jj) {
      // j th particle; always put neighmask two bits are reserved for LAMMPS
      int j = (jlist[jj] & NEIGHMASK);
      int jtype = _type[j];
      int jkind = _pair->getParticleKind(jtype);

      if (!_filter_op.yes(ikind, jkind))
        continue;
      
      // compute distance between partile i and j
      double rsq = 0.0, rij[3] = {};
      for (int k=0;k<_dim;++k) {
        rij[k] = _x[i][k] - _x[j][k];
        rsq += rij[k]*rij[k];
      }

      // integrate
      double rth = min(LOOKUP(_pair->cutsq, itype, jtype), _shiftcutsq);
      if (rsq < rth) {
        double r = sqrt(rsq) + ISPH_EPSILON, rir2 = pow(ri/r,2);
        double beta = _alpha/r*rir2*(1.0 + (!(jkind & PairIsph::Fluid))*_nonfluidweight*rir2);
	
        // alpha = C*dt*Umax
        for (int k=0;k<_dim;++k) 
          _dr[i][k] += beta*rij[k];
      }
    }
  }
  
}

#endif

