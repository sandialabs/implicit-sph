#pragma once
#ifndef __FUNCTOR_ADVANCE_TIME_END_H__
#define __FUNCTOR_ADVANCE_TIME_END_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph>
  class FunctorOuterAdvanceTimeEnd : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterAdvanceTimeEnd(PairIsph *isph, 
                               double dt, 
                               double **v, double *p,
                               int extend = 0) 
      : FunctorOuter<PairIsph>(isph),
        _dt(dt),
        _v(v),_vnp1(isph->vstar),
        _p(p),_dp(isph->dp),
        _extend(extend) { }
    void enterFor();
    void operator()(const int ii);
    size_t getNumberOfWork() const;

  protected:
    double _dt,**_v,**_vnp1,*_p,*_dp;
    int _extend;

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  template<class PairIsph> inline size_t
  FunctorOuterAdvanceTimeEnd<PairIsph>::getNumberOfWork() const { 
    return (FunctorOuter<PairIsph>::getNumberOfWork() + _extend);
  }    

  template<class PairIsph> inline void
  FunctorOuterAdvanceTimeEnd<PairIsph>::enterFor() {
    // if other functors are properly established, no need for specialization
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorAdvanceTimeEnd:: Filter is not allowed");
  }

  template<class PairIsph> inline void
  FunctorOuterAdvanceTimeEnd<PairIsph>::operator()(const int ii) {
    int i = ii; //_ilist[ii]; 
    int itype = _type[i];
    
    if (_pair->isParticleFixed(itype)) {
      for (int k=0;k<_dim;++k) 
	_v[i][k] = _vnp1[i][k];
      return;
    }

    _p[i] += _dp[i];    
    for (int k=0;k<_dim;++k) {
      double delta = 0.5*_dt*(_vnp1[i][k] + _v[i][k]);  
      
      _x[i][k] += delta;
      _v[i][k]  = _vnp1[i][k];
    }
  }

}

#endif

