#pragma once
#ifndef __FUNCTOR_CORRECT_PRESSURE_H__
#define __FUNCTOR_CORRECT_PRESSURE_H__


#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph>
  class FunctorOuterCorrectPressure : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterCorrectPressure(PairIsph *isph, double *p, double *dp, int extend = 0)
      : FunctorOuter<PairIsph>(isph),_p(p),_dp(dp),_extend(extend) { }
    void operator()(const int ii);
    size_t getNumberOfWork() const;

  protected:

    double *_p,*_dp;
    int _extend;

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  template<class PairIsph> inline void
  FunctorOuterCorrectPressure<PairIsph>::operator()(const int ii) {
    // update pressure for all types of particles.
    // in the poisson problem, dp is set to be either initial 
    // pressure (for timestep =0) or zero pressure difference.
    // so this correction either update pressure to initial condition
    // or zero for solid particles

    if (_pair->ns.is_incremental_pressure_used) {
      _p[ii] += _dp[ii];
    } else {
      _p[ii] = _dp[ii];
    }
      
  }
  template<class PairIsph> inline size_t 
  FunctorOuterCorrectPressure<PairIsph>::getNumberOfWork() const {
    return (FunctorOuter<PairIsph>::getNumberOfWork() + _extend);
  }

}

#endif

