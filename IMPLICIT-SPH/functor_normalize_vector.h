#pragma once
#ifndef __FUNCTOR_NORMALIZE_VECTOR_H__
#define __FUNCTOR_NORMALIZE_VECTOR_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph>
  class FunctorOuterNormalizeVector : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterNormalizeVector(PairIsph *isph, double **vec, double *mag, int extend = 0)
      : FunctorOuter<PairIsph>(isph),
      _vec(vec),_mag(mag),_extend(extend) { }
    void operator()(const int ii);
    size_t getNumberOfWork() const;

  protected:
    double **_vec,*_mag;
    int _extend;

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  template<class PairIsph> inline void
  FunctorOuterNormalizeVector<PairIsph>::operator()(const int ii) {
    const int i = _ilist[ii];
    
    double *vec = _vec[i]; 
    const double mag = sqrt(util.dotVectors(_dim, vec, vec));

    if (mag != 0.0)
      for (int k=0;k<_dim;++k)
        vec[k] /= mag;
    
    if (_mag != NULL)
      _mag[i] = mag;
  }
  
  template<class PairIsph> inline size_t 
  FunctorOuterNormalizeVector<PairIsph>::getNumberOfWork() const {
    return (FunctorOuter<PairIsph>::getNumberOfWork() + _extend);
  }

}

#endif

