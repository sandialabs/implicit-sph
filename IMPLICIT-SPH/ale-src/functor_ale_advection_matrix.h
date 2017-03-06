#pragma once
#ifndef __FUNCTOR_ALE_ADVECTION_MATRIX_H__
#define __FUNCTOR_ALE_ADVECTION_MATRIX_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class GradientOperator>
  class FunctorOuterAleAdvectionMatrix: public FunctorOuter<PairIsph> {
  public:
    FunctorOuterAleAdvectionMatrix(PairIsph *isph, 
                                   double **v, double **xdot, 
                                   double alpha = 1.0) 
      : FunctorOuter<PairIsph>(isph), 
        _op_grad(isph, alpha),
        _v(v), 
        _xdot(xdot) { }
    
    void enterFor();
    void operator()(const int ii);

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;

    GradientOperator<PairIsph> _op_grad;
    double **_v, **_xdot;

    Epetra_IntSerialDenseVector _idx;
    Epetra_SerialDenseVector _val;  
  };
  
  template<class PairIsph,
           template<typename> class GradientOperator>
  inline void
  FunctorOuterAleAdvectionMatrix<PairIsph,GradientOperator>::enterFor() {
    int n = _pair->tags_in_cut->MaxNumIndices();
    _idx.Size(n);
    _val.Size(n);

    // use the same filter
    _op_grad.setFilter(_filter);
  }
  template<class PairIsph,
           template<typename> class GradientOperator>
  inline void 
  FunctorOuterAleAdvectionMatrix<PairIsph,GradientOperator>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];
    
    if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
      return; 
    
    // compute gradients
    _op_grad(ii);
    
    // extract pre-allocated pointers
    double *val = _val.Values();
    int *idx = _idx.Values();

    // this already includes the self contribution
    int jnum = _op_grad.getSize();
    for (int jj=0;jj<jnum;++jj) {
      int j = _op_grad.getNeighbor(jj);

      idx[jj] = _tag[j];
      val[jj] = 0.0;

      for (int k=0;k<_dim;++k)
        val[jj] += (_v[i][k] - _xdot[i][k])*_op_grad.getValue(jj, k);
    }
    
    // add the contribution to the matrix
    _pair->A.crs->SumIntoGlobalValues(_tag[i], jnum, val, idx);
  }
}
#endif

