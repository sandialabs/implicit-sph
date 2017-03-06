#pragma once
#ifndef __KERNEL_WENDLAND_H__
#define __KERNEL_WENDLAND_H__

#include <iostream>
#include "math.h"
#include "float.h"
#include "kernel.h"

// Wendland quintic kernel

namespace LAMMPS_NS {

  class KernelFuncWendland : public KernelFunction {
  protected:
    void setSmoothingLength(double h);  
  public:
    KernelFuncWendland(unsigned int dim) : KernelFunction(dim) { }
    KernelFuncWendland(unsigned int dim, double h) : KernelFunction(dim) { 
      setSmoothingLength(h); 
    }
    double val (double r) const;
    double dval(double r) const;
  };

  using namespace std;

  inline void
  KernelFuncWendland::setSmoothingLength(double h) {
    // if the given h is not same as cached _h, then _C is recomputed.
    if (fabs(_h - h) > DBL_EPSILON || _h == 0.0) {
      // smoothing parameter is cached
      _h = h;

      // set normalization parameter
      if (_dim == 3) 
        _C = 21.0/(16*M_PI*pow(_h,3));
      else 
        _C =  7.0/(4*M_PI*pow(_h,2));
    }
  }
  
  inline double
  KernelFuncWendland::val(double r) const {
    double s = fabs(r/_h);
    const int n = 4;
    double r_val = pow(1-0.5 *s,n)*(2*s+1.)*(s<2);

    r_val *= _C;

    return r_val;
  }

  inline double
  KernelFuncWendland::dval(double r) const {
    double s = fabs(r/_h);
    const int n = 3;
    double r_val = -5.0*s*pow(1-0.5*s,n)*(s<2);

    r_val *= _C/_h;

    return r_val;
  }

}

#endif
