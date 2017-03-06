#pragma once
#ifndef __KERNEL_CUBIC_H__
#define __KERNEL_CUBIC_H__

#include <iostream>
#include "math.h"
#include "float.h"
#include "kernel.h"

// Cubic spline kernel

namespace LAMMPS_NS {

  class KernelFuncCubic: public KernelFunction {  
  protected:
    void setSmoothingLength(double h);
  public:  
    KernelFuncCubic(unsigned int dim) : KernelFunction(dim) {}
    KernelFuncCubic(unsigned int dim, double h) : KernelFunction(dim) {
      setSmoothingLength(h); 
    }
    double val (double r) const;
    double dval(double r) const;
  };

  using namespace std;

  inline void
  KernelFuncCubic::setSmoothingLength(double h) {
    // if the given h is not same as cached _h, then _C is recomputed.
    if (fabs(_h - h) > DBL_EPSILON || _h == 0.0) {
      // smoothing parameter is cached
      _h = h;

      // set normalization parameter
      if (_dim == 3) 
        _C = 1.0/(pow(_h,3)*M_PI);
      else 
        _C =  10.0/(pow(_h,2)* 7.0*M_PI);
    }
  }
  
  inline double 
  KernelFuncCubic::val(double r) const {
    double r_val = 0.0, s = fabs(r/_h);

    switch (static_cast<int>(floor(s))) {
    case 0:    r_val = 1.0 - 0.75*(2-s)*s*s; break;
    case 1:    r_val = 0.25*pow(2.0-s,3);
    }

    r_val *= _C;

    return r_val;
  }

  inline double
  KernelFuncCubic::dval(double r) const {
    double r_val = 0.0, s = fabs(r/_h);

    switch (static_cast<int>(floor(s))) {
    case 0:    r_val = (2.25*s-3)*s; break;
    case 1:    r_val = -0.75*pow(2-s,2);
    }

    r_val *= _C/_h;

    return r_val;
  }

}

#endif
