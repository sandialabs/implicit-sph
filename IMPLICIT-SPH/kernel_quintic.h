#pragma once
#ifndef __KERNEL_QUINTIC_H__
#define __KERNEL_QUINTIC_H__

#include <iostream>
#include "math.h"
#include "float.h"
#include "kernel.h"

// Quintic spline kernel

namespace LAMMPS_NS {

  class KernelFuncQuintic : public KernelFunction {
  protected:
    void setSmoothingLength(double h);  
  public:
    KernelFuncQuintic(unsigned int dim) : KernelFunction(dim) { }
    KernelFuncQuintic(unsigned int dim, double h) : KernelFunction(dim) {
      setSmoothingLength(h);
    }
    double val (double r) const;
    double dval(double r) const;
  };

  using namespace std;

  inline void 
  KernelFuncQuintic::setSmoothingLength(double h) {
    // if the given h is not same as cached _h, then _C is recomputed.
    if (fabs(_h - h) > DBL_EPSILON || _h == 0.0) {
      // smoothing parameter is cached
      _h = h;
      
      // set normalization parameter
      if (_dim == 3) 
        _C = 14.0/(pow(_h,3)*1745.0*M_PI);
      else 
        _C =  7.0/(pow(_h,2)* 478.0*M_PI);
    }
  }
  
  inline double 
  KernelFuncQuintic::val(double r) const {
    double r_val = 0.0, s = fabs(r/_h);
    const int n = 5;
    /*
      if (s <= 1.0)
      r_val = pow(3.0-s,n) - 6.0*pow(2.0-s,n) + 15.0*pow(1.0-s,n);
      else if (s <= 2.0)
      r_val = pow(3.0-s,n) - 6.0*pow(2.0-s,n);
      else if (s <= 3.0)
      r_val = pow(3.0-s,n);
    */
    switch (static_cast<int>(floor(s))) {
    case 0:    r_val += (15.0*pow(1.0-s,n));
    case 1:    r_val -= (6.0*pow(2.0-s,n));
    case 2:    r_val += (     pow(3.0-s,n));
    }

    r_val *= _C;

    return r_val;
  }

  inline double 
  KernelFuncQuintic::dval(double r) const {
    double r_val = 0.0, s = fabs(r/_h);
    const int n = 4;

    switch (static_cast<int>(floor(s))) {
    case 0:    r_val -= (75.0*pow((1-s),n));
    case 1:    r_val += (30.0*pow((2-s),n));
    case 2:    r_val -= (  5*pow((3-s),n));
    }

    r_val *= _C/_h;

    return r_val;
  }

}

#endif
