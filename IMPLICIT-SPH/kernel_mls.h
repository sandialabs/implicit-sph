#pragma once
#ifndef __KERNEL_MLS_H__
#define __KERNEL_MLS_H__

#include <iostream>
#include "math.h"
#include "float.h"

// MLS simple 1D kernel

namespace LAMMPS_NS {

  using namespace std;

  class KernelFuncMLS {
  public:
    KernelFuncMLS() { }
    virtual~KernelFuncMLS() { }

    virtual double val(double r, double rth) const {
      const int n = 6;
      return pow(1-fabs(r/rth),n)*(r<rth);
    }
  };

}

#endif
