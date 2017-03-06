#pragma once
#ifndef __KERNEL_H__
#define __KERNEL_H__

// 1D SPH kernel base class

namespace LAMMPS_NS {

  class KernelFunction {
  protected:  
    unsigned int _dim;
    double _h, _C; 

    virtual void setSmoothingLength(double h) = 0;  
  public:
    KernelFunction(unsigned int dim) : _dim(dim), _h(0.0), _C(0.0) {}
    virtual~KernelFunction() { }

    virtual double val (double r) const = 0;
    virtual double dval(double r) const = 0;

    virtual double val (double r, double h) { setSmoothingLength(h); return val(r);  }
    virtual double dval(double r, double h) { setSmoothingLength(h); return dval(r); }
  };

}

#endif
