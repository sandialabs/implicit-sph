#pragma once
#ifndef __PAIRWISE_FORCE_H__
#define __PAIRWISE_FORCE_H__

// Pairwise force interaction model
#include <iostream>
#include "math.h"
#include "float.h"
#include "kernel.h"

namespace LAMMPS_NS {

  class PairwiseForceFunction {
  protected:  
    unsigned int _dim;
    double _A, _h, _h0, _eps, _eps0;

    virtual void setCutLength(double h) = 0;  
  public:
    PairwiseForceFunction(unsigned int dim) 
      : _dim(dim), 
        _A(0.0),
        _h(0.0),
        _h0(0.0),
        _eps(0.0),
        _eps0(0.0)
    { }
    virtual~PairwiseForceFunction() { }

    virtual double val (double sij, double r) const = 0;
    virtual double val (double sij, double r, double h) { 
      setCutLength(h); return val(sij, r);  
    }
    virtual double lambda() const { return 0.0; }
	virtual double xi() const { return 0.0; }
  };

  class PairwiseForceFunction_TartakovskyMeakin : public PairwiseForceFunction {
  protected:
    void setCutLength(double h) {
      if (fabs(_h - h) > DBL_EPSILON || _h == 0.0) {
        // smoothing parameter is cached
        _h = h;
      }
    }

  public:
    PairwiseForceFunction_TartakovskyMeakin(unsigned int dim) 
      : PairwiseForceFunction(dim) { }

    double val(const double s, const double r) const {
      // F(s, r) = -s cos (3pi/2h r) * (r<=h)
      return -s*cos(4.71238898038469/_h*r)*(r<=_h);
    }
  };

  class PairwiseForceFunction_TartakovskyPanchenkoVar1 : public PairwiseForceFunction {
  protected:
    void setCutLength(double h) {
      if (fabs(_h - h) > DBL_EPSILON || _h == 0.0) {
        // smoothing parameter is cached
        _h = h;

        _eps  = _h/3.5;
        _eps0 = _eps/2.0;

        // set parameters
        if (_dim == 3) {
          _A = 8.0;
        } else {
          _A = 4.0;
        }
      }
    }
  public:
    PairwiseForceFunction_TartakovskyPanchenkoVar1(unsigned int dim) 
      : PairwiseForceFunction(dim) { }

    double psi(const double r, const double eps) const {
      // psi(r, eps) = exp(-r^2/2eps^2)
      return exp(-pow(r,2)/pow(eps,2)/2.0);
    }
    double val(const double s, const double r) const {
      // F(s, r) = s (-A psi(r, eps0) + psi(r, eps));
      return s*(-_A*psi(r,_eps0) + psi(r,_eps));
    }
  };

  class PairwiseForceFunction_TartakovskyPanchenkoVar2 : public PairwiseForceFunction {
  protected:
    void setCutLength(double h) {
      if (fabs(_h - h) > DBL_EPSILON || _h == 0.0) {
        // smoothing parameter is cached
        _h = h;

        _eps  = _h/3.5;
        _eps0 = _eps/2.0;

        // set parameters
        if (_dim == 3) {
          _A = 16.0;
        } else {
          _A = 8.0;
        }
      }
    }
  public:
    PairwiseForceFunction_TartakovskyPanchenkoVar2(unsigned int dim) 
      : PairwiseForceFunction(dim) { }
    double psi(const double r, const double eps) const {
      // psi(r, eps) = exp(-r^2/2eps^2)
      return exp(-pow(r,2)/pow(eps,2)/2.0);
    }
    double val(const double s, const double r) const {
      // F(s, r) = s r (-A psi(r, eps0) + psi(r, eps));
      return s*r*(-_A*psi(r,_eps0) + psi(r,_eps));
    }
  };

}

#endif
