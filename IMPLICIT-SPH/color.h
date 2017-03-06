#pragma once
#ifndef __COLOR_H__
#define __COLOR_H__

// phase coloring

namespace LAMMPS_NS {

  typedef class ColorFunctionAdamiVar2 ColorFunc;

  class ColorFunction {
  public:
    static const int CorrectedGrad = 1;
    static const int AdamiGrad = 2;

    virtual double val(const int iphase, const double irho,
                       const int jphase, const double jrho) = 0;
    virtual bool yes(const int iphase, const int jphase) = 0;
    virtual int sign(const int iphase, const int jphase) = 0;
    virtual int gradtype() = 0;
  };

  // normal is defined toward the interface
  
  // cij = 1 if particles are different phases
  // cij = 0 if particles are in the same phase
  
  class ColorFunctionCorrected : public ColorFunction {
  public:
    double val(const int iphase, const double irho,
               const int jphase, const double jrho) {
      return (iphase != jphase);
    }
    
    bool yes(const int iphase, const int jphase) {
      return true;
    }
    
    int sign(const int iphase, const int jphase) {
      return (iphase == jphase ? 1.0 : -1.0);       
    }
    
    int gradtype() { 
      return CorrectedGrad;
    }
  };

  class ColorFunctionAdami : public ColorFunction {
  public:
    double val(const int iphase, const double irho,
               const int jphase, const double jrho) {
      return (iphase != jphase)*irho/(irho+jrho);
    }
    
    bool yes(const int iphase, const int jphase) {
      return true;
    }
    
    int sign(const int iphase, const int jphase) {
      return (iphase == jphase ? 1.0 : -1.0);       
    }

    int gradtype() { 
      return AdamiGrad;
    }
  };
  
}

#endif
