#pragma once
#ifndef __MIRROR_MORRIS_HOLMES_H__
#define __MIRROR_MORRIS_HOLMES_H__

#include "utils.h"
#include "mirror.h"

namespace LAMMPS_NS {

  class MirrorMorrisHolmes : public MirrorFunction {
  public:  
    MirrorMorrisHolmes(unsigned int dim, 
                       int *type, double **h,
                       double safe,
                       double *pnd, double *vfrac) 
      : MirrorFunction(dim, type, h),
        _safe(safe),
        _pnd(pnd),
        _vfrac(vfrac),
        _xi_at_i(0),
        _xi_at_j(0),
        _type_at_i(0),
        _type_at_j(0) { }

    void setFreeParticle(const int i);
    void setMirrorParticle(const int j);
    double computeMirrorCoefficient(const double r);
    
  protected:
    double _safe;
    double *_pnd, *_vfrac;
    double _xi_at_i, _xi_at_j;
    int _type_at_i, _type_at_j;

    using MirrorFunction::_type; // type array
    using MirrorFunction::_h;    // smoothing length
  };
  
  inline void MirrorMorrisHolmes::setFreeParticle(const int i) { 
    _xi_at_i = _pnd[i]*_vfrac[i];
    _type_at_i = _type[i];
  }
  inline void MirrorMorrisHolmes::setMirrorParticle(const int j) { 
    _xi_at_j = _pnd[j]*_vfrac[j];
    _type_at_j = _type[j];
  }
  inline double MirrorMorrisHolmes::computeMirrorCoefficient(const double r) { 
    double d_i = 2.0*r*(_xi_at_i-0.5) + ISPH_EPSILON;
    double d_j = 2.0*r*(_xi_at_j-0.5) + ISPH_EPSILON;  
    double h   = LOOKUP(_h, _type_at_i, _type_at_j);

    return (1.0+d_j/max(d_i, _safe*h));
  }

}

#endif
