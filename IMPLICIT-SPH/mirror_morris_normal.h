#pragma once
#ifndef __MIRROR_MORRIS_NORMAL_H__
#define __MIRROR_MORRIS_NORMAL_H__

#include "utils.h"
#include "mirror.h"

namespace LAMMPS_NS {

  class MirrorMorrisNormal : public MirrorFunction {
  public:
    MirrorMorrisNormal(unsigned int dim, 
                       int *type, double **h,
                       double safe,
                       double **coords, double *bd_coord, double **normal)
      : MirrorFunction(dim, type, h),  
        _safe(safe),
        _x(coords), 
        _bd_coord(bd_coord),
        _n(normal), 
        _xi_at_i(0), 
        _xi_at_j(0), 
        _i(0),
        _type_at_i(0),
        _type_at_j(0) { }

    void setFreeParticle(const int i);
    void setMirrorParticle(const int j);
    double computeMirrorCoefficient(const double r);
    
  protected:
    double _safe;
    double **_x, *_bd_coord, **_n;
    double _xi_at_i, _xi_at_j;
    int _i, _type_at_i, _type_at_j;
    
    using MirrorFunction::_type; // type array
    using MirrorFunction::_h;    // smoothing length
  };
  
  inline void MirrorMorrisNormal::setFreeParticle(const int i) {
    _i = i;
    _type_at_i = _type[i];
    _xi_at_i = util.dotVectors(_dim, &_x[i][0], &_n[i][0]);
  }
  inline void MirrorMorrisNormal::setMirrorParticle(const int j) {
    //cout << _bd_coord[_i] << " "<<endl;
    _type_at_j = _type[j];
    _xi_at_j = util.dotVectors(_dim, &_x[j][0], &_n[_i][0]);
  }
  inline double MirrorMorrisNormal::computeMirrorCoefficient(const double r) {
    double d_i = fabs(_xi_at_i - _bd_coord[_i]) + r*1e-8;
    double d_j = fabs(_xi_at_j - _bd_coord[_i]);
    double h   = LOOKUP(_h, _type_at_i, _type_at_j);

    return (1.0+d_j/max(d_i, _safe*h));
  }

}

#endif
