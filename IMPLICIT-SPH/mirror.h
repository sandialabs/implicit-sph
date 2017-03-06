#pragma once
#ifndef __MIRROR_H__
#define __MIRROR_H__


namespace LAMMPS_NS {

  typedef class MirrorFunction MirrorNothing;

  class MirrorFunction {
  public:
    MirrorFunction(unsigned int dim, int *type, double **h) 
      : _dim(dim), _type(type), _h(h) {}
    virtual~MirrorFunction() { }
  
    virtual void setFreeParticle(const int i) { /* do nothing */ }
    virtual void setMirrorParticle(const int j) { /* do nothing */ }
    virtual double computeMirrorCoefficient(const double r) { return 1.0; }
  
  protected:
    unsigned int _dim;
    int *_type;
    double **_h;
  };

}

#endif
