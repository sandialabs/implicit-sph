#pragma once
#ifndef __FUNCTOR_CURL_H__
#define __FUNCTOR_CURL_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"

namespace LAMMPS_NS {

  // TODO:: This shoudl be further improved; for now happy enough
  template<class PairIsph, 
           template<typename,bool> class Gradient,
           bool AntiSymmetric = false>
  class FunctorOuterCurl : public Gradient<PairIsph,AntiSymmetric> {
  public:
    enum CurlType { ScalarCurl2D, 
                    VectorCurl2D, 
                    Curl3D };

    // ** Scalar Curl 2D
    FunctorOuterCurl(PairIsph *isph, 
                     double **f, 
                     double alpha,
                     double *curl) 
      : Gradient<PairIsph,AntiSymmetric>(isph, f, alpha),
        _curl_type(ScalarCurl2D),
        _scalar_curl(curl) {
      if (_dim != 2)
        _pair->error->all(FLERR, "FunctorOuterCurl:: This constructor is restricted to 2D");
    }

    // ** Vector Curl 2D
    FunctorOuterCurl(PairIsph *isph, 
                     double *f, 
                     double alpha,
                     double **curl) 
      : Gradient<PairIsph,AntiSymmetric>(isph, f, alpha),
        _curl_type(VectorCurl2D),
        _vector_curl(curl) { 
      if (_dim != 2)
        _pair->error->all(FLERR, "FunctorOuterCurl:: This constructor is restricted to 2D");
    }

    // ** Curl 3D
    FunctorOuterCurl(PairIsph *isph, 
                     double **f, 
                     double alpha,
                     double **curl) 
      : Gradient<PairIsph,AntiSymmetric>(isph, f, alpha),
        _curl_type(Curl3D),
        _vector_curl(curl) { 
      //if (_dim != 3)
      //  _pair->error->all(FLERR, "FunctorOuterCurl:: This constructor is restricted to 3D");
    }
      
    void operator()(const int ii);

    double getScalarCurl();
    double* getVectorCurl();

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;

    CurlType _curl_type;

    double *_scalar_curl,_scalar_curl_at_i;
    double **_vector_curl,_vector_curl_at_i[3];
  };

  template<class PairIsph,
           template<typename,bool> class Gradient,
           bool AntiSymmetric>
  inline double
  FunctorOuterCurl<PairIsph,Gradient,AntiSymmetric>::getScalarCurl() {
    return _scalar_curl_at_i;
  }
  template<class PairIsph,
           template<typename,bool> class Gradient,
           bool AntiSymmetric>
  inline double*
  FunctorOuterCurl<PairIsph,Gradient,AntiSymmetric>::getVectorCurl() {
    return _vector_curl_at_i;
  }
    
  template<class PairIsph,
           template<typename,bool> class Gradient,
           bool AntiSymmetric>
  inline void
  FunctorOuterCurl<PairIsph,Gradient,AntiSymmetric>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];

    // initialize with zero
    _scalar_curl_at_i = 0.0;
    memset(_vector_curl_at_i, 0, sizeof(double)*_dim);

    if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
      return;
    
    // compute gradients; 
    // this will compute either scalar gradient or vector gradient (componenet-wise)
    Gradient<PairIsph,AntiSymmetric>::operator()(ii);

    switch (_curl_type) {
    case ScalarCurl2D: {
      double *grad[2] = { this->getVectorGradient(0),
                          this->getVectorGradient(1) };

      _scalar_curl_at_i = grad[1][0] - grad[0][1];

      if (_scalar_curl != NULL)
        _scalar_curl[i] = _scalar_curl_at_i;

      break;
    }
    case VectorCurl2D: {
      double *grad = this->getScalarGradient();

      _vector_curl_at_i[0] =  grad[1];
      _vector_curl_at_i[1] = -grad[0];

      if(_vector_curl != NULL)
        memcpy(_vector_curl[i], _vector_curl_at_i, sizeof(double)*2);

      break;
    }
    case Curl3D: {
      double *grad[3] = { this->getVectorGradient(0),
                          this->getVectorGradient(1),
                          this->getVectorGradient(2) };
    
      _vector_curl_at_i[0] = grad[2][1] - grad[1][2];
      _vector_curl_at_i[1] = grad[0][2] - grad[2][0];
      _vector_curl_at_i[2] = grad[1][0] - grad[0][1];

      if(_vector_curl != NULL)
        memcpy(_vector_curl[i], _vector_curl_at_i, sizeof(double)*3);

      break;
    }
    }
    
  }

}
#endif

