#pragma once
#ifndef __FUNCTOR_CURLCURL_H__
#define __FUNCTOR_CURLCURL_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"
#include "pair_for.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename,bool> class Gradient,
           template<typename,bool> class GradientInner = Gradient,
           bool AntiSymmetric = false>
  class FunctorOuterCurlCurl : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterCurlCurl(PairIsph *isph, 
                         double **f, 
                         double alpha = 1.0,
                         double **curlcurl = NULL)
      : FunctorOuter<PairIsph>(isph),
        _op_curl(NULL),        
        _f(f),
        _alpha(alpha),
        _curlcurl(curlcurl) { }
    virtual~FunctorOuterCurlCurl() {
      if (_op_curl != NULL)
        delete _op_curl;
    }

    void enterFor();
    void operator()(const int ii);
    void exitFor();
    double* getCurlCurl();
    
  protected:
    FUNCTOR_REMOVE_THIS_POINTER;    
    FunctorOuterCurl<PairIsph,Gradient,AntiSymmetric> *_op_curl;
    
    double **_f,_alpha,**_curlcurl;
  };

  template<class PairIsph,
           template<typename,bool> class Gradient,
           template<typename,bool> class GradientInner,
           bool AntiSymmetric>
  inline double*
  FunctorOuterCurlCurl<PairIsph,Gradient,GradientInner,AntiSymmetric>::getCurlCurl() {
    return _op_curl->getVectorCurl();
  }

  template<class PairIsph, 
           template<typename,bool> class Gradient,
           template<typename,bool> class GradientInner,
           bool AntiSymmetric>
  inline void
  FunctorOuterCurlCurl<PairIsph,Gradient,GradientInner,AntiSymmetric>::enterFor() {

    // use the all-all filter for the first curl;
    // first culr should be universaly defined to compute the second curl
    FilterAny filter_op_curl;

    switch (_dim) {
    case 2: {
      // need to compute for all-all 
      FunctorOuterCurl<PairIsph,GradientInner,AntiSymmetric> op_curl(_pair, _f, _alpha, _pair->work);
      op_curl.setFilter(&filter_op_curl);

      PairFor(op_curl, op_curl.getNumberOfWork());

      _pair->comm_variable = PairIsph::WorkScalar;
      _pair->comm_forward = 1;
      _pair->comm->forward_comm_pair(_pair);

      _op_curl = new FunctorOuterCurl<PairIsph,Gradient,AntiSymmetric>(_pair, _pair->work, 1.0, _curlcurl);
      break;
    }
    case 3: {
      // call 3D curl
      FunctorOuterCurl<PairIsph,GradientInner,AntiSymmetric> op_curl(_pair, _f, _alpha, _pair->work3);
      op_curl.setFilter(&filter_op_curl);

      PairFor(op_curl, op_curl.getNumberOfWork());
      
      _pair->comm_variable = PairIsph::WorkVector;
      _pair->comm_forward = 3;
      _pair->comm->forward_comm_pair(_pair);

      _op_curl = new FunctorOuterCurl<PairIsph,Gradient,AntiSymmetric>(_pair, _pair->work3, 1.0, _curlcurl);
      break;
    }
    }

    _op_curl->setFilter(_filter);
  }

  template<class PairIsph,
           template<typename,bool> class Gradient,
           template<typename,bool> class GradientInner,
           bool AntiSymmetric>
  inline void
  FunctorOuterCurlCurl<PairIsph,Gradient,GradientInner,AntiSymmetric>::operator()(const int ii) {
    (*_op_curl)(ii);
  }      

  template<class PairIsph, 
           template<typename,bool> class Gradient,
           template<typename,bool> class GradientInner,
           bool AntiSymmetric>
  inline void
  FunctorOuterCurlCurl<PairIsph,Gradient,GradientInner,AntiSymmetric>::exitFor() {
    if (_op_curl != NULL) 
      delete _op_curl;
    _op_curl = NULL;
  }
}
#endif

