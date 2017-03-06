#pragma once
#ifndef __FUNCTOR_H__
#define __FUNCTOR_H__

#include "filter.h"

namespace LAMMPS_NS {

  template<bool AntiSymmetric>
  double sphOperator(const double fi, const double fj);

  template<>
  inline double sphOperator<true>(const double fi, const double fj) {
    return (fi + fj);
  }

  template<>
  inline double sphOperator<false>(const double fi, const double fj) {
    return (fj - fi);
  }
  
#define FUNCTOR_REMOVE_THIS_POINTER                  \
  using FunctorOuter<PairIsph>::_pair;               \
  using FunctorOuter<PairIsph>::_filter;             \
  using FunctorOuter<PairIsph>::_vfrac;              \
  using FunctorOuter<PairIsph>::_x;                  \
  using FunctorOuter<PairIsph>::_type;               \
  using FunctorOuter<PairIsph>::_tag;                \
  using FunctorOuter<PairIsph>::_inum;               \
  using FunctorOuter<PairIsph>::_ilist;              \
  using FunctorOuter<PairIsph>::_numneigh;           \
  using FunctorOuter<PairIsph>::_firstneigh;         \
  using FunctorOuter<PairIsph>::_dim;                \
  using FunctorOuter<PairIsph>::_dimsq;              \
  using FunctorOuter<PairIsph>::_dimL;             
  
  template<class PairIsph>
  class FunctorOuter {
  public:

    FunctorOuter(PairIsph *isph, Filter *filter = NULL);
    virtual~FunctorOuter();
    
    virtual void setFilter(Filter *filter);
    virtual Filter* getFilter() const;

    virtual size_t getNumberOfWork() const;
    virtual void enterFor();
    virtual void exitFor();
    virtual void operator()(const int ii) = 0;

  protected:

    // pointer to pair isph class
    PairIsph *_pair;

    // filter
    Filter *_filter;
    
    // atom variabels
    double *_vfrac, **_x;
    int *_type, *_tag;
    
    // list structure for i
    int _inum, *_ilist, *_numneigh, **_firstneigh;

    // domain dimension
    int _dim, _dimsq, _dimL;
  };
  
  template<class PairIsph> inline 
  FunctorOuter<PairIsph>::FunctorOuter(PairIsph *isph, Filter *filter) {
    // need error check for the pointer
    _pair       = isph;

    _filter     = filter;

    _vfrac      = _pair->atom->vfrac;
    _x          = _pair->atom->x;
    _type       = _pair->atom->type;
    _tag        = _pair->atom->tag;
    
    _inum       = _pair->list->inum;
    _ilist      = _pair->list->ilist;
    _numneigh   = _pair->list->numneigh;
    _firstneigh = _pair->list->firstneigh;

    _dim        = _pair->domain->dimension;
    _dimsq      = _dim*_dim;
    _dimL       = _dim*(_dim+1)/2;
  }
  template<class PairIsph> inline 
  FunctorOuter<PairIsph>::~FunctorOuter() { }

  template<class PairIsph> inline size_t 
  FunctorOuter<PairIsph>::getNumberOfWork() const { return _inum; }

  template<class PairIsph> inline void 
  FunctorOuter<PairIsph>::enterFor() { }

  template<class PairIsph> inline void 
  FunctorOuter<PairIsph>::exitFor() { }

  template<class PairIsph> inline void 
  FunctorOuter<PairIsph>::setFilter(Filter *filter) { 
    _filter = filter;
  }
  template<class PairIsph> inline Filter* 
  FunctorOuter<PairIsph>::getFilter() const { 
    return _filter;
  }
}

#endif

