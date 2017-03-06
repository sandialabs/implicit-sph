#pragma once
#ifndef __FUNCTOR_ADVANCE_TIME_EXACT_SOLUTION_H__
#define __FUNCTOR_ADVANCE_TIME_EXACT_SOLUTION_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"
#include "pair_for.h"

namespace LAMMPS_NS {

  // ---------------------------------------------------------------------
  class FunctorOuterExactSolution_Vortex {
  public:
    FunctorOuterExactSolution_Vortex(int *ilist, 
                                     double **xp, double **v, double *p) 
      : _ilist(ilist),_xp(xp),_v(v),_p(p) { }
    void enterFor() { }
    void operator()(const int ii);
    void exitFor() { }    

  protected:
    int *_ilist;
    double **_xp,**_v,*_p;
  };
  
  inline void 
  FunctorOuterExactSolution_Vortex::operator()(const int ii) {
    int i = _ilist[ii];

    double *v = _v[i], *x = _xp[i];
    const double c = 1.0;
    
    v[0] =  sin(0.5*x[0]*c)*cos(0.5*x[1]*c);
    v[1] = -cos(0.5*x[0]*c)*sin(0.5*x[1]*c);
    v[2] =  0.0;
    
    double &p = _p[i];
    p = 0.0;
  }


  // ---------------------------------------------------------------------
  typedef FunctorOuterExactSolution_Vortex FunctorExactSolution;
  // ---------------------------------------------------------------------
  template<class PairIsph>
  class FunctorOuterAdvanceTimeExactSolutionBegin : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterAdvanceTimeExactSolutionBegin(PairIsph *isph, 
                                              double dt, 
                                              double **v, 
                                              double **vnp1,
                                              double *p,
                                              double **xtmp) 
      : FunctorOuter<PairIsph>(isph),
        _dt(dt),
        _v(v),
        _vnp1(vnp1),
        _p(p),
        _xtmp(xtmp) { }
    void enterFor();
    void operator()(const int ii);
    void exitFor();

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;
    double _dt,**_v,**_vnp1,*_p,**_xtmp;
  };

  template<class PairIsph>
  inline void
  FunctorOuterAdvanceTimeExactSolutionBegin<PairIsph>::enterFor() {
    FunctorExactSolution functor(_ilist, _x, _v, _p);
    PairFor(functor, this->getNumberOfWork());
  }  

  template<class PairIsph>
  inline void
  FunctorOuterAdvanceTimeExactSolutionBegin<PairIsph>::operator()(const int ii) {
    int i = _ilist[ii]; 
    int itype = _type[i];

    if (_pair->getParticleKind(itype) == PairIsph::Fluid)
      for (int k=0;k<_dim;++k)
        _xtmp[i][k] = _x[i][k] + _dt*_v[i][k];
  }

  template<class PairIsph>
  inline void
  FunctorOuterAdvanceTimeExactSolutionBegin<PairIsph>::exitFor() {
    FunctorExactSolution functor(_ilist, _xtmp, _vnp1, _p);
    PairFor(functor, this->getNumberOfWork());
  }  

  // ---------------------------------------------------------------------
  template<class PairIsph>
  class FunctorOuterAdvanceTimeExactSolutionEnd : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterAdvanceTimeExactSolutionEnd(PairIsph *isph, 
                                            double dt, 
                                            double **v, 
                                            double **vnp1,
                                            double *p)
      : FunctorOuter<PairIsph>(isph),
        _dt(dt),
        _v(v),
        _vnp1(vnp1),
        _p(p) { }
    void operator()(const int ii);

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;
    double _dt,**_v,**_vnp1,*_p;
  };

  template<class PairIsph>
  inline void
  FunctorOuterAdvanceTimeExactSolutionEnd<PairIsph>::operator()(const int ii) {
    int i = _ilist[ii]; 
    int itype = _type[i];

    if (_pair->getParticleKind(itype) == PairIsph::Fluid)
      for (int k=0;k<_dim;++k)
        _x[i][k] += 0.5*_dt*(_v[i][k] + _vnp1[i][k]);
  }

}

#endif

