#pragma once
#ifndef __FUNCTOR_HELMHOLTZ_ANALYTIC_F_H__
#define __FUNCTOR_HELMHOLTZ_ANALYTIC_F_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph>
  class FunctorOuterHelmholtzAnalyticF : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterHelmholtzAnalyticF(PairIsph *isph, double **f)
      : FunctorOuter<PairIsph>(isph), _f(f) { }
    void enterFor();
    void operator()(const int ii);
    
  protected:

    double **_f;

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  template<class PairIsph> inline void
  FunctorOuterHelmholtzAnalyticF<PairIsph>::enterFor() {
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorOuterHelmholtzAnalyticF:: Filter is not allowed");
  }        

  template<class PairIsph> inline void 
  FunctorOuterHelmholtzAnalyticF<PairIsph>::operator()(const int ii) {
    int i = _ilist[ii];

    const double& x = _x[i][0];
    const double& y = _x[i][1];
    const double& z = _x[i][2];
    const double t = _pair->update->dt*_pair->update->ntimestep;
    const double& nu = _pair->atom->viscosity[i];
    double sx(std::sin(x)), sy(std::sin(y)), sz(std::sin(z)), st(std::sin(t));
    double cx(std::cos(x)), cy(std::cos(y)), cz(std::cos(z)), ct(std::cos(t));


    _f[i][0] = sx*(cy-cz)*(ct+st*(2*nu+cx*(cy-cz)*st)) + sx*st*st*(sy*sy*(cx-cz) + sz*sz*(cx-cy)) + cx*sy*sz;
    _f[i][1] = sy*(cx-cz)*-(ct+st*(2*nu-cy*(cx-cz)*st))+ sy*st*st*(sx*sx*(cy-cz) - sz*sz*(cx-cy)) + sx*cy*sz;
    _f[i][2] = sz*(cx-cy)*(ct+st*(2*nu+cz*(cx-cy)*st)) + sz*st*st*(sx*sx*(cz-cy) + sy*sy*(cz-cx)) + sx*sy*cz;

  }

}

#endif

