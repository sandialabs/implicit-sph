#pragma once
#ifndef __FUNCTOR_POISSON_BOLTZMAN_EXTRA_F_H__
#define __FUNCTOR_POISSON_BOLTZMAN_EXTRA_F_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "RTC_FunctionRTC.hh"  

namespace LAMMPS_NS {

  template<class PairIsph>
  class FunctorOuterPoissonBoltzmannExtraF : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterPoissonBoltzmannExtraF(PairIsph *isph, 
                                      double *psi, double *psi0, double *eps,
                                      Epetra_Vector *f,
                                      Teuchos::ParameterList param) 
      : FunctorOuter<PairIsph>(isph),_psi(psi),_psi0(psi0), _eps(eps), _f(f),_param(param) { }
    void enterFor();
    void operator()(const int ii);
    
  protected:

    double *_psi,*_psi0, *_eps;
    Epetra_Vector *_f; 

    Teuchos::ParameterList _param, _vars, _func;
    Teuchos::RCP<PG_RuntimeCompiler::Function> _extra_f;

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  template<class PairIsph> inline void
  FunctorOuterPoissonBoltzmannExtraF<PairIsph>::enterFor() {
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorPoissonBoltzmannExtraExtraF:: Filter is not allowed");

    // create the function object
    _extra_f = Teuchos::rcp(new PG_RuntimeCompiler::Function(0, "Extra F")); 

    // error check
    if (!_param.isSublist("Variable List") ||
        !_param.isSublist("Function List")) 
      _pair->error->all(FLERR, "FunctorPoissonBoltzmannExtraF:: Variable and/or Function lists are empty");

    // construct variable list
    _vars = _param.sublist("Variable List");

    long long int cnt = 0;
    for (auto iter=_vars.begin();iter!=_vars.end();++iter) {
      string name = _vars.name(iter);
      string type = _vars.get(name, "");
      _extra_f->addVar(type, name);
      _vars.set(name, to_string(cnt));

      // initialize all variables zero; necessary to work with RTC
      _extra_f->varValueFill(stoi(_vars.get(name,"")), 0);
      ++cnt;
    }

    // load funciton list and error check
    _func = _param.sublist("Function List");
    if (_func.get("f","none") == "none") {
      _pair->error->all(FLERR, "FunctorPoissonBoltzmannExtraF:: Function list does not contain an expression for f");
    } else {
      _extra_f->addBody(_func.get("f",""));
    }      
  }        

  template<class PairIsph> inline void 
  FunctorOuterPoissonBoltzmannExtraF<PairIsph>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];

    if (_pair->getParticleKind(itype) & PairIsph::Solid) {
      // do nothing
    } else {
      _extra_f->varValueFill(stoi(_vars.get("pt.x","")), _x[i][0]);
      _extra_f->varValueFill(stoi(_vars.get("pt.y","")), _x[i][1]);
      _extra_f->varValueFill(stoi(_vars.get("pt.z","")), _x[i][2]);  
      _extra_f->execute();
      (*_f)[i] += _extra_f->getValueOfVar("f");  
    }
  }

}

#endif

