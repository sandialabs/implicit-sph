#pragma once
#ifndef __PRECOND_IFPACK_H__
#define __PRECOND_IFPACK_H__

#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h" 

#include "precond.h"

namespace LAMMPS_NS {

  using namespace std;

  class PrecondWrapper_Ifpack : public PrecondWrapper {
  public:
    PrecondWrapper_Ifpack(MPI_Comm comm) : PrecondWrapper(comm) { }
    virtual~PrecondWrapper_Ifpack() { }

    virtual Teuchos::ParameterList* setParameters(Teuchos::ParameterList *param = NULL);
    virtual void create();
    virtual void free();
    virtual Epetra_Operator* getPrecondOperator();

  protected:
    Teuchos::RCP<Ifpack_Preconditioner> _prec;
  };

  inline Teuchos::ParameterList*
  PrecondWrapper_Ifpack::setParameters(Teuchos::ParameterList *param) {
    if (param == NULL) {
      _param = Teuchos::rcp(new Teuchos::ParameterList);
      
      // Ifpack native params
      _param->set("fact: drop tolerance", 1e-9);
      _param->set("fact: level-of-fill", 0);
      // the combine mode is on the following:
      // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
      // Their meaning is as defined in file Epetra_CombineMode.h
      _param->set("schwarz: combine mode", "Add");

      // wrapping param
      _param->set("Precond Type","ILU");
      _param->set("Overlap Level",1); // must be >= 0
      
    } else {
      if (_param.get() != param)
        _param = Teuchos::rcp(param, false);
    }
    return _param.get();
  }

  inline void
  PrecondWrapper_Ifpack::create() {
    // check if A is set
    IS_VALID(_A.get() != NULL, _comm.MyPID(), ">> A is null");
    
    int ierr;
    IS_VALID(_A->Filled(), _comm.MyPID(), ">> A is not filled yet; call FillComplete()");

    Ifpack Factory; 
    _prec = Teuchos::rcp(Factory.Create(_param->get("Precond Type","ILU"), 
                                        _A.get(), 
                                        _param->get("Overlap Level",1)));
    IS_VALID(_prec != Teuchos::null, _comm.MyPID(), ">> Error in Ifpack, Factor.Create() return null: ");

    // load default or do nothing if it sets up
    setParameters(_param.get());

    // send native ifpack param to ifpack preconditioner
    ierr = _prec->SetParameters(*_param); 
    IS_VALID(ierr >=0, _comm.MyPID(), ">> Error in Ifpack, SetParameters(): " << ierr);

    ierr = _prec->Initialize(); IS_VALID(ierr >=0, _comm.MyPID(), ">> Error in Ifpack, Initialize(): " << ierr);
    ierr = _prec->Compute();    IS_VALID(ierr >=0, _comm.MyPID(), ">> Error in Ifpack, Compute(): " << ierr);
  }

  inline void
  PrecondWrapper_Ifpack::free() {
    _prec = Teuchos::null;
  }
   
  inline Epetra_Operator*
  PrecondWrapper_Ifpack::getPrecondOperator() {
    return _prec.get();
  }
  
}

#endif

