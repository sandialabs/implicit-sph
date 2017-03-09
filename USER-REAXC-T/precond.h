#pragma once
#ifndef __PRECOND_H__
#define __PRECOND_H__

#include "Epetra_MpiComm.h" 
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_ParameterList.hpp" 
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"

namespace LAMMPS_NS {

  using namespace std;

  class PrecondWrapper { 
  protected:
    Epetra_MpiComm _comm;   

    Teuchos::RCP<Teuchos::ParameterList> _param;
    Teuchos::RCP<Epetra_CrsMatrix> _A;
    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > _Athyra;

  public:
    PrecondWrapper(MPI_Comm comm) : _comm(comm) { }
    virtual~PrecondWrapper() { }

    virtual void setMatrix(Epetra_CrsMatrix *A) { 
      if (A != NULL)
        _A = Teuchos::rcp(A, false);
    }
    virtual void setBlockMatrix(Thyra::PhysicallyBlockedLinearOpBase<double> *Athyra) {
      if (Athyra != NULL)
        _Athyra = Teuchos::rcp(Athyra, false);
    }
    virtual Teuchos::ParameterList* setParameters(Teuchos::ParameterList *param = NULL) { 
      return _param.get(); 
    }
    virtual void setNullVector(double *n) { return; }
    virtual void create() { return; }
    virtual void create(const int dim) { return; } 
    virtual void free() { return; } 
    virtual Epetra_Operator* getPrecondOperator() { return NULL; }
    virtual Thyra::LinearOpBase<double>* getBlockPrecondOperator() { return NULL; }
  };

  // -------------------------------------------------------------------
  class SwapApply: public Epetra_Operator {      
  public:
    SwapApply(Epetra_Operator *op) 
      : _op(op) {;}
    virtual~SwapApply() {};

    int SetUseTranspose(bool UseTranspose);
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;    
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
    double NormInf() const;

    const char* Label() const;
    bool UseTranspose() const;
    bool HasNormInf() const;

    const Epetra_Comm& Comm() const;
    const Epetra_Map& OperatorDomainMap() const;
    const Epetra_Map& OperatorRangeMap() const;

  private:
    Epetra_Operator *_op;

  };

  inline int 
  SwapApply::SetUseTranspose(bool UseTranspose) {
    return _op->SetUseTranspose(UseTranspose);
  }
  inline int 
  SwapApply::Apply(const Epetra_MultiVector& X, 
                   Epetra_MultiVector& Y) const {
    return _op->ApplyInverse(X,Y);
  }
  inline int 
  SwapApply::ApplyInverse(const Epetra_MultiVector& X, 
                          Epetra_MultiVector& Y) const {
    return _op->Apply(X,Y);
  }
  inline double 
  SwapApply::NormInf() const {
    return _op->NormInf();
  }

  inline const char* 
  SwapApply::Label() const { 
    return "Swap Apply";
  }
  inline bool 
  SwapApply::UseTranspose() const {
    return _op->UseTranspose();
  }
  inline bool 
  SwapApply::HasNormInf() const {
    return _op->HasNormInf();
  }

  inline const Epetra_Comm& 
  SwapApply::Comm() const {
    return _op->Comm();
  }
  inline const Epetra_Map& 
  SwapApply::OperatorDomainMap() const {
    return _op->OperatorDomainMap();
  }
  inline const Epetra_Map& 
  SwapApply::OperatorRangeMap() const {
    return _op->OperatorRangeMap();
  }

}

#endif

