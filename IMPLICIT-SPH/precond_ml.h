#pragma once
#ifndef __PRECOND_ML_H__
#define __PRECOND_ML_H__

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "precond.h"

#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"

namespace LAMMPS_NS {

  using namespace std;

  class PrecondWrapper_ML : public PrecondWrapper {
  public:
    PrecondWrapper_ML(MPI_Comm comm) : PrecondWrapper(comm) { }
    virtual~PrecondWrapper_ML() { }

    virtual Teuchos::ParameterList* setParameters(Teuchos::ParameterList *param = NULL);

    // for zoltan re-partition
    void setCoordinates(const int dim, double *x, double *y, double *z);

    virtual void setNullVector(double *n);
    virtual void create();
    virtual void create(const int dim);
    virtual void free();
    virtual Epetra_Operator* getPrecondOperator();
    virtual Thyra::LinearOpBase<double>* getBlockPrecondOperator();

  protected:
    Teuchos::RCP<Teuchos::ParameterList> _param_null;      
    Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> _prec;
    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > _prec_thyra;
  };

  inline Teuchos::ParameterList* 
  PrecondWrapper_ML::setParameters(Teuchos::ParameterList *param) {
    if (param == NULL) {

      _param = Teuchos::rcp(new Teuchos::ParameterList);
      
      _param->set("ML output",10);
      _param->set("max levels",5);
      _param->set("increasing or decreasing","increasing"); 
      _param->set("aggregation: type","Uncoupled");
      _param->set("smoother: type","symmetric Gauss-Seidel");                      
      _param->set("smoother: sweeps",1);                                        
      _param->set("smoother: pre or post","both");
      _param->set("coarse: type","Amesos-KLU");
      //_param->set("coarse: type","Amesos-MUMPS");
    } else {
      if (_param.get() != param)
        _param = Teuchos::rcp(param, false);
    }
    return _param.get();
  }

  inline void
  PrecondWrapper_ML::setCoordinates(const int dim, double *x, double *y, double *z) {
    IS_VALID(_param.get() != NULL, _comm.MyPID(), ">> ParameterList is not yet created");

    // check this pointer for all process ids
    IS_VALID(x != NULL, 0, ">> Coordinates x is null");
    IS_VALID(y != NULL, 0, ">> Coordinates y is null");
    IS_VALID(z != NULL, 0, ">> Coordinates z is null");

    if (dim < 2) {
      _param->set("repartition: enable",0);
      return;
    }
    
    if (_param->isParameter("repartition: enable")) {
      // trust that user provides full parameterlist for repartitioning
      _param->set("repartition: Zoltan dimensions",dim);
    } else {
      _param->set("repartition: enable",1);
      _param->set("repartition: max min ratio",1.3);
      _param->set("repartition: min per proc",500);
      _param->set("repartition: partitioner","Zoltan");
      _param->set("repartition: Zoltan dimensions",dim);
    }
    
    switch (dim) {
    case 3:
      _param->set("z-coordinates",z);
    case 2:
      _param->set("x-coordinates",x);
      _param->set("y-coordinates",y);
    }
  }
  
  inline void
  PrecondWrapper_ML::setNullVector(double *n) {
    if (n == NULL) {
      // recover the original one when it is unset with null
      if (_param_null == Teuchos::null) {
        // do nothing; but at least needs to give error message as it
        // tries to reset params
      } else {
        // recover the original param and delete the temporary
        _param = _param_null;
        _param_null = Teuchos::null;
      }
    } else {
      // copy the _param
      _param_null = Teuchos::rcp(new Teuchos::ParameterList(*_param)); 

      // added or change the param
      _param_null->set("null space: type", "pre-computed");
      _param_null->set("null space: dimension", 1);
      _param_null->set("null space: vectors", n);
      _param_null->set("null space: add default vectors", false);

      // use the same smoother in the coarse level
      _param_null->set("coarse: type", _param->get("smoother: type","symmetric Gauss-Seidel"));
      _param_null->set("coarse: sweeps", _param->get("smoother: sweeps",4));

      // swap rcp's
      auto tmp = _param;
      _param = _param_null;
      _param_null = tmp;
    }
  }

  inline void
  PrecondWrapper_ML::create() {
    // check if A is set
    IS_VALID(_A.get() != NULL, _comm.MyPID(), ">> A is null");
    IS_VALID(_A->Filled(), _comm.MyPID(), ">> A is not filled yet; call FillComplete()");

    _prec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*_A, *_param, true));
  }

  inline void
  PrecondWrapper_ML::create(const int dim) {
    // create diagonal preconditioners; assume that all diagonals are same now
    create(); 

    // once this works, then we need to separate diagonal precond operators
    // according to _Athyra member blocks

    Teuchos::RCP<SwapApply> op_swap = Teuchos::rcp(new SwapApply(_prec.get()));
    Teuchos::RCP<const Thyra::EpetraLinearOp> prec = Thyra::epetraLinearOp(op_swap);

    _prec_thyra = Thyra::defaultBlockedLinearOp<double>();
    _prec_thyra->beginBlockFill(dim, dim);
    // diagonal only, do not put null; it becomes an error in asserting
    for (int k=0;k<dim;++k)
      _prec_thyra->setBlock(k, k, prec);
    _prec_thyra->endBlockFill();
  }
  
  inline void
  PrecondWrapper_ML::free() {
    _prec = Teuchos::null;
    _prec_thyra = Teuchos::null;
  }
   
  inline Epetra_Operator* 
  PrecondWrapper_ML::getPrecondOperator() {
    return (_prec.get());
  }

  inline Thyra::LinearOpBase<double>*
  PrecondWrapper_ML::getBlockPrecondOperator() {
    return (_prec_thyra.get());
  }
  
}

#endif

