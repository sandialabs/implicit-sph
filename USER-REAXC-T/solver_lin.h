#pragma once
#ifndef __SOLVER_LIN_H__
#define __SOLVER_LIN_H__

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h" 
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

#include "Teuchos_ParameterList.hpp" 
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"

#include "precond.h"

namespace LAMMPS_NS {

  // -------------------------------------------------------------------
  class SolverLin {
  public:
    enum SolutionInitType { Random, Zero, Value };

  public:
    SolverLin(MPI_Comm &comm);
    virtual~SolverLin();
    
    // solver setting
    int createLinearMap(int num_global_nodes, int index_base);
    int createNodalMap(int num_local_nodes, int *gID);
    
    int createLoadMultiVector(int num_vectors);
    int createLoadMultiVector(double *b, int lda, int num_vectors);
    
    int createSolutionMultiVector(int num_vectors);
    int createSolutionMultiVector(double *x, int lda, int num_vectors);
    
    int createNullVector();

    int createBlockMatrix(const int dim, const char *name);
    int freeBlockMatrix();

    void setNullVectorMask(Epetra_IntSerialDenseVector *mask);

    void setMatrixIsBlocked(const bool is_blocked);
    void setMatrixIsSingular(const bool is_singular);
    
    void setNodalMap(Epetra_Map *map);
    void setMatrix(Epetra_CrsMatrix *A);

    void setBlockBegin();
    void setBlock(const int i, const int j, const Epetra_CrsMatrix *A) ;
    void setBlockEnd();

    void setInitialSolution(SolutionInitType init, double val = 0.0);
    
    Teuchos::RCP<Epetra_Map> getNodalMap() const;
    Teuchos::RCP<Epetra_MultiVector> getLoadMultiVector();
    Teuchos::RCP<Epetra_MultiVector> getSolutionMultiVector();    
    Teuchos::RCP<Epetra_Vector> getNullVector();    
    
    virtual void setParameters(Teuchos::ParameterList *param = NULL) { }
    virtual int solveProblem(PrecondWrapper *prec = NULL, const char* name = NULL) { return 0; };
    virtual int solveBlockProblem(PrecondWrapper *prec = NULL, const char* name = NULL) { return 0; };

  protected:
    virtual int createBlockVector(const Teuchos::RCP<Epetra_MultiVector> v,
                                  Teuchos::RCP<Thyra::ProductMultiVectorBase<double> > &vthyra);

    // base communicator for Epetra
    Epetra_MpiComm _comm;
    
    // map is either created or given; 
    // when they are set from outside, rcp does not have the ownership 
    Teuchos::RCP<Epetra_Map> _map;
    
    // main parameter and sublists
    Teuchos::RCP<Teuchos::ParameterList> _param;
    
    // delivered from problem and set
    Teuchos::RCP<Epetra_CrsMatrix> _A;  

    // delivered from problem and set; 3x3 block matrix form
    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > _Athyra;
    int _dim;
    bool _is_blocked;

    // interface to solution vector (rhs) used in Lin
    Teuchos::RCP<Epetra_MultiVector> _x, _b; 

    // null vector for projection; use a single vector here
    Teuchos::RCP<Epetra_Vector> _n; 
    Teuchos::RCP<Epetra_IntSerialDenseVector> _null_mask;
    bool _is_singular;
  };

  // -------------------------------------------------------------------
  class PoissonProjection : public Epetra_Operator {
  public:
    
    PoissonProjection(Epetra_CrsMatrix &A, Epetra_Vector &n) 
      : _A(A), _n(n) { } 
    virtual~PoissonProjection() { }

    int SetUseTranspose(bool trans);
    int Apply(const Epetra_MultiVector &x, Epetra_MultiVector &y) const;
    int ApplyInverse(const Epetra_MultiVector &x, Epetra_MultiVector &y) const;
    double NormInf() const;

    const char* Label()  const;
    bool  UseTranspose() const;
    bool  HasNormInf()   const;

    const Epetra_Comm& Comm()              const;
    const Epetra_Map&  OperatorDomainMap() const;
    const Epetra_Map&  OperatorRangeMap()  const;
    
  private:
    Epetra_CrsMatrix &_A;
    Epetra_Vector &_n;
  };

  inline int 
  PoissonProjection::SetUseTranspose(bool trans) { 
    return 1;
  }
  inline int 
  PoissonProjection::Apply(const Epetra_MultiVector &x, 
                           Epetra_MultiVector &y) const {
    int r_val = _A.Apply(x, y);
    
    double val = 0.0;
    y.Dot(_n, &val);
    y.Update(-val, _n, 1.0);
    
    return r_val;
  }
  inline int 
  PoissonProjection::ApplyInverse(const Epetra_MultiVector &x, 
                                  Epetra_MultiVector &y) const { 
    return 1;
  }
  inline double 
  PoissonProjection::NormInf() const { 
    return 1.0; 
  }
  inline const char* 
  PoissonProjection::Label()  const { 
    return "PoissonProjection"; 
  }
  inline bool 
  PoissonProjection::UseTranspose() const { 
    return false;
  }
  inline bool 
  PoissonProjection::HasNormInf() const { 
    return false; 
  }

  inline const Epetra_Comm& 
  PoissonProjection::Comm() const { 
    return _A.Comm();
  }
  inline const Epetra_Map& 
  PoissonProjection::OperatorDomainMap() const { 
    return _A.OperatorDomainMap();
  }
  inline const Epetra_Map& 
  PoissonProjection::OperatorRangeMap() const { 
    return _A.OperatorRangeMap(); 
  }


}

#endif

