#pragma once
#ifndef __SOLVER_NOX_H__
#define __SOLVER_NOX_H__


// basic headers

// NOX basic headers
#include "NOX.H"
#include "NOX_Epetra.H"

// Epetra basic headers
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"

// Epetra2NOX base clasees (contains pure virtual functions)
#include "NOX_Epetra_Interface_Required.H" 
#include "NOX_Epetra_Interface_Jacobian.H" 
#include "NOX_Epetra_Interface_Preconditioner.H" 


namespace LAMMPS_NS {

  template<class Problem>
  class SolverNOX {

  public:

    enum SolutionInitType { Random, Zero, Value };
    enum JacobianType { Analytic, MatrixFree, FiniteDifference };

  public:

    SolverNOX(Problem &problem, MPI_Comm &comm);
    SolverNOX(int jac, Problem &problem, MPI_Comm &comm);
    virtual~SolverNOX();

    // solver setting
    int createLinearMap(int num_global_nodes, int index_base);
    int createNodalMap(int num_local_nodes, int *gID);
    int createSolutionVector(double *x);

    void setNodalMap(Epetra_Map *map);
    void setJacobianType(int jacobian_type);
    void setJacobianMatrix(Epetra_RowMatrix *J);
    void setParameters(Teuchos::ParameterList *param = NULL);
    void setConvergenceTests(NOX::StatusTest::Combo *test = NULL);
    void setInitialSolution(SolutionInitType init, double val = 0.0);

    Teuchos::RCP<Epetra_Map> getNodalMap() const;

    // non-linear solve with given linear solver from a derived class
    // AztecOO or Belos
    virtual int solveProblem(const char *name) { return 0; };

    // derived class needs to load default linear solver setting 
    virtual void addLinearSolverParamUnderNewtonParam() = 0; 

  protected:

    // base communicator for Epetra
    Epetra_MpiComm _comm;
    int _pid;

    // map is either created or given; 
    // when they are set from outside, rcp does not have the ownership 
    Teuchos::RCP<Epetra_Map> _map;
    Teuchos::RCP<Problem> _problem;

    // Analytic, MatrixFree, FiniteDifference
    int _jacobian_type;
    
    // main parameter and sublists
    Teuchos::RCP<Teuchos::ParameterList> _param;
    
    // these are delivered from problem and set
    Teuchos::RCP<Epetra_Operator> _J;  

    // abstract interface to solution vector used in NOX
    Teuchos::RCP<NOX::Epetra::Vector> _x; 

    // composable stopping criteria
    Teuchos::RCP<NOX::StatusTest::Combo> _test;

    // nox print utility
    Teuchos::RCP<NOX::Utils> _print;

  };

}

#endif

