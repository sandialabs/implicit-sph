#pragma once
#ifndef __SOLVER_NOX_AZTECOO_H__
#define __SOLVER_NOX_AZTECOO_H__


#include "solver_nox.h"
#include "solver_nox_impl.h"

#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

namespace LAMMPS_NS {

  template<class Problem>
  class SolverNOX_AztecOO : public SolverNOX<Problem> {
  public:


    SolverNOX_AztecOO(Problem &problem, MPI_Comm &comm);
    SolverNOX_AztecOO(int jac, Problem &problem, MPI_Comm &comm);
    
    int solveProblem(const char *name);
    void addLinearSolverParamUnderNewtonParam();

  };

  template<class Problem> inline 
  SolverNOX_AztecOO<Problem>::SolverNOX_AztecOO(Problem &problem, MPI_Comm &comm) 
    : SolverNOX<Problem>(problem, comm) { }
  template<class Problem> inline
  SolverNOX_AztecOO<Problem>::SolverNOX_AztecOO(int jac, Problem &problem, MPI_Comm &comm) 
    : SolverNOX<Problem>(jac, problem, comm) { }
  
  template<class Problem> int
  SolverNOX_AztecOO<Problem>::solveProblem(const char *name) {
    auto& print_param      = this->_param->sublist("Printing");
    auto& lin_solver_param = this->_param->sublist("Direction").sublist("Newton").sublist("Linear Solver");
    auto& problem          = this->_problem;
    auto& print            = this->_print;

    if (name != NULL)
      print->out() << ">> SolverNOX<Problem>::Label - " << name << endl;

    Teuchos::RCP<NOX::Epetra::Interface::Required>      required = problem;
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> precond = problem;
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linear_system;

    auto& J = this->_J;
    auto& x = this->_x;

    // over-ride jacobian_type if the operator is null
    auto& jacobian_type = this->_jacobian_type;
    if (jacobian_type == SolverNOX<Problem>::Analytic && J == Teuchos::null) {
      print->out() << ERROR(">> SolverNOX::Error - User must provide Jacobian") << endl;
      throw -1;
    }    

    switch (jacobian_type) {
    case SolverNOX<Problem>::Analytic: {
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> jacobian = problem;
      linear_system = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(print_param,
                                                                        lin_solver_param,
                                                                        required,
                                                                        jacobian,
                                                                        J,
                                                                        *x));
      break;
    }
    case SolverNOX<Problem>::MatrixFree: {
      auto op_matrix_free = Teuchos::rcp(new NOX::Epetra::MatrixFree(print_param, problem, *x));
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> jacobian = op_matrix_free;
      linear_system = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(print_param,
                                                                        lin_solver_param,
                                                                        required,
                                                                        jacobian,
                                                                        op_matrix_free,
                                                                        *x));
      break;
    }
    case SolverNOX<Problem>::FiniteDifference: {
      auto op_finite_difference = Teuchos::rcp(new NOX::Epetra::FiniteDifference(print_param, problem, *x));
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> jacobian = op_finite_difference;
      linear_system = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(print_param,
                                                                        lin_solver_param,
                                                                        required,
                                                                        jacobian,
                                                                        op_finite_difference,
                                                                        *x));
      break;
    }
    }
    auto group = Teuchos::rcp(new NOX::Epetra::Group(print_param,
                                                     required,
                                                     *x,
                                                     linear_system));
    
    // set convergence tests
    this->setConvergenceTests(this->_test.getRawPtr());

    // put all together
    auto solver = NOX::Solver::buildSolver(group, this->_test, this->_param);

    if (solver->solve() == NOX::StatusTest::Converged)
      print->out() << ">> NOX::Utils - Passed!" << endl;
    else
      print->out() << ">> NOX::Utils - The solver failed to converge!" << endl;
    
    auto& final_group = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    *x = dynamic_cast<const NOX::Epetra::Vector&>(final_group.getX());
    
    return LAMMPS_SUCCESS;
  }
  template<class Problem> inline void
  SolverNOX_AztecOO<Problem>::addLinearSolverParamUnderNewtonParam() {
    auto& linear = this->_param->sublist("Direction").sublist("Newton").sublist("Linear Solver");
    
    linear.set("Aztec Solver", "GMRES");
    linear.set("Max Iterations", 800);
    linear.set("Tolerance", 1e-4);
    linear.set("Preconditioner", "None");
    //linear.set("Preconditioner", "Ifpack");                                                                
    linear.set("Max Age Of Prec", 5);
  }
  
}

#endif

