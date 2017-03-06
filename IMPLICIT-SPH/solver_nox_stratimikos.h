#pragma once
#ifndef __SOLVER_NOX_STRATIMIKOS_H__
#define __SOLVER_NOX_STRATIMIKOS_H__


#include "solver_nox.h"
#include "solver_nox_impl.h"

#include "Epetra_LinearProblem.h"
#include "NOX_Epetra_LinearSystem_Stratimikos.H"

namespace LAMMPS_NS {

  template<class Problem>
  class SolverNOX_Stratimikos : public SolverNOX<Problem> {
  public:

    SolverNOX_Stratimikos(Problem &problem, MPI_Comm &comm);
    SolverNOX_Stratimikos(int jac, Problem &problem, MPI_Comm &comm);
    
    int solveProblem(const char *name);
    void addLinearSolverParamUnderNewtonParam();
    
  };
  
  // Stratimikos always require a user-provided Jacobian operator
  template<class Problem> inline
  SolverNOX_Stratimikos<Problem>::SolverNOX_Stratimikos(Problem &problem, MPI_Comm &comm)
    : SolverNOX<Problem>(SolverNOX<Problem>::Analytic, problem, comm) { } 
  template<class Problem> inline 
  SolverNOX_Stratimikos<Problem>::SolverNOX_Stratimikos(int jac, Problem &problem, MPI_Comm &comm)
    : SolverNOX<Problem>(SolverNOX<Problem>::Analytic, problem, comm) { } 

  template<class Problem> inline int
  SolverNOX_Stratimikos<Problem>::solveProblem(const char *name) {
    auto& print_param      = this->_param->sublist("Printing");
    auto& newton_param     = this->_param->sublist("Direction").sublist("Newton");
    auto& lin_solver_param = newton_param.sublist("Stratimikos Linear Solver");
    auto& problem          = this->_problem;
    auto& print            = this->_print;

    if (name != NULL)
      print->out() << ">> SolverNOX<Problem>::Label - " << name << endl;

    Teuchos::RCP<NOX::Epetra::Interface::Required>      required = problem;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian>      jacobian = problem;
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> precond = problem;
    
    if (this->_J == Teuchos::null) {
      print->out() << ">> SolverNOX::Error - User must provide Jacobian" << endl;
      throw -1;
    }
    auto& J = this->_J;
    auto& x = this->_x;
    
    auto linear_system = Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(print_param,
                                                                               lin_solver_param,
                                                                               jacobian,
                                                                               J,
                                                                               *x));
    
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
  SolverNOX_Stratimikos<Problem>::addLinearSolverParamUnderNewtonParam() {
    auto& newton = this->_param->sublist("Direction").sublist("Newton");
    auto& linear = newton.sublist("Stratimikos Linear Solver");
    {
      auto& nox = linear.sublist("NOX Stratimikos Options");
      
      // following options are avail for preconditioner 
      nox.set("Preconditioner Reuse Policy","Reuse");
      // nox.set("Preconditioner Reuse Policy","Rebuild");
      // nox.set("Preconditioner Reuse Policy","Recompute");

      nox.set("Max Age Of Prec", 10);
      
      // Inexact Newton must be set in a second sublist when using
      // Stratimikos: This code snippet sets it automatically
      bool inexact = (newton.get("Forcing Term Method", "Constant") != "Constant");
      nox.set("Use Linear Solve Tolerance From NOX", inexact);
    }
    {
      auto& strat = linear.sublist("Stratimikos");
      
      strat.set("Linear Solver Type", "Belos");

      //strat.set("Preconditioner Type", "Ifpack");
      strat.set("Preconditioner Type", "ML");

      {
        auto& belos = strat.sublist("Linear Solver Types").sublist("Belos");
        belos.set("Solver Type","Block GMRES");
        {
          auto& belos_type =  belos.sublist("Solver Types").sublist("Block GMRES");
          belos_type.set("Convergence Tolerance", 1.0e-6);
          belos_type.set("Maximum Iterations", 80);
          belos_type.set("Verbosity",0);
        }
      }
    }
  }
  
}

#endif

