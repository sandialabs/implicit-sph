#pragma once
#ifndef __SOLVER_NOX_IMPL_H__
#define __SOLVER_NOX_IMPL_H__


// basic headers
#define HAVE_MPI
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// implementation specific headers
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_ParameterList.hpp"

// user header
#include "solver_nox.h"
#include "utils.h"

namespace LAMMPS_NS {

  using namespace std;

  // ------------------------------------------------------------------------------
  template<class Problem> inline
  SolverNOX<Problem>::SolverNOX(Problem &problem, MPI_Comm &comm)
    : _comm(comm), _problem(&problem, false), _jacobian_type(JacobianType::MatrixFree) { }
  template<class Problem> inline
  SolverNOX<Problem>::SolverNOX(int jac, 
                                Problem &problem, MPI_Comm &comm) 
    : _comm(comm), _problem(&problem, false), _jacobian_type(jac) { }
  template<class Problem> inline
  SolverNOX<Problem>::~SolverNOX() { }

  // ------------------------------------------------------------------------------
  template<class Problem> inline int
  SolverNOX<Problem>::createLinearMap(int num_global_nodes, int index_base) {
    _map = Teuchos::rcp(new Epetra_Map(num_global_nodes, index_base, _comm));
    return LAMMPS_SUCCESS;
  }
  template<class Problem> inline int
  SolverNOX<Problem>::createNodalMap(int num_local_nodes, int *gID) {
    auto idx = Teuchos::rcp(new Epetra_IntSerialDenseVector(Copy, gID, num_local_nodes)); 
    _map = Teuchos::rcp(new Epetra_Map(-1, num_local_nodes, idx->Values(), 1, _comm));

    return LAMMPS_SUCCESS;
  }
  template<class Problem> inline int  
  SolverNOX<Problem>::createSolutionVector(double *x) {
    _x = Teuchos::rcp(new NOX::Epetra::Vector(Teuchos::rcp(new Epetra_Vector(View, *_map, x)),
                                              NOX::Epetra::Vector::CreateView));
    return LAMMPS_SUCCESS;
  }
  // ------------------------------------------------------------------------------
  template<class Problem> inline void
  SolverNOX<Problem>::setNodalMap(Epetra_Map *map) {
    if (map != NULL) 
      _map = Teuchos::rcp(map, false);
  }
  template<class Problem> inline void
  SolverNOX<Problem>::setJacobianMatrix(Epetra_RowMatrix *J) {
    if (J != NULL)
      _J = Teuchos::rcp(J, false);
  }
  template<class Problem> inline void
  SolverNOX<Problem>::setJacobianType(int jacobian_type) {
    _jacobian_type = jacobian_type;
  }
  template<class Problem> inline void
  SolverNOX<Problem>::setParameters(Teuchos::ParameterList *param) {
    if (param == NULL) {

      // create top level parameter list 
      _param = Teuchos::rcp(new Teuchos::ParameterList);
    
      // set default parameters for NOX
      _param->set("Nonlinear Solver", "Line Search Based");
      // _param->set("Nonlinear Solver", "Trust Region Based");
    
      { // set printing parameters in the "Printing" sublist
        auto& print_param = _param->sublist("Printing");
        print_param.set("Output Precision", 6);
        print_param.set("Output Processor", 0);

        print_param.set("Output Information", 
                        NOX::Utils::OuterIteration + 
                        NOX::Utils::OuterIterationStatusTest + 
                        NOX::Utils::InnerIteration +
                        NOX::Utils::LinearSolverDetails +
                        NOX::Utils::Parameters + 
                        NOX::Utils::Details + 
                        NOX::Utils::Warning +
                        NOX::Utils::Debug + 
                        NOX::Utils::TestDetails +
                        NOX::Utils::Error);
        
        // create printer
        _print = Teuchos::rcp(new NOX::Utils(print_param));
      }
    
      { // set line search parameters
        auto& search_param = _param->sublist("Line Search");
        search_param.set("Method", "Full Step");
      }
    
      { // set direction parameters
        auto& direction_param = _param->sublist("Direction");
        {
          // when it attempts to retrieve newton parameter, it will use "Newton"
          auto& newton_param = direction_param.sublist("Newton");
          newton_param.set("Forcing Term Method", "Constant");
          addLinearSolverParamUnderNewtonParam();
        }
      }
    } else {
      _param = Teuchos::rcp(param, false);
    }
  }
  template<class Problem> inline void
  SolverNOX<Problem>::setConvergenceTests(NOX::StatusTest::Combo *test) {
    if (test == NULL) { // set default stopping criteria
      _test = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
      { 
        auto maxiters  = Teuchos::rcp(new NOX::StatusTest::MaxIters(100));
        auto fv        = Teuchos::rcp(new NOX::StatusTest::FiniteValue);
        auto converged = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

        _test->addStatusTest(fv);
        _test->addStatusTest(maxiters);

        {
          auto abs_res = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
          //auto rel_res = Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), 1.0e-2));
          auto update  = Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
          //auto wrms    = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));

          converged->addStatusTest(abs_res);
          // converged->addStatusTest(rel_res);
          //converged->addStatusTest(wrms);
          converged->addStatusTest(update);
        }
        _test->addStatusTest(converged);
      }
    } 
    else if (_test.get() != test) {
      _test = Teuchos::rcp(test, false);
    }
  }
  template<class Problem> inline void
  SolverNOX<Problem>::setInitialSolution(SolutionInitType init, double val) {
    switch (init) {
    case SolutionInitType::Random:   _x->random();  break;
    case SolutionInitType::Zero  :   _x->init(0.0); break;
    case SolutionInitType::Value :   _x->init(val); break;
    }
  }
  template<class Problem> inline Teuchos::RCP<Epetra_Map>
  SolverNOX<Problem>::getNodalMap() const {
    return _map.get();
  }

}

#endif


