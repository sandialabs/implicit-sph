#pragma once
#ifndef __SOLVER_LIN_BELOS_H__
#define __SOLVER_LIN_BELOS_H__


#include "utils.h"
#include "solver_lin.h"

#include "BelosEpetraAdapter.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosGCRODRSolMgr.hpp"

#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Teuchos_Array.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "BelosThyraAdapter.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"

#include "precond.h"
#include "precond_ml.h"

namespace LAMMPS_NS {
  
  
  class SolverLin_Belos : public SolverLin {
  public:
    
    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator    OP;
    
    typedef Thyra::MultiVectorBase<double> MVt;
    typedef Thyra::LinearOpBase<double>    OPt;
    
    SolverLin_Belos(MPI_Comm &comm) : SolverLin(comm) { }
    int solveProblem(PrecondWrapper *prec = NULL, const char *name = NULL);
    int solveBlockProblem(PrecondWrapper *prec = NULL, const char *name = NULL);
    
    void setParameters(Teuchos::ParameterList *param = NULL);
  };

  using namespace std;

  inline int 
  SolverLin_Belos::solveBlockProblem(PrecondWrapper *prec, const char* name) {
    // if (_comm.MyPID() == 0 && name != NULL)
    //   cout << ">> Belos(Block)::Label - "  << name << endl;

    IS_VALID(_b->NumVectors() == _dim, _comm.MyPID(), 
             ">> SolverLin_Belos::solveBlockProblem, dimension of rhs does not match to the block matrix");
    IS_VALID(!_is_singular, _comm.MyPID(), 
             ">> SolverLin_Belos::solveBlockProblem does not support singular problems");
    
    // if null, load a default parameterlist
    setParameters(_param.get());
   
    // diagonal preconditioner
    Teuchos::RCP<Thyra::LinearOpBase<double> > prec_belos;  
    if (prec != NULL) {
      prec->create(_dim);
      if (!prec->getBlockPrecondOperator()) 
        throw runtime_error("SolverLin_Belos::solveProblem getBlockPrecondOperator failed");
      
      prec_belos = Teuchos::rcp(prec->getBlockPrecondOperator(), false);
    }
    
    // Build the blocked vectors (note: these are views) & wrap to thyra
    Teuchos::RCP<Thyra::ProductMultiVectorBase<double> > bthyra, xthyra;
    createBlockVector(_x, xthyra);
    createBlockVector(_b, bthyra);
    
    // put together Ax = b
    Teuchos::RCP<Belos::LinearProblem<double, MVt, OPt> > problem;
    problem = Teuchos::rcp(new Belos::LinearProblem<double,MVt,OPt>(_Athyra, xthyra, bthyra));
    
    //problem->setLeftPrec(prec_belos);
    problem->setRightPrec(prec_belos);
    problem->setProblem();
    
    // create a solver object
    Teuchos::RCP< Belos::SolverManager<double,MVt,OPt> > solver;
    {
      auto solver_type = _param->get("Solver Type","Block GMRES");
      if      (solver_type == "Block GMRES") 
        solver = Teuchos::rcp(new Belos::BlockGmresSolMgr<double,MVt,OPt>(problem, _param));
      else if (solver_type == "Block CG") 
        solver = Teuchos::rcp(new Belos::BlockCGSolMgr<double,MVt,OPt>(problem, _param));
    }
    
    auto status = solver->solve();
    
    if (prec != NULL)
      prec->free();
    
    // solve the problem
    if (status == Belos::Converged) {
      // if (_comm.MyPID() == 0) 
      //   cout << ">> Belos::Status - Passed! " << (name == NULL ? " " : name) << endl;
    } else {
      if (_comm.MyPID() == 0)
        cout << ">> Belos::Status - Failed to converge! " << (name == NULL ? " " : name) << endl;;

      // Print final residual information 
      MV resid(_A->RowMap(), true);
      _A->Apply(*_x, resid);
      resid.Update(1.0, *_b, -1.0);
    
      double r_norm = 0.0;
      resid.Norm2(&r_norm);
    
      double b_norm = 0.0;
      _b->Norm2(&b_norm);
    
      if (_comm.MyPID() == 0)
        printf(">> Belos:: ||r|| / ||b|| = %6.4e\n", r_norm/ b_norm);
    }
  
    return LAMMPS_SUCCESS;
  }

  inline int
  SolverLin_Belos::solveProblem(PrecondWrapper *prec, const char *name) {
    // if (_comm.MyPID() == 0 && name != NULL)
    //   cout << ">> Belos::Label - "  << name << endl;;

    // if null, load a default parameterlist
    setParameters(_param.get());

    if (_is_singular) {
      this->createNullVector();

      double val = 0.0;
      _b->Dot(*_n, &val);
      _b->Update(-val, *_n, 1.0);
    }

    // preconditioner
    Teuchos::RCP<Belos::EpetraPrecOp> prec_belos;
    if (prec != NULL) {
      // for now let us consider poisson problem only
      if (_is_singular) 
        prec->setNullVector(_n->Values());

      prec->create();
      auto op = Teuchos::rcp(prec->getPrecondOperator(), false);
      prec_belos = Teuchos::rcp(new Belos::EpetraPrecOp(op));
    }

    // put together Ax = b
    // if singular, create null vector and set it with 1/norm
    // for poisson problem on pressure field only now
    Teuchos::RCP<Belos::LinearProblem<double, MV, OP> > problem;
    if (_is_singular) {
      Teuchos::RCP<PoissonProjection> P = Teuchos::rcp(new PoissonProjection(*_A, *_n));
      problem = Teuchos::rcp(new Belos::LinearProblem<double,MV,OP>(P,_x,_b));
    } else {    
      problem = Teuchos::rcp(new Belos::LinearProblem<double,MV,OP>(_A,_x,_b));
    }
    //problem->setLeftPrec(prec_belos);
    problem->setRightPrec(prec_belos);
    problem->setProblem();

    // create a solver object
    Teuchos::RCP< Belos::SolverManager<double,MV,OP> > solver;
    {
      auto solver_type = _param->get("Solver Type","Block GMRES");
      if      (solver_type == "Block GMRES") 
        solver = Teuchos::rcp(new Belos::BlockGmresSolMgr<double,MV,OP>(problem, _param));
      else if (solver_type == "Recycling GMRES") 
        solver = Teuchos::rcp(new Belos::GCRODRSolMgr<double,MV,OP>(problem, _param));
      else if (solver_type == "Block CG") 
        solver = Teuchos::rcp(new Belos::BlockCGSolMgr<double,MV,OP>(problem, _param));
    }

    auto status = solver->solve();

    if (prec != NULL) {
      if (_is_singular)
        prec->setNullVector(NULL);

      prec->free();
    }
    
    // solve the problem
    if (status == Belos::Converged) {
      // if (_comm.MyPID() == 0) 
      //   cout << ">> Belos::Status - Passed! " << (name == NULL ? " " : name) << endl;
    } else {
      if (_comm.MyPID() == 0)
        cout << ">> Belos::Status - Failed to converge! " << (name == NULL ? " " : name) << endl;;

      MV resid(_A->RowMap(), true);
      _A->Apply(*_x, resid);
      resid.Update(1.0, *_b, -1.0);

      double r_norm = 0.0;
      resid.Norm2(&r_norm);

      double b_norm = 0.0;
      _b->Norm2(&b_norm);

      if (_comm.MyPID() == 0)
        printf(">> Belos:: ||r|| / ||b|| = %6.4e\n", r_norm/ b_norm);
    }

    if (_is_singular) {
      double val = 0.0;
      _x->Dot(*_n, &val);
      _x->Update(-val, *_n, 1.0);
    }

    return LAMMPS_SUCCESS;
  }

  inline void 
  SolverLin_Belos::setParameters(Teuchos::ParameterList *param) {
    if (param == NULL) {

      // create top-level parameter
      _param = Teuchos::rcp(new Teuchos::ParameterList);  
      
      //_param->set("Flexible Gmres", true);           // Flexible Gmres
      _param->set("Num Blocks", 50);                 // Maximum number of blocks in Krylov factorization
      _param->set("Block Size", 1);                  // Blocksize to be used by iterative solver
      _param->set("Maximum Iterations", 500);        // Maximum number of iterations allowed
      _param->set("Maximum Restarts", 15 );          // Maximum number of restarts allowed                
      _param->set("Convergence Tolerance", 1.0e-6);  // Relative convergence tolerance requested
      _param->set("Orthogonalization", "IMGS");    // Orthogonalization method
      //_param->set("Orthogonalization", "DGKS");      // Orthogonalization method

      // Solver choice : somehow CG does not work on me... 
      _param->set("Solver Type","Block CG"); 
      //_param->set("Solver Type","Block GMRES");
      //_param->set("Solver Type","Recycling GMRES"); 
      _param->set("Num Recycled Blocks", 30);
     
      // output 
      _param->set ("Output Frequency", 10);
      _param->set ("Output Style", 1);
      _param->set ("Verbosity", 33);

      // _param->set("Verbosity", 
      //             Belos::Errors + 
      //             Belos::Warnings + 
      //             Belos::FinalSummary);
      
      // _param->set("Verbosity", 
      //             Belos::Errors +
      //             Belos::Warnings +                                           
      //             Belos::TimingDetails + 
      //             Belos::FinalSummary + 
      //             Belos::StatusTestDetails);
      
      _param->set("Output Frequency", 5);                                                       
    } else {
      if (_param.get() != param) 
        _param = Teuchos::rcp(param, false); 
    }
  }
  
}

#endif

