// user header
#include "precond.h"
#include "solver_lin.h"
#include "utils.h"

#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Teuchos_Array.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

// these are developer headers
//#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
// #include "Thyra_DefaultBlockedLinearOp_decl.hpp"
// #include "Thyra_DefaultBlockedLinearOp_def.hpp"

namespace LAMMPS_NS {

  using namespace std;

  // ------------------------------------------------------------------------------
  SolverLin::SolverLin(MPI_Comm &comm) 
    : _comm(comm),_dim(0),_is_singular(false) { }
  SolverLin::~SolverLin() { }

  // ------------------------------------------------------------------------------
  int SolverLin::createLinearMap(int num_global_nodes, int index_base) {
    _map = Teuchos::rcp(new Epetra_Map(num_global_nodes, index_base, _comm));
    return LAMMPS_SUCCESS;
  }
  int SolverLin::createNodalMap(int num_local_nodes, int *gID) {
    auto idx = Teuchos::rcp(new Epetra_IntSerialDenseVector(Copy, gID, num_local_nodes)); 
    _map = Teuchos::rcp(new Epetra_Map(-1, num_local_nodes, idx->Values(), 1, _comm));
    
    return LAMMPS_SUCCESS;
  }
  int SolverLin::createLoadMultiVector(double *b, int lda, int num_vectors) {
    if (b == NULL)
      _b = Teuchos::rcp(new Epetra_MultiVector(*_map, num_vectors));
    else
      _b = Teuchos::rcp(new Epetra_MultiVector(View, *_map, b, lda, num_vectors));
    return LAMMPS_SUCCESS;
  }
  int SolverLin::createSolutionMultiVector(double *x, int lda, int num_vectors) {
    if (x == NULL)
      _x = Teuchos::rcp(new Epetra_MultiVector(*_map, num_vectors));
    else
      _x = Teuchos::rcp(new Epetra_MultiVector(View, *_map, x, lda, num_vectors));
    return LAMMPS_SUCCESS;
  }
  int SolverLin::createNullVector() {
    _n = Teuchos::rcp(new Epetra_Vector(*_map, true));
    
    if (_null_mask == Teuchos::null) {
      _n->PutScalar(1.0);
    } else {
      const int length = min(_null_mask->Length(), _n->MyLength());
      const int *mask = _null_mask->Values();
      double *value = _n->Values();
      for (int i=0;i<length;++i)
        value[i] = mask[i];
    }

    double norm;
    _n->Norm2(&norm);
    _n->Scale(1.0/norm);

    return LAMMPS_SUCCESS;
  }
  int SolverLin::createBlockMatrix(const int dim, const char *name) {
    _Athyra = Thyra::defaultBlockedLinearOp<double>();
    _Athyra->setObjectLabel(name);
    _dim = dim;

    return LAMMPS_SUCCESS;
  }
  int SolverLin::freeBlockMatrix() {
    _Athyra = Teuchos::null;
    _dim = 0;

    return LAMMPS_SUCCESS;
  }

  // this is internal function
  int SolverLin::createBlockVector(const Teuchos::RCP<Epetra_MultiVector> v,
                                   Teuchos::RCP<Thyra::ProductMultiVectorBase<double> > &vthyra) {
    
    Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<double> > > vv(_dim);
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > vs(_dim);
    for (int k=0;k<_dim;++k) {
      vv[k] = Thyra::create_Vector(Teuchos::rcp(new Epetra_Vector(View, *v, k)),
                                   Thyra::create_VectorSpace(_map));
      vs[k] = vv[k]->range();
    }

    vthyra = Thyra::defaultProductVector<double>(Thyra::productVectorSpace<double>(vs), vv);

    return LAMMPS_SUCCESS;
  }
  // ------------------------------------------------------------------------------
  void SolverLin::setNullVectorMask(Epetra_IntSerialDenseVector *mask) {
    if (mask != NULL) 
      _null_mask = Teuchos::rcp(mask, false);
  }
  void SolverLin::setMatrixIsBlocked(const bool is_blocked) {
    _is_blocked = is_blocked;
  }
  void SolverLin::setMatrixIsSingular(const bool is_singular) {
    _is_singular = is_singular;
  }
  void SolverLin::setNodalMap(Epetra_Map *map) {
    if (map != NULL) 
      _map = Teuchos::rcp(map, false);
  }
  void SolverLin::setMatrix(Epetra_CrsMatrix *A) {
    if (A != NULL)
      _A = Teuchos::rcp(A, false);
  }
  void SolverLin::setBlockBegin() {
    _Athyra->beginBlockFill(_dim, _dim);
  }
  void SolverLin::setBlock(const int i, const int j, 
                           const Epetra_CrsMatrix *A) {
    // I do not check whether Athyra is created or not; let us clean later
    if ((A != NULL) && (i < _dim && j < _dim))
      _Athyra->setBlock(i, j, Thyra::epetraLinearOp(Teuchos::rcp(A, false)));
  }
  void SolverLin::setBlockEnd() {
    _Athyra->endBlockFill();
  }
  // ------------------------------------------------------------------------------
  void SolverLin::setInitialSolution(SolutionInitType init, double val) {
    switch (init) {
    case SolutionInitType::Random:   _x->Random();       break;
    case SolutionInitType::Zero  :   _x->PutScalar(0.0); break;
    case SolutionInitType::Value :   _x->PutScalar(val); break;
    }
  }
  // ------------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Map> SolverLin::getNodalMap() const {
    return _map;
  }
  Teuchos::RCP<Epetra_MultiVector> SolverLin::getLoadMultiVector() {
    return _b;
  }
  Teuchos::RCP<Epetra_MultiVector> SolverLin::getSolutionMultiVector() {
    return _x;
  }
  Teuchos::RCP<Epetra_Vector> SolverLin::getNullVector() {
    return _n;
  }
}

