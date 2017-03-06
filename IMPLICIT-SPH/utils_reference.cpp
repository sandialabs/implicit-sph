#include <iostream>
#include <bitset>
#include "string.h"

#include "utils.h"
#include "utils_reference.h"

using namespace std;
using namespace LAMMPS_NS;

#define A_(i,j) VIEW2(A, dim, i, j)
#define B_(i,j) VIEW2(B, dim, i, j)

#include "Teuchos_LAPACK.hpp"
static Teuchos::LAPACK<int, double> lapack;

UtilsReference::UtilsReference() { }
UtilsReference::~UtilsReference() { }

double UtilsReference::factorial(const unsigned int n) {
  double r_val = 1.0;
  for (int i=2;i<=n;++i)
    r_val *= i;

  return r_val;
}

void UtilsReference::factorial(const unsigned int n, double *fact) {
  IS_VALID(fact!=NULL, 0, ">> Input array (fact) is NULL");

  double tmp = 1.0;

  fact[0] = 1.0;
  for (int i=1;i<=n;++i) {
    tmp *= i;
    fact[i] = tmp;
  }
}

void UtilsReference::symmetrize(const UploType uplo,
                                const unsigned int m,
                                double *A) {
  IS_VALID(A!=NULL, 0, ">> A is NULL");
  IS_VALID(uplo == UpperTriangular, 0, ">> UploType::UpperTriangular is only supported");

  for (int j=0;j<m;++j)
    for (int i=0;i<j;++i)
      VIEW2(A, m, j, i) = VIEW2(A, m, i, j);
}

void UtilsReference::printDenseMatrix(const char *name,
                                      const unsigned int m,
                                      const unsigned int n,
                                      const double *A) {
  cout << name << scientific << endl;
  for (int i=0;i<m;++i) {
    for (int j=0;j<n;++j)
      cout << VIEW2(A, m, i, j) << "  ";
    cout << endl;
  }
}

void UtilsReference::copyDenseMatrix(const unsigned int m,
                                     const unsigned int n,
                                     double *A,
                                     const bool is_transpose,
                                     double *B) {

  if (is_transpose)
    for (int j=0;j<n;++j)
      for (int i=0;i<m;++i)
        VIEW2(B, n, j, i) = VIEW2(A, m, i, j);
  else
    memcpy(B, A, sizeof(double)*m*n);
}

void UtilsReference::transposeDenseMatrix(const unsigned int m,
                                          const unsigned int n,
                                          double *A) {
  int size = (m*n - 1);

  vector<bool> b;
  b.assign(size + 1, false);
  b[0] = b[size] = 1;

  int i = 0;
  while (++i < size) {
    if (b[i])
      continue;

    int j = i;
    do {
      j = (j == size ? size : (n*j)%size);
      swap(A[j], A[i]);
      b[j] = true;
    } while (j != i);
  }
}


#define A_(i,j) VIEW2(A, dim, i, j)
#define B_(i,j) VIEW2(B, dim, i, j)

double UtilsReference::dotVectors(const unsigned int dim,
                                  const double *x, const double *A, const double *y) {
  QUICK_RETURN(dim == 0, 0.0);

  IS_VALID(x!=NULL, 0, ">> x is NULL");
  IS_VALID(y!=NULL, 0, ">> y is NULL");
  IS_VALID(A!=NULL, 0, ">> y is NULL");

  double r_val = 0.0;
  for (int j=0;j<dim;++j)
    for (int i=0;i<dim;++i)
      r_val += x[i]*A_(i,j)*y[j];

  return r_val;
}


double UtilsReference::dotVectors(const unsigned int dim, const double *x, const double *y) {
  QUICK_RETURN(dim == 0, 0.0);

  IS_VALID(x!=NULL, 0, ">> x is NULL");
  IS_VALID(y!=NULL, 0, ">> y is NULL");

  double r_val = 0.0;
  for (int i=0;i<dim;++i)
    r_val += x[i]*y[i];

  return r_val;
}

double UtilsReference::dotVectors(const unsigned int dim,
                                  const double *x, const int incx,
                                  const double *y, const int incy) {
  QUICK_RETURN(dim == 0, 0.0);

  IS_VALID(x!=NULL, 0, ">> x is NULL");
  IS_VALID(y!=NULL, 0, ">> y is NULL");

  double r_val = 0.0;
  for (int i=0;i<dim;++i)
    r_val += x[i*incx]*y[i*incy];

  return r_val;
}

double UtilsReference::computeDetDenseMatrix(const unsigned int dim, const double *A) {
  double r_val = 0.0;

  QUICK_RETURN(dim == 0, 0.0);
  IS_VALID(A!=NULL, 0, ">> A is NULL");

  switch (dim) {
  case 2:
    r_val = (A_(0,0)*A_(1,1) - A_(0,1)*A_(1,0));
    break;
  case 3:
    r_val = (A_(0,0) * A_(1,1) * A_(2,2) +
             A_(1,0) * A_(2,1) * A_(0,2) +
             A_(2,0) * A_(0,1) * A_(1,2) -
             A_(2,0) * A_(1,1) * A_(0,2) -
             A_(0,0) * A_(2,1) * A_(1,2) -
             A_(1,0) * A_(0,1) * A_(2,2));
    break;
  default:
    // error at this moment
    break;
  }
  return r_val;
}

int UtilsReference::scaleDenseMatrix(const double alpha,
                                     const UploType uplo,
                                     const unsigned int m,
                                     const unsigned int n,
                                     double *A) {
  QUICK_RETURN(m == 0 || n ==0, LAMMPS_SUCCESS);

  IS_VALID(A!=NULL, 0, ">> A is NULL");

  if (alpha == 0.0)
    memset(A, 0, sizeof(double)*m*n);
  else {
    unsigned int size = n*m;
    for (int i=0;i<size;++i)
      A[i] *= alpha;
  }
  return LAMMPS_SUCCESS;
}

int UtilsReference::updateSymRankDenseMatrix(const double alpha,
                                             const unsigned int m,
                                             double *x,
                                             const UploType uplo,
                                             const unsigned int lda,
                                             double *A) {
  IS_VALID(x!=NULL, 0, ">> x is NULL");
  IS_VALID(A!=NULL, 0, ">> A is NULL");
  IS_VALID(uplo != Full, 0, ">> UploType::Full is not supported");

  switch (uplo) {
  case UpperTriangular: {
    for (int j=0;j<m;++j)
      for (int i=0;i<(j+1);++i) // include diagonal
        VIEW2(A, lda, i, j) += alpha*x[i]*x[j];
    break;
  }
  case LowerTriangular: {
    for (int j=0;j<m;++j)
      for (int i=j;i<m;++i) // include diagonal
        VIEW2(A, lda, i, j) += alpha*x[i]*x[j];
    break;
  }
  }

  return LAMMPS_SUCCESS;
}

int UtilsReference::updateSymRankDenseMatrix(const double alpha,
                                             const unsigned int s,
                                             int *idx, double *x,
                                             const UploType uplo,
                                             const unsigned int lda,
                                             double *A) {
  IS_VALID(idx!=NULL, 0, ">> idx is NULL");
  IS_VALID(x!=NULL, 0, ">> x is NULL");
  IS_VALID(A!=NULL, 0, ">> A is NULL");
  IS_VALID(uplo != Full, 0, ">> UploType::Full is not supported");

  switch (uplo) {
  case UpperTriangular: {
    for (int j=0;j<s;++j)
      for (int i=0;i<s;++i) // include diagonal
        if (idx[i] <= idx[j])
          VIEW2(A, lda, idx[i], idx[j]) += alpha*x[i]*x[j];
    break;
  }
  case LowerTriangular: {
    for (int j=0;j<s;++j)
      for (int i=j;i<s;++i) // include diagonal
        VIEW2(A, lda, i, j) += alpha*x[i]*x[j];
    break;
  }
  }

  return LAMMPS_SUCCESS;
}

int UtilsReference::invertDenseMatrix(const unsigned int dim,
                                      double *A,
                                      double *B,
                                      const unsigned int lw,
                                      double *W) {
  QUICK_RETURN(dim == 0, LAMMPS_SUCCESS);

  IS_VALID(A != NULL, 0, ">> A is NULL");
  IS_VALID(B != NULL, 0, ">> B (inverse of A) is NULL");
  IS_VALID(A != B,    0, ">> A and B are the same pointer");

  int r_val = LAMMPS_SUCCESS;

  switch (dim) {
  case 1: {
    if (*A == 0.0) r_val = LAMMPS_FAILURE;

    (*B) = 1.0 / (*A);
    break;
  }
  case 2: {
    const double val = computeDetDenseMatrix(dim, A);
    if (val == 0.0) r_val = LAMMPS_FAILURE;

    B_(0,0) =   A_(1,1)/val;
    B_(1,1) =   A_(0,0)/val;

    B_(1,0) = - A_(1,0)/val;
    B_(0,1) = - A_(0,1)/val;

    break;
  }
  case 3: {
    const double val = computeDetDenseMatrix(dim, A);
    if (val == 0.0) r_val = LAMMPS_FAILURE;

    double val0, val1, val2;

    val0 =   A_(1,1)*A_(2,2) - A_(2,1)*A_(1,2);
    val1 = - A_(1,0)*A_(2,2) + A_(2,0)*A_(1,2);
    val2 =   A_(1,0)*A_(2,1) - A_(2,0)*A_(1,1);

    B_(0,0) = val0/val;
    B_(1,0) = val1/val;
    B_(2,0) = val2/val;

    val0 =   A_(2,1)*A_(0,2) - A_(0,1)*A_(2,2);
    val1 =   A_(0,0)*A_(2,2) - A_(2,0)*A_(0,2);
    val2 = - A_(0,0)*A_(2,1) + A_(2,0)*A_(0,1);

    B_(0,1) = val0/val;
    B_(1,1) = val1/val;
    B_(2,1) = val2/val;

    val0 =   A_(0,1)*A_(1,2) - A_(1,1)*A_(0,2);
    val1 = - A_(0,0)*A_(1,2) + A_(1,0)*A_(0,2);
    val2 =   A_(0,0)*A_(1,1) - A_(1,0)*A_(0,1);

    B_(0,2) = val0/val;
    B_(1,2) = val1/val;
    B_(2,2) = val2/val;
    break;
  }
  default: {
    // set B unit matrix
    memset(B, 0, sizeof(double)*dim*dim);
    for (int i=0;i<dim;++i)
      B_(i,i) = 1.0;

    // invert
    r_val = solveDenseMatrix(Full, dim, dim, A, dim, dim, B, lw, W);
    break;
  }
  }
  return r_val;
}

int UtilsReference::invertDenseMatrixSVD(const unsigned int dim,
                                         double *A,
                                         double *B,
                                         const unsigned int lw,
                                         double *W) {
  QUICK_RETURN(dim == 0, LAMMPS_SUCCESS);

  IS_VALID(A != NULL, 0, ">> A is NULL");
  IS_VALID(B != NULL, 0, ">> B (inverse of A) is NULL");
  IS_VALID(A != B,    0, ">> A and B are the same pointer");

  int r_val = LAMMPS_SUCCESS;

  // set B unit matrix
  memset(B, 0, sizeof(double)*dim*dim);
  for (int i=0;i<dim;++i)
    B_(i,i) = 1.0;

  // invert
  r_val = solveDenseMatrixSVD(Full, dim, dim, A, dim, dim, B, lw, W);

  return r_val;
}

double UtilsReference::checkInversion(const unsigned int dim,
                                      double *A,
                                      double *Ainv) {
  double r_val = 0.0;
  for (int i=0;i<dim;++i) {
    for (int j=0;j<dim;++j) {
      double val = 0;
      for (int k=0;k<dim;++k)
        val += VIEW2(A, dim, i, k)*VIEW2(Ainv, dim, k, j);

      double diff = abs(val) - (i==j);
      r_val += abs(diff);
    }
  }

  return r_val;
}


int UtilsReference::solveDenseMatrix(const UploType uplo,
                                     const unsigned int mA,
                                     const unsigned int nA,
                                     double *A,
                                     const unsigned int mB,
                                     const unsigned int nB,
                                     double *B,
                                     const unsigned int lw,
                                     double *W) {
  QUICK_RETURN(mA == 0 || nA == 0, LAMMPS_SUCCESS);
  QUICK_RETURN(mB == 0 || nB == 0, LAMMPS_SUCCESS);

  IS_VALID(A!=NULL, 0, ">> A is NULL");
  IS_VALID(B!=NULL, 0, ">> B (inverse of A) is NULL");
  IS_VALID(A!=B,    0, ">> A and B are the same pointer");

  IS_VALID(uplo==Full, 0, ">> UploType::Full is only supported now");

  int r_val = LAMMPS_SUCCESS;
  int val = mA - nA;
  int branch = ((0 < val) - (val < 0));

  switch (branch) {
  case 0: {
    IS_VALID(mA == mB,
             0, ">> A and B has a different numer of rows:" << mA << ", " << mB);

    if (uplo == Full) {
      // square - LU with partial pivoting
      IS_VALID(lw>=mA, 0, ">> Work array size is too small");
      int *ipiv = (int*)W;

      lapack.GESV(mA, nB, A, mA, ipiv, B, mA, &r_val);
      if (r_val != LAMMPS_SUCCESS) {
        cout << ">> DGESV failed: info " << r_val << endl;
        r_val = LAMMPS_FAILURE;
      }
    } else {
      // spd - Cholesky
      char ch_uplo = (uplo == UpperTriangular ? 'U' : 'L');

      lapack.POSV(ch_uplo, mA, nB, A, mA, B, mB, &r_val);

      if (r_val != LAMMPS_SUCCESS) {
        cout << ">> DPOSV failed: info " << r_val << endl;
        r_val = LAMMPS_FAILURE;
      }
    }
    break;
  }
  case 1: {
    // over-determined - QR
    cout << ">> QR solver is not yet implemented" << endl;
    r_val = LAMMPS_FAILURE;
    break;
  }
  case -1: {
    // under-determined - LQ
    cout << ">> LQ solver is not yet implemented" << endl;
    r_val = LAMMPS_FAILURE;
    break;
  }
  }
  return r_val;
}

int UtilsReference::solveDenseMatrixSVD(const UploType uplo,
                                        const unsigned int mA,
                                        const unsigned int nA,
                                        double *A,
                                        const unsigned int mB,
                                        const unsigned int nB,
                                        double *B,
                                        const unsigned int lw,
                                        double *W) {
  QUICK_RETURN(mA == 0 || nA == 0, LAMMPS_SUCCESS);
  QUICK_RETURN(mB == 0 || nB == 0, LAMMPS_SUCCESS);

  IS_VALID(A!=NULL, 0, ">> A is NULL");
  IS_VALID(B!=NULL, 0, ">> B (inverse of A) is NULL");
  IS_VALID(A!=B,    0, ">> A and B are the same pointer");

  IS_VALID(uplo==Full, 0, ">> UploType::Full is only allowed");

  int r_val = LAMMPS_SUCCESS, shift;

  const int lw_min = min(mA, nA)+ 3*min(mA, nA) + max(max(2*min(mA, nA), max(mA, nA)), nB);
  IS_VALID(lw>=lw_min, 0, ">> Work array size is too small");

  double *S = W; shift = min(mA, nA);

  const unsigned int lww = lw - shift;
  W = &S[shift];

  const double rcond = 1.0e-12;
  int rank = 0;

  lapack.GELSS(mA, nA, nB, A, mA, B, mB, S, rcond, &rank, W, lww, &r_val);
  if (r_val != LAMMPS_SUCCESS) {
    cout << ">> DGELSS failed: info " << r_val << endl;
    r_val = LAMMPS_FAILURE;
  }

  return r_val;
}

// Build a rotation matrix Q such that   Q * n = [1 0]'  in 2d and Q * n = [1 0 0]' in 3d.
// The vector n is expected to be unitary (i.e. norm(n) == 1).
int UtilsReference::computeRotationMatrix(const unsigned int dim,
                                          const double *n,
                                          double *Q) {
  if (dim == 3) {
    double n0(n[0]), n1(n[1]), n2(n[2]);
    if (n2*n2 < 0.5) { //needed to make sure that c is not close to 0.
      double c = sqrt(1-n2*n2);
      Q[0] = n0;       Q[3] = n1;       Q[6] = n2;
      Q[1] = -n1/c;    Q[4] = n0/c;     Q[7] = 0;
      Q[2] = -n0*n2/c; Q[5] = -n1*n2/c; Q[8] = c;
    } else {   // n1*n1 <= 0.5
      double c = sqrt(1-n1*n1);
      Q[0] = n0;       Q[3] = n1;       Q[6] = n2;
      Q[1] = -n2/c;    Q[4] = 0;        Q[7] = n0/c;
      Q[2] = -n0*n1/c; Q[5] = c;        Q[8] = -n1*n2/c;
    }
  } else if(dim == 2) {
    Q[0] = n[0];  Q[2] = n[1];
    Q[1] = -n[1]; Q[3] = n[0];
  } else if(dim == 1) //trivial case
    Q[0] = n[0];
  else
    return LAMMPS_FAILURE;

  return LAMMPS_SUCCESS;
}
