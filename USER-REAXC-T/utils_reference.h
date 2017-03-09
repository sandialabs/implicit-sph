#ifndef __UTILS_REFERENCE_H__
#define __UTILS_REFERENCE_H__

// 1D Quintic spline kernel

namespace LAMMPS_NS {

  class UtilsReference {
  public:

    enum UploType { Full, UpperTriangular, LowerTriangular };

  public:

    UtilsReference();
    virtual~UtilsReference();

    double factorial(const unsigned int n);
    void factorial(const unsigned int n, double *fact);

    // let's not consider lda; assume all matrices are fully packed (m = lda)

    // Given uplo'ed matrix A, return a full symmetrized matrix
    void symmetrize(const UploType uplo,
                    const unsigned int m,
                    double *A);

    void printDenseMatrix(const char *name,
                          const unsigned int m,
                          const unsigned int n,
                          const double *A);

    // Copy matrix with transposition
    void copyDenseMatrix(const unsigned int m,
                         const unsigned int n,
                         double *A,
                         const bool is_transpose,
                         double *B);

    // In-place matrix transposition
    void transposeDenseMatrix(const unsigned int m,
                              const unsigned int n,
                              double *A);

    // Vector utils
    double dotVectors(const unsigned int dim,
                      const double *x, const double *A, const double *y);

    double dotVectors(const unsigned int dim,
                      const double *x,
                      const double *y);

    double dotVectors(const unsigned int dim,
                      const double *x, const int incx,
                      const double *y, const int incy);

    // Matrix utils
    double computeDetDenseMatrix(const unsigned int dim, const double *A);
    int scaleDenseMatrix(const double alpha,
                         const UploType uplo,
                         const unsigned int m,
                         const unsigned int n,
                         double *A);

    int updateSymRankDenseMatrix(const double alpha,
                                 const unsigned int m,
                                 double *x,
                                 const UploType uplo,
                                 const unsigned int lda,
                                 double *A);

    int updateSymRankDenseMatrix(const double alpha,
                                 const unsigned int s,
                                 int *idx, double *x,
                                 const UploType uplo,
                                 const unsigned int m,
                                 double *A);

    int invertDenseMatrix(const unsigned int dim,
                          double *A,
                          double *B,
                          const unsigned int lw,
                          double *W);

    int invertDenseMatrixSVD(const unsigned int dim,
                             double *A,
                             double *B,
                             const unsigned int lw,
                             double *W);

    double checkInversion(const unsigned int dim,
                          double *A,
                          double *Ainv);

    int solveDenseMatrix(const UploType uplo,
                         const unsigned int mA,
                         const unsigned int nA,
                         double *A,
                         const unsigned int mB,
                         const unsigned int nB,
                         double *B,
                         const unsigned int lw,
                         double *W);

    int solveDenseMatrixSVD(const UploType uplo,
                            const unsigned int mA,
                            const unsigned int nA,
                            double *A,
                            const unsigned int mB,
                            const unsigned int nB,
                            double *B,
                            const unsigned int lw,
                            double *W);

    int computeRotationMatrix(const unsigned int dim,
                              const double * n,
                              double *Q);

  };
}

#endif
