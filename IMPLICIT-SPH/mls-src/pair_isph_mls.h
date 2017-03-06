/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(isph/mls,PairISPH_MLS)

#else

#ifndef LMP_PAIR_ISPH_MLS_H
#define LMP_PAIR_ISPH_MLS_H

#include "pair_isph.h"

// kernel
#include "kernel_mls.h"

// basis 
#include "scaled_taylor_monomial.h"

// setup
#include "functor_graph.h"
#include "functor_mls_mass_matrix.h"
#include "functor_mls_mass_matrix_compact_poisson.h"
#include "functor_mls_normal.h"

// math operators
#include "functor_mls_helper.h"
#include "functor_mls_helper_compact_poisson.h"

#include "functor_mls_gradient.h"
#include "functor_mls_gradient_operator.h"
#include "functor_mls_gradient_compact_poisson.h"

#include "functor_mls_divergence.h"

#include "functor_mls_laplacian.h"
#include "functor_mls_laplacian_compact_poisson.h"
#include "functor_mls_laplacian_matrix.h"
#include "functor_mls_laplacian_matrix_compact_poisson.h"

#include "functor_mls_curl.h"
#include "functor_mls_curlcurl.h"

// boundary techniques
#include "functor_mls_boundary_neumann.h"

// gradient based operators
#include "functor_gradient_dot_operator_matrix.h"

#include "functor_curl.h"
#include "functor_curlcurl.h"

// application: pbe
#include "functor_poisson_boltzmann_f.h"
#include "functor_poisson_boltzmann_extra_f.h"
#include "functor_poisson_boltzmann_jacobian.h"

// application: ns - pressure correction
#include "functor_incomp_navier_stokes_helmholtz.h"
#include "functor_incomp_navier_stokes_block_helmholtz.h"
#include "functor_incomp_navier_stokes_poisson.h"
#include "functor_correct_velocity.h"
#include "functor_correct_pressure.h"

#include "functor_advance_time_begin.h"
#include "functor_advance_time_end.h"
#include "functor_compute_shift.h"
#include "functor_apply_shift.h"

// application: ns - velocity correction
#include "functor_ale_advection.h"
#include "functor_ale_advection_matrix.h"
#include "functor_ale_predict_velocity.h"
#include "functor_ale_predict_velocity_curlcurl.h"
#include "functor_ale_correct_velocity.h"
#include "functor_ale_incomp_navier_stokes_poisson_boundary.h"
#include "functor_ale_incomp_navier_stokes_compact_poisson_boundary.h"
#include "functor_ale_incomp_navier_stokes_helmholtz.h"
#include "functor_ale_incomp_navier_stokes_helmholtz_curlcurl.h"

#include "functor_ale_track_particles.h"
#include "functor_ale_apply_shift.h"

// application: ns - post process
#include "functor_traction_vector.h"

// test suite
#include "test_mls.h"

namespace LAMMPS_NS {

  class PairISPH_MLS : public PairISPH {
  public:
    // ** short name for this class 
    typedef class PairISPH_MLS Self; 

    // ** test suite
    friend class MLS::TestSuite;

    // ** common
    friend class FunctorOuter<class PairISPH_MLS>; 

    // ** algorithm specific; have a separate namespace "MLS"
    friend class MLS::FunctorOuterMassMatrix<class PairISPH_MLS>;
    friend class MLS::FunctorOuterMassMatrixCompactPoisson<class PairISPH_MLS>;
    friend class MLS::FunctorOuterNormal<class PairISPH_MLS>;

    // operators;
    friend class MLS::FunctorOuterHelper<class PairISPH_MLS>; 
    friend class MLS::FunctorOuterHelperCompactPoisson<class PairISPH_MLS>; 

    friend class MLS::FunctorOuterGradient<class PairISPH_MLS>; 
    friend class MLS::FunctorOuterGradientOperator<class PairISPH_MLS>;
    friend class MLS::FunctorOuterGradientCompactPoisson<class PairISPH_MLS>; 
    friend class MLS::FunctorOuterDivergence<class PairISPH_MLS>;
    friend class MLS::FunctorOuterLaplacian<class PairISPH_MLS>;
    friend class MLS::FunctorOuterLaplacianCompactPoisson<class PairISPH_MLS>; 
    friend class MLS::FunctorOuterLaplacianMatrix<class PairISPH_MLS>;
    friend class MLS::FunctorOuterLaplacianMatrixCompactPoisson<class PairISPH_MLS>;
    friend class MLS::FunctorOuterCurl<class PairISPH_MLS>;
    friend class MLS::FunctorOuterCurlCurl<class PairISPH_MLS>;

    friend class MLS::FunctorOuterBoundaryNeumann<class PairISPH_MLS>;    

    // gradient base operators
    friend class FunctorOuterGradientDotOperatorMatrix<class PairISPH_MLS,MLS::FunctorOuterGradientOperator>;

    friend class FunctorOuterCurl<class PairISPH_MLS,MLS::FunctorOuterGradient>;
    friend class FunctorOuterCurlCurl<class PairISPH_MLS,MLS::FunctorOuterGradient>;

    // ** application specific; 
    // must have operator parameters as its their template arguments

    // - PB -
    template<class ArgSelf, 
             template<typename> class Laplacian>
    friend class FunctorOuterPoissonBoltzmannF;

    template<class ArgSelf>
    friend class FunctorOuterPoissonBoltzmannExtraF;

    template<class ArgSelf, 
             template<typename> class LaplacianMatrix>
    friend class FunctorOuterPoissonBoltzmannJacobian;

    // - NS - Pressure correction
    template<class ArgSelf, 
             template<typename> class LaplacianMatrix, 
             template<typename> class Gradient>
    friend class FunctorOuterIncompNavierStokesHelmholtz;

    template<class ArgSelf, 
             template<typename> class LaplacianMatrix, 
             template<typename> class Divergence,
             template<typename> class GradientOperator>
    friend class FunctorOuterIncompNavierStokesPoisson;

    template<class ArgSelf, 
             template<typename> class Gradient>
    friend class FunctorOuterCorrectVelocity;

    template<class ArgSelf>
    friend class FunctorOuterCorrectPressure;

    template<class ArgSelf, 
             template<typename> class Gradient>
    friend class FunctorOuterAdvanceTimeBegin;

    template<class ArgSelf>
    friend class FunctorOuterAdvanceTimeEnd;

    template<class ArgSelf>
    friend class FunctorOuterComputeShift;

    template<class ArgSelf, template<typename> class Gradient>
    friend class FunctorOuterApplyShift;

    // - NS - Velocity correction
    template<class ArgSelf,
             template<typename> class GradientOperator>
    friend class FunctorOuterAleAdvection;

    template<class ArgSelf,
             template<typename> class GradientOperator>
    friend class FunctorOuterAleAdvectionMatrix;

    template<class ArgSelf,
             template<typename> class GradientOperator,
             template<typename> class Gradient,
             template<typename> class GradientInner>
    friend class FunctorOuterAlePredictVelocity;

    template<class ArgSelf,
             template<typename> class GradientOperator,
             template<typename> class CurlCurl>
    friend class FunctorOuterAlePredictVelocityCurlCurl;

    template<class ArgSelf,
             template<typename> class Gradient>
    friend class FunctorOuterAleCorrectVelocity;

    template<class ArgSelf, 
             template<typename> class LaplacianMatrix, 
             template<typename> class Divergence,
             template<typename> class GradientOperator,
             template<typename> class NeumannBoundary>
    friend class FunctorOuterAleIncompNavierStokesPoissonBoundary;

    template<class ArgSelf, 
             template<typename> class LaplacianMatrixCompactPoisson>
    friend class FunctorOuterAleIncompNavierStokesCompactPoissonBoundary;

    template<class ArgSelf,
             template<typename> class LaplacianMatrix,
             template<typename> class GradientOperator,
             template<typename> class Gradient,
             template<typename> class GradientInner>
    friend class FunctorOuterAleIncompNavierStokesHelmholtz;

    template<class ArgSelf,
             template<typename> class LaplacianMatrix,
             template<typename> class GradientOperator,
             template<typename> class CurlCurl>
    friend class FunctorOuterAleIncompNavierStokesHelmholtzCurlCurl;

    template<class ArgSelf>
    friend class FunctorOuterAleTrackParticles;

    template<class ArgSelf>
    friend class FunctorOuterAleApplyShift;

    // NS - Post process
    template<class ArgSelf,
             template<typename> class Gradient>
    friend class FunctorOuterTractionVector;

  public:

    PairISPH_MLS(class LAMMPS *);
    virtual~PairISPH_MLS();

    // interface is declared in base pair
    // void compute(int, int); 
    void coeff(int, char **);

    // overload comm routines to specialize them for this pairing
    //int pack_forward_comm(int, int *, double *, int, int *);
    //void unpack_forward_comm(int, int, double *);

    void grow(double extend = 1.0);

    void computePre(const int mask = 0);
    void computeMassMatrix();
    void computeMassMatrixCompactPoisson();
    void computeNormals();

    void enterCompactPoisson();
    void exitCompactPoisson();
    bool isCompactPoisson();

    void computeUnitTests();

    // PB NOX
    bool computeF(const Epetra_Vector &x, Epetra_Vector &f,
                  const FillType flag = Residual);
    bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &J);
    bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M,
                               Teuchos::ParameterList *param = NULL);

    void computePsiGradient();
    
    // NS - Projection 
    void computeHelmholtz(double *b);
    void computePoisson(double *b);
    void correctVelocity();
    void correctPressure();

    // NS - ALE
    void predictAleVelocity();
    void computeAlePoisson(double *b);
    void correctAleVelocity();
    void computeAleHelmholtz(double *b);

    // NS - debugging
    void computeAleDebug(double **v, double *div, double **curl);

    // NS - Time
    void advanceTime();
    void shiftParticles(const double shift,
                        const double shiftcut,
                        const double boundaryweight);

    // NS - Post processing
    void computeVelocityCurl(double **v, double **curl);
    void computeVelocityDivergence(double **v, double *div);
    void computeTractionVector(double *p, double **v, double **traction);

  protected:

    // ----------------
    // user attributes
    // ----------------

    // mass matrix: a - without penalty, b - with penalty for poisson solver
    // use M for interface
    double **M, **Ma, **Mb;

    // MLS; approx order
    int np;

    // kernel function (created when pair is created)
    KernelFuncMLS *kernel;

    // basis function
    ScaledTaylorMonomial *basis;

    // interpolation flag
    bool is_interpolation_enabled;

    // compact mls flag for poisson solver
    struct {
      bool is_enabled;
      double tau_interior, tau_boundary;

      double interior_penalty_coeff(const double rth) const { return tau_interior*pow(rth, 4); }
      double boundary_penalty_coeff(const double rth) const { return tau_boundary*pow(rth, 2); }
    } cp;

  };

}

#endif
#endif
