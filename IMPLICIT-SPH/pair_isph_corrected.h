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

PairStyle(isph/corrected,PairISPH_Corrected)

#else

#ifndef LMP_PAIR_ISPH_CORRECTED_H
#define LMP_PAIR_ISPH_CORRECTED_H

#include "pair_isph.h"

// kernel
#include "kernel.h"

// setup
#include "functor_graph.h"
#include "functor_volume.h"
#include "functor_gradient_correction.h"
#include "functor_laplacian_correction.h"
#include "functor_normal.h"

// math operators
#include "functor_smooth_field.h"
#include "functor_gradient.h"
#include "functor_gradient_operator.h"

#include "functor_divergence.h"

#include "functor_laplacian.h"
#include "functor_laplacian_matrix.h"

// boundary techniques
//#include "functor_boundary_morris_normal.h"
#include "functor_boundary_morris_holmes.h"
#include "functor_boundary_navier_slip.h"
#include "functor_boundary_dirichlet.h"
#include "functor_boundary_neumann.h"

// gradient based operators
//#include "functor_curl.h"
//#include "functor_curlcurl.h"

#include "functor_gradient_dot_operator_matrix.h"

// application: pbe
#include "functor_poisson_boltzmann_f.h"
#include "functor_poisson_boltzmann_extra_f.h"
#include "functor_poisson_boltzmann_jacobian.h"

// application: surface tension
#include "functor_phase_gradient.h"
#include "functor_phase_divergence.h"
#include "functor_phase_divergence_adami.h"
#include "functor_correct_phase_normal.h"
#include "functor_continuum_surface_force.h"

#include "functor_pairwise_force.h"

// application: random stress
#include "functor_random_stress.h"

// application: solute transport
#include "functor_solute_transport.h"

// application: non-uniform applied electric potential
#include "functor_applied_electric_potential.h"

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
// #include "functor_ale_advection.h"
// #include "functor_ale_advection_matrix.h"
// #include "functor_ale_predict_velocity.h"
// #include "functor_ale_correct_velocity.h"
// #include "functor_ale_incomp_navier_stokes_poisson.h"
// #include "functor_ale_incomp_navier_stokes_helmholtz.h"

// #include "functor_ale_track_particles.h"
// #include "functor_ale_apply_shift.h"

namespace LAMMPS_NS {
  
  class PairISPH_Corrected : public PairISPH {
  private:
    // ** short name for this class
    typedef class PairISPH_Corrected Self;

    // ** common; 
    friend class FunctorOuter<class PairISPH_Corrected>;

    // ** algorithm specific; have a separate namespace "Corrected"
    friend class Corrected::FunctorOuterVolume<class PairISPH_Corrected>;
    friend class Corrected::FunctorOuterNormal<class PairISPH_Corrected>;
    friend class Corrected::FunctorOuterGradientCorrection<class PairISPH_Corrected>;
    friend class Corrected::FunctorOuterLaplacianCorrection<class PairISPH_Corrected>;

    friend class Corrected::FunctorOuterPhaseGradient<class PairISPH_Corrected>;
    friend class Corrected::FunctorOuterPhaseDivergence<class PairISPH_Corrected>;
    friend class Corrected::FunctorOuterPhaseDivergenceAdami<class PairISPH_Corrected>;
    friend class Corrected::FunctorOuterCorrectPhaseNormal<class PairISPH_Corrected>;
    
    // operators; 
    // ----------------------------------------------------------------
    friend class Corrected::FunctorOuterSmoothField<class PairISPH_Corrected>;

    template<class ArgSelf, bool AntiSymmetric>
    friend class Corrected::FunctorOuterGradient;

    template<class ArgSelf, bool AntiSymmetric>
    friend class Corrected::FunctorOuterGradient_MorrisHolmes;

    friend class Corrected::FunctorOuterGradientOperator<class PairISPH_Corrected>;
    friend class Corrected::FunctorOuterGradientOperator_MorrisHolmes<class PairISPH_Corrected>;

    template<class ArgSelf, bool AntiSymmetric>
    friend class Corrected::FunctorOuterDivergence;

    template<class ArgSelf, bool AntiSymmetric>
    friend class Corrected::FunctorOuterDivergence_MorrisHolmes;

    friend class Corrected::FunctorOuterLaplacian<class PairISPH_Corrected>;
    friend class Corrected::FunctorOuterLaplacian_MorrisHolmes<class PairISPH_Corrected>;
    //friend class Corrected::FunctorOuterLaplacian_MorrisNormal<class PairISPH_Corrected>;

    friend class Corrected::FunctorOuterLaplacianHelper<class PairISPH_Corrected>;

    template<class ArgSelf, bool AntiSymmetric>
    friend class Corrected::FunctorOuterLaplacianMatrix;

    template<class ArgSelf, bool AntiSymmetric>
    friend class Corrected::FunctorOuterLaplacianMatrix_MorrisHolmes;

    //template<class ArgSelf, bool AntiSymmetric>
    //friend class Corrected::FunctorOuterLaplacianMatrix_MorrisNormal;

    // TODO:: change the boundary condition technique name to include where they 
    //        are applied; these boundary technique is not general to cope with 
    //        all possible situations
    friend class Corrected::FunctorOuterBoundaryNavierSlip<class PairISPH_Corrected>;
    friend class Corrected::FunctorOuterBoundaryDirichlet<class PairISPH_Corrected>;
    friend class Corrected::FunctorOuterBoundaryNeumann<class PairISPH_Corrected>;

    // gradient based operator
    // ----------------------------------------------------------------
    template<class ArgSelf, 
             template<typename> class GradientOperator>
    friend class FunctorOuterGradientDotOperatorMatrix;

    // template<class ArgSelf, 
    //          template<typename> class Gradient>
    // friend class FunctorOuterCurl;

    // template<class ArgSelf, 
    //          template<typename> class Gradient,
    //          template<typename> class GradientInner>
    // friend class FunctorOuterCurlCurl;

    // ** application specific; 
    // must have operator parameters as its their template arguments

    // - PB - Poisson Boltzmann
    template<class ArgSelf, 
             template<typename> class Laplacian> 
    friend class FunctorOuterPoissonBoltzmannF;

    template<class ArgSelf>
    friend class FunctorOuterPoissonBoltzmannExtraF; 

    template<class ArgSelf, 
             template<typename> class LaplacianMatrix>
    friend class FunctorOuterPoissonBoltzmannJacobian;

    // - ST - Surface Tension
    template<class ArgSelf, 
             template<typename> class Divergence>
    friend class FunctorOuterContinuumSurfaceForce;

    template<class ArgSelf>
    friend class FunctorOuterPairwiseForce;

    // RS - Random stress
    template<class ArgSelf, 
             template<typename> class Divergence>
    friend class FunctorOuterRandomStress;

    // - RT - Solute Transport
    template<class ArgSelf, 
             template<typename> class LaplacianMatrix,
             template<typename> class GradientOperator>
    friend class FunctorOuterSoluteTransport;

    // - EAPP - Non-uniform applied electric field
    template<class ArgSelf, 
             template<typename> class LaplacianMatrix>
    friend class FunctorOuterAppliedElectricPotential;

    // - NS - Pressure correction 
    template<class ArgSelf, 
             template<typename> class LaplacianMatrix, 
             template<typename> class Gradient>
    friend class FunctorOuterIncompNavierStokesHelmholtz;

    template<class ArgSelf,
             template<typename> class LaplacianMatrix,
             template<typename> class Gradient,
             template<typename> class NaturalBoundary>
    friend class FunctorOuterIncompNavierStokesBlockHelmholtz;

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

    template<class ArgSelf, 
             template<typename> class Gradient>
    friend class FunctorOuterApplyShift;

    // - NS - Velocity correction 
    // template<class ArgSelf, 
    //          template<typename> class GradientOperator>
    // friend class FunctorOuterAleAdvection;

    // template<class ArgSelf, 
    //          template<typename> class GradientOperator>
    // friend class FunctorOuterAleAdvectionMatrix;

    // template<class ArgSelf, 
    //          template<typename> class GradientOperator,
    //          template<typename> class Gradient, 
    //          template<typename> class GradientInner>
    // friend class FunctorOuterAlePredictVelocity;

    // template<class ArgSelf, 
    //          template<typename> class GradientOperator,
    //          template<typename> class Laplacian> 
    // friend class FunctorOuterAlePredictVelocityLaplacian;

    // template<class ArgSelf, 
    //          template<typename> class Gradient>
    // friend class FunctorOuterAleCorrectVelocity;

    // template<class ArgSelf, 
    //          template<typename> class LaplacianMatrix, 
    //          template<typename> class Divergence,
    //          template<typename> class NeumannBoundary>
    // friend class FunctorOuterAleIncompNavierStokesPoisson;

    // template<class ArgSelf, 
    //          template<typename> class LaplacianMatrix,
    //          template<typename> class GradientOperator,
    //          template<typename> class Gradient,
    //          template<typename> class GradientInner>
    // friend class FunctorOuterAleIncompNavierStokesHelmholtz;

    // template<class ArgSelf, 
    //          template<typename> class LaplacianMatrix,
    //          template<typename> class GradientOperator,
    //          template<typename> class Laplacian>
    // friend class FunctorOuterAleIncompNavierStokesHelmholtzLaplacian;

    // template<class ArgSelf>
    // friend class FunctorOuterAleTrackParticles;

    // template<class ArgSelf>
    // friend class FunctorOuterAleApplyShift;

  public:

    PairISPH_Corrected(class LAMMPS *);
    virtual~PairISPH_Corrected();

    // interface is declared in base pair
    //void compute(int, int);
    void coeff(int, char **);

    // overload comm routines to specialize them for this pairing
    int pack_forward_comm(int, int *, double *, int, int *);
    void unpack_forward_comm(int, int, double *);

    void grow(double extend = 1.0);

    void computePre(const int mask = 0);
    void computeVolumes();
    void computeGradientCorrection();
    void computeLaplacianCorrection();
    void computeNormals();

    // PB NOX 
    bool computeF(const Epetra_Vector &x, Epetra_Vector &f,
                  const FillType flag = Residual);
    bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &J);
    bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M,
                               Teuchos::ParameterList *param = NULL);
    void computePsiGradient();
    void computeAppliedElectricPotential(double *b);
    void computePhiGradient();

    // ST - Surface tension
    void computeSurfaceTension();
    void computeSurfaceTension_ContinuumSurfaceForce();
    void computeSurfaceTension_PairwiseForce();

    // RS - Random stress
    void computeForceFromRandomStress();  

    // RT - Solute transport
    void computeSoluteTransportSpecies(const double dcoeff, double *b);

    // NS - Projection 
    void computeHelmholtz(double *b);
    void computeBlockHelmholtz(double *b);
    void computePoisson(double *b);
    void correctVelocity();
    void correctPressure();

    // NS - ALE
    void predictAleVelocity();
    void computeAlePoisson(double *b);
    void correctAleVelocity();
    void computeAleHelmholtz(double *b);

    // NS - Time 
    void advanceTime();
    void shiftParticles(const double shift,
                        const double shiftcut,
                        const double solidweight);

    // NS - Post process
    void computeVelocityDivergence(double **v, double *div);

  protected:

    // ----------------
    // user attributes
    // ----------------

    // correction tensors for particles
    double **Gc, **Lc, *Gi, *Li;   

    // normal to fluid/structure interface and particle number density
    double *pnd, *bd_coord;

    // kernel function (created when pair is created)
    KernelFunction *kernel;

    // when morris mirroring is applied to operators, it chooses 
    // max(d_a, morris_safe_coeff*h) to prevent incident blow-up of mirror variables
    double morris_safe_coeff; 
  };

}

#endif
#endif
