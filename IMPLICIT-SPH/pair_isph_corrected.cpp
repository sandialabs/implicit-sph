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

#include "math.h"
#include "stdlib.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h" 
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "random_mars.h"

// user
#include <iostream>
#include <fstream>
#include <string>
#include "utils.h"

#include "kernel.h"
#include "kernel_cubic.h"
#include "kernel_quintic.h"
#include "kernel_wendland.h" 

#include "pair_isph.h"
#include "pair_isph_corrected.h"

#include "pair_for.h"

// symmetric is consistent form
// anti symmetirc preserve momentum
#define FunctorSymType(SymType)
#define FunctorSelect                           \
  if (ns.is_momentum_preserve_operator_used) {  \
    FunctorSymType(AntiSymmetric);              \
  } else {                                      \
    FunctorSymType(Symmetric);                  \
  }                                         

using namespace std;
using namespace LAMMPS_NS;

typedef class PairISPH_Corrected Self;

// common operators
// ================
typedef class FunctorOuterNormalizeVector<Self> FunctorNormalizeVector;
typedef class FunctorOuterGraph<Self>           FunctorGraph;

// math-related operators for "corrected" scheme
// =============================================
typedef class Corrected::FunctorOuterVolume<Self>              FunctorVolume;
typedef class Corrected::FunctorOuterNormal<Self>              FunctorNormal;
typedef class Corrected::FunctorOuterGradientCorrection<Self>  FunctorGradientCorrection;
typedef class Corrected::FunctorOuterLaplacianCorrection<Self> FunctorLaplacianCorrection;

typedef class Corrected::FunctorOuterPhaseGradient<Self>        FunctorPhaseGradient;
//typedef class Corrected::FunctorOuterPhaseDivergence<Self>      FunctorPhaseDivergence;
//typedef class Corrected::FunctorOuterPhaseDivergenceAdami<Self> FunctorPhaseDivergence;
typedef class Corrected::FunctorOuterCorrectPhaseNormal<Self>   FunctorCorrectPhaseNormal;

typedef class Corrected::FunctorOuterSmoothField<Self>     FunctorSmoothField;
typedef class Corrected::FunctorOuterGradient<Self,false>        FunctorGradient;
typedef class Corrected::FunctorOuterDivergence<Self,false>      FunctorDivergence;
typedef class Corrected::FunctorOuterLaplacian<Self>       FunctorLaplacian;
typedef class Corrected::FunctorOuterLaplacianMatrix<Self,false> FunctorLaplacianMatrix;

// sph-related boundary operators for "corrected" scheme
// =====================================================
typedef class Corrected::FunctorOuterDivergence_MorrisHolmes<Self> FunctorDivergence_MorrisHolmes;

typedef class Corrected::FunctorOuterLaplacian_MorrisHolmes<Self>       FunctorLaplacian_MorrisHolmes;
typedef class Corrected::FunctorOuterLaplacianMatrix_MorrisHolmes<Self> FunctorLaplacianMatrix_MorrisHolmes;

//typedef class Corrected::FunctorOuterLaplacian_MorrisNormal<Self>       FunctorLaplacian_MorrisNormal;
//typedef class Corrected::FunctorOuterLaplacianMatrix_MorrisNormal<Self> FunctorLaplacianMatrix_MorrisNormal;

typedef class Corrected::FunctorOuterBoundaryNavierSlip<Self> FunctorBoundaryNavierSlip;
typedef class Corrected::FunctorOuterBoundaryDirichlet<Self>  FunctorBoundaryDirichlet;
typedef class Corrected::FunctorOuterBoundaryNeumann<Self>    FunctorBoundaryNeumann;  

// Poisson-Boltzmann NOX-related operators for "corrected" scheme
// ==============================================================
typedef class FunctorOuterPoissonBoltzmannF<Self,
                                            Corrected::FunctorOuterLaplacian
                                            > FunctorPoissonBoltzmannF;
typedef class FunctorOuterPoissonBoltzmannF<Self,
                                            Corrected::FunctorOuterLaplacian_MorrisHolmes
                                            > FunctorPoissonBoltzmannF_MorrisHolmes;
// typedef class FunctorOuterPoissonBoltzmannF<Self,
//                                             Corrected::FunctorOuterLaplacian_MorrisNormal
//                                             > FunctorPoissonBoltzmannF_MorrisNormal;

typedef class FunctorOuterPoissonBoltzmannExtraF<Self> FunctorPoissonBoltzmannExtraF;

typedef class FunctorOuterPoissonBoltzmannJacobian<Self,
                                                   Corrected::FunctorOuterLaplacianMatrixSymmetric
                                                   > FunctorPoissonBoltzmannJacobian;
typedef class FunctorOuterPoissonBoltzmannJacobian<Self,
                                                   Corrected::FunctorOuterLaplacianMatrixSymmetric_MorrisHolmes
                                                   > FunctorPoissonBoltzmannJacobian_MorrisHolmes;
// typedef class FunctorOuterPoissonBoltzmannJacobian<Self,
//                                                    Corrected::FunctorOuterLaplacianMatrixSymmetric_MorrisNormal
//                                                    > FunctorPoissonBoltzmannJacobian_MorrisNormal;


// Surface tension scheme
// ======================
typedef class FunctorOuterContinuumSurfaceForce<Self,
                                                Corrected::FunctorOuterPhaseDivergenceAdami
                                                > FunctorContinuumSurfaceForce;
typedef class FunctorOuterPairwiseForce<Self> FunctorPairwiseForce;

// Random stress tensor
// ====================
typedef class FunctorOuterRandomStress<Self,
                                       Corrected::FunctorOuterDivergenceAntiSymmetric
                                       > FunctorRandomStress;

// Solute transport
// ==================
typedef class FunctorOuterSoluteTransport<Self,
                                            Corrected::FunctorOuterLaplacianMatrixSymmetric,
                                            Corrected::FunctorOuterGradientOperator
                                            > FunctorSoluteTransport;
// Non-uniform appled electric potential
// =====================================
typedef class FunctorOuterAppliedElectricPotential<Self,
                                                   Corrected::FunctorOuterLaplacianMatrixSymmetric
                                                   > FunctorAppliedElectricPotential;

// Navier-Stokes pressure correction scheme
// ========================================

typedef class FunctorOuterAdvanceTimeEnd<Self> FunctorAdvanceTimeEnd;
typedef class FunctorOuterComputeShift<Self> FunctorComputeShift;


// ** consistency family
typedef class FunctorOuterIncompNavierStokesHelmholtz<Self,
                                                      Corrected::FunctorOuterLaplacianMatrixAntiSymmetric,
                                                      Corrected::FunctorOuterGradientAntiSymmetric
                                                      > FunctorIncompNavierStokesHelmholtz_AntiSymmetric;
typedef class FunctorOuterIncompNavierStokesHelmholtz<Self,
                                                      Corrected::FunctorOuterLaplacianMatrixAntiSymmetric_MorrisHolmes,
                                                      Corrected::FunctorOuterGradientAntiSymmetric
                                                      > FunctorIncompNavierStokesHelmholtz_MorrisHolmes_AntiSymmetric;

typedef class FunctorOuterIncompNavierStokesBlockHelmholtz<Self,
                                                      Corrected::FunctorOuterLaplacianMatrixAntiSymmetric_MorrisHolmes,
                                                      Corrected::FunctorOuterGradientAntiSymmetric,
                                                      Corrected::FunctorOuterBoundaryNavierSlip
                                                      > FunctorIncompNavierStokesBlockHelmholtz_MorrisHolmes_AntiSymmetric;

typedef class FunctorOuterIncompNavierStokesPoisson<Self,
                                                    Corrected::FunctorOuterLaplacianMatrixAntiSymmetric,
                                                    Corrected::FunctorOuterDivergenceAntiSymmetric,
                                                    Corrected::FunctorOuterGradientOperator
                                                    > FunctorIncompNavierStokesPoisson_AntiSymmetric;
typedef class FunctorOuterIncompNavierStokesPoisson<Self,
                                                    Corrected::FunctorOuterLaplacianMatrixAntiSymmetric,
                                                    Corrected::FunctorOuterDivergenceAntiSymmetric_MorrisHolmes,
                                                    Corrected::FunctorOuterGradientOperator
                                                    > FunctorIncompNavierStokesPoisson_MorrisHolmes_AntiSymmetric;

typedef class FunctorOuterCorrectVelocity<Self,
                                          Corrected::FunctorOuterGradientAntiSymmetric
                                          > FunctorCorrectVelocity_AntiSymmetric;
typedef class FunctorOuterCorrectPressure<Self> FunctorCorrectPressure_AntiSymmetric;

typedef class FunctorOuterAdvanceTimeBegin<Self,
                                           Corrected::FunctorOuterGradientAntiSymmetric
                                           > FunctorAdvanceTimeBegin_AntiSymmetric;

typedef class FunctorOuterApplyShift<Self,
                                     Corrected::FunctorOuterGradientAntiSymmetric
                                     > FunctorApplyShift_AntiSymmetric;


// ** momentum preserve family
typedef class FunctorOuterIncompNavierStokesHelmholtz<Self,
                                                      Corrected::FunctorOuterLaplacianMatrixSymmetric,
                                                      Corrected::FunctorOuterGradientSymmetric
                                                      > FunctorIncompNavierStokesHelmholtz_Symmetric;
typedef class FunctorOuterIncompNavierStokesHelmholtz<Self,
                                                      Corrected::FunctorOuterLaplacianMatrixSymmetric_MorrisHolmes,
                                                      Corrected::FunctorOuterGradientSymmetric
                                                      > FunctorIncompNavierStokesHelmholtz_MorrisHolmes_Symmetric;

typedef class FunctorOuterIncompNavierStokesBlockHelmholtz<Self,
                                                      Corrected::FunctorOuterLaplacianMatrixSymmetric_MorrisHolmes,
                                                      Corrected::FunctorOuterGradientSymmetric,
                                                      Corrected::FunctorOuterBoundaryNavierSlip
                                                      > FunctorIncompNavierStokesBlockHelmholtz_MorrisHolmes_Symmetric;

typedef class FunctorOuterIncompNavierStokesPoisson<Self,
                                                    Corrected::FunctorOuterLaplacianMatrixSymmetric,
                                                    Corrected::FunctorOuterDivergenceSymmetric,
                                                    Corrected::FunctorOuterGradientOperator
                                                    > FunctorIncompNavierStokesPoisson_Symmetric;
typedef class FunctorOuterIncompNavierStokesPoisson<Self,
                                                    Corrected::FunctorOuterLaplacianMatrixSymmetric,
                                                    Corrected::FunctorOuterDivergenceSymmetric_MorrisHolmes,
                                                    Corrected::FunctorOuterGradientOperator
                                                    > FunctorIncompNavierStokesPoisson_MorrisHolmes_Symmetric;

typedef class FunctorOuterCorrectVelocity<Self,
                                          Corrected::FunctorOuterGradientSymmetric
                                          > FunctorCorrectVelocity_Symmetric;
typedef class FunctorOuterCorrectPressure<Self> FunctorCorrectPressure_Symmetric;

typedef class FunctorOuterAdvanceTimeBegin<Self,
                                           Corrected::FunctorOuterGradientSymmetric
                                           > FunctorAdvanceTimeBegin_Symmetric;

typedef class FunctorOuterApplyShift<Self,
                                     Corrected::FunctorOuterGradientSymmetric
                                     > FunctorApplyShift_Symmetric;

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::PairISPH_Corrected"
PairISPH_Corrected::PairISPH_Corrected(LAMMPS *lmp) 
  : PairISPH(lmp), 
    Gc(NULL), Lc(NULL), Gi(NULL), Li(NULL),
    pnd(NULL), bd_coord(NULL) {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor initalizationTimerMonitor(*g_timer_init);

  discretization = Corrected;

  nmax = 0;

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::~PairISPH_Corrected"
PairISPH_Corrected::~PairISPH_Corrected() {
  FUNCT_ENTER(comm->me);

  delete kernel;

  memory->destroy(Gc);
  memory->destroy(Lc);

  memory->destroy(Gi);
  memory->destroy(Li);

  memory->destroy(pnd);
  memory->destroy(bd_coord);

  FUNCT_EXIT(comm->me);
}

// BEGIN::COMMON
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::grow"
void PairISPH_Corrected::grow(double extend) {
  FUNCT_ENTER(comm->me);

  if (atom->nmax > nmax) {
    // nmax is updated here
    PairISPH::grow(extend);

    int dim  = domain->dimension;
    int dimL = dim*(dim+1)/2;

    memory->grow(Gc,       nmax, dim*dim, "pair:Gc");
    memory->grow(Lc,       nmax,    dimL, "pair:Lc");

    memory->grow(Gi,       dim*dim,       "pair:Gi");
    memory->grow(Li,       dimL,          "pair:Li");

    memory->grow(pnd,      nmax,          "pair:pnd");
    memory->grow(bd_coord, nmax,          "pair:bd_coord");

  }
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computePre"
void PairISPH_Corrected::computePre(const int mask) {
  FUNCT_ENTER(comm->me);

  computeVolumes();
  computeGradientCorrection();
  computeLaplacianCorrection();

  if (boundary_particle != NoParticle)
    computeNormals();

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeVolumes"
void PairISPH_Corrected::computeVolumes() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_volumes);

  // compute for all interacting with all neighbors
  FunctorVolume functor(this);
  PairFor(functor, functor.getNumberOfWork());

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeGradientCorrection"
void PairISPH_Corrected::computeGradientCorrection() {
  FUNCT_ENTER(comm->me);
  
  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_gradient_correction);

  // compute for all interacting with all neighbors
  FunctorGradientCorrection functor(this);
  PairFor(functor, functor.getNumberOfWork());

  // identity matrix
  const int dim = domain->dimension;
  for (int k2=0;k2<dim;++k2)
    for (int k1=0;k1<dim;++k1)
      VIEW2(Gi, dim, k1, k2) = (k1 == k2);

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeLaplacianCorrection"
void PairISPH_Corrected::computeLaplacianCorrection() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_laplacian_correction);

  // compute for all interacting with all neighbors
  FunctorLaplacianCorrection functor(this);
  PairFor(functor, functor.getNumberOfWork());

  const int dim = domain->dimension;
  for (int k2=0, op=0;k2<dim;++k2)
    for (int k1=0;k1<(k2+1);++k1,++op)
      Li[op] = (k1 == k2);

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeNormals"
void PairISPH_Corrected::computeNormals() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_normals);
  
  // this should be more improved 
  // orientaiton is set from solid to fluid particles
  double orient[MaxParticleKind] = {};
  orient[Fluid]           = -1.0;
  orient[BufferDirichlet] = -1.0;
  orient[BufferNeumann]   = -1.0;
  orient[Solid]           =  1.0;
  orient[Boundary]        =  1.0;

  const bool use_bd_coord = ( (pb.boundary == Dirichlet || pb.boundary == MorrisNormal) ||
                              (ns.boundary == Dirichlet || ns.boundary == MorrisNormal) );

  switch (boundary_particle) {
  case Boundary: {
    // do nothing
    break;
  }
  case Solid: {
    if (use_part) {
      FilterExclusiveOr filter;
      FunctorNormal functor(this, atom->part,
                            &orient[0], normal, pnd, (use_bd_coord ? bd_coord : NULL));
      functor.setFilter(&filter);
      PairFor(functor, functor.getNumberOfWork());
    } else {
      {
        FilterBinary filter;
        filter.setPairYes(Fluid, Solid);
        FunctorNormal functor(this, NULL,
                              &orient[0], normal, pnd, (use_bd_coord ? bd_coord : NULL));
        functor.setFilter(&filter);
        PairFor(functor, functor.getNumberOfWork());
      }
      {
        FilterBinary filter;
        filter.setPairYes(Solid, Fluid);
        FunctorNormal functor(this, NULL,
                              &orient[0], normal, pnd, (use_bd_coord ? bd_coord : NULL));
        functor.setFilter(&filter);
        PairFor(functor, functor.getNumberOfWork());
      }
    }
    break;
  }
  }

  comm_variable = NormalVector;
  comm_forward = 4;
  comm->forward_comm_pair(this);

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
// END::COMMON

// BEGIN::PBE NOX
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeF"
bool PairISPH_Corrected::computeF(const Epetra_Vector &psitmp,
                                  Epetra_Vector &f,
                                  const FillType flag) {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_f_poisson_boltzmann);

  // update psi; this is iteratively called from NOX 
  psitmp.ExtractCopy(atom->psi);
  comm_variable = Psi;
  comm_forward = 1;
  comm->forward_comm_pair(this);  

  switch (pb.boundary) {
  case HomogeneousNeumann:  
    // this considers only Fluid and Fluid interactions, 
    // others are for Fluid and All interactions
  case NoBoundaryCond:     
  case ConstExtension: {
    FunctorPoissonBoltzmannF functor(this, atom->psi, atom->psi0, atom->eps, &f);
    PairFor(functor, functor.getNumberOfWork());    
    break;
  }
  case MorrisHolmes: {
    FunctorPoissonBoltzmannF_MorrisHolmes functor(this, atom->psi, atom->psi0, atom->eps, &f);
    PairFor(functor, functor.getNumberOfWork());    
    break;
  }
  case MorrisNormal: {
    //FunctorPoissonBoltzmannF_MorrisNormal functor(this, atom->psi, atom->psi0, atom->eps, &f);
    //PairFor(functor, functor.getNumberOfWork());
    break;
  }
  }

  // this is for verification with a manufactured solution
  if (pb.has_extra_f) {
    FunctorPoissonBoltzmannExtraF functor(this, 
                                          atom->psi, atom->psi0, atom->eps, &f,
                                          pb.extra_f);
    PairFor(functor, functor.getNumberOfWork());
  }

  FUNCT_EXIT(comm->me);
  return true;
} 

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeJacobian"
bool PairISPH_Corrected::computeJacobian(const Epetra_Vector &psitmp,
                                         Epetra_Operator &J) {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_jacobian_poisson_boltzmann);

  // no need to communicate psi; ghosts are not used 

  // psitmp and J are not explicitly used;
  // we already have those references; psi and matrix
  switch (pb.boundary) {
  case HomogeneousNeumann:  
    // this considers only Fluid and Fluid interactions, 
    // others are for Fluid and All interactions
  case NoBoundaryCond:     
  case ConstExtension: {
    FunctorPoissonBoltzmannJacobian functor(this, atom->psi, atom->eps);
    PairFor(functor, functor.getNumberOfWork());    
    break;
  }
  case MorrisHolmes: {
    FunctorPoissonBoltzmannJacobian_MorrisHolmes functor(this, atom->psi, atom->eps);
    PairFor(functor, functor.getNumberOfWork());    
    break;
  }
  case MorrisNormal: {
    //FunctorPoissonBoltzmannJacobian_MorrisNormal functor(this, atom->psi, atom->eps);
    //PairFor(functor, functor.getNumberOfWork());
    break;
  }
  }

  FUNCT_EXIT(comm->me);
  return true; 
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computePreconditioner"
bool PairISPH_Corrected::computePreconditioner(const Epetra_Vector &x,
                                               Epetra_Operator &M,
                                               Teuchos::ParameterList *param) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "We use Jacobian only at this moment: isph");  
  FUNCT_EXIT(comm->me);
  return false;
}  

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computePsiGradient"
void PairISPH_Corrected::computePsiGradient() {
  FUNCT_ENTER(comm->me);

  // cleanup
  memset(&atom->psigrad[0][0], 0, sizeof(double)*3*atom->nlocal);

  FilterBinary filter;
  filter.setPairYes(Fluid, All);

  switch (pb.boundary) {
  case MorrisHolmes: {
    Corrected::FunctorOuterGradient_MorrisHolmes<Self> functor(this, atom->psi, 1.0, atom->psigrad);
    functor.setFilter(&filter);
    PairFor(functor, functor.getNumberOfWork());
    break;
  }
  default: {
    Corrected::FunctorOuterGradient<Self> functor(this, atom->psi, 1.0, atom->psigrad);
    functor.setFilter(&filter);
    PairFor(functor, functor.getNumberOfWork());
    break;
  }
  }

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeAppliedElectricPotential"
void PairISPH_Corrected::computeAppliedElectricPotential(double *b) {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_applied_electric_potential);
  
  double *sigmatmp = NULL, *sigma = NULL;
  if (ae.smooth_phi) {
    const int nmax = atom->nmax, nlocal = atom->nlocal;
    memory->grow(sigmatmp, nmax, "pair:computeAppliedElectricPotential:sigmatmp");
    memcpy(sigmatmp, atom->sigma, sizeof(double)*nmax);

    for (int iter=0;iter<ae.smooth_phi;++iter) {
      // smooth interface between solid and fluid (all, all)
      FunctorSmoothField functor(this, sigmatmp, work);    
      PairFor(functor, functor.getNumberOfWork());    

      memcpy(sigmatmp, work, sizeof(double)*nlocal);

      comm_variable = WorkScalar;
      comm_forward = 1;
      comm->forward_comm_pair(this);
    }
    memory->destroy(sigmatmp);
    sigma = work;
  } else {
    sigma = atom->sigma;
  }

  // consider all solid particles as fluid particles
  // - compute non-homogeneous laplace equation with buffer zone
  {
    if (ae.is_solid2fluid)
      mapParticleKind(Solid, Fluid);

    // homogeneous neumann boundary on solid
    FunctorAppliedElectricPotential functor(this, 
                                            sigma, atom->phi,
                                            normal,
                                            b);

    PairFor(functor, functor.getNumberOfWork());
    
    if (ae.is_solid2fluid)
      mapParticleKind(Solid, Solid);
  }

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computePhiGradient"
void PairISPH_Corrected::computePhiGradient() {
  FUNCT_ENTER(comm->me);

  const int nlocal = atom->nlocal;
  const int *type  = atom->type;

  // cleanup
  memset(&atom->phigrad[0][0], 0, sizeof(double)*3*nlocal);

  {
    FilterBinary filter;
    filter.setPairYes(Fluid, Fluid);
    
    Corrected::FunctorOuterGradient<Self> functor(this, atom->phi, 1.0, atom->phigrad);
    functor.setFilter(&filter);

    PairFor(functor, functor.getNumberOfWork());
  }
  {
    double **phigrad = atom->phigrad;

    // in the buffer area, assign the constant phigrad from eapp 
    for (int i=0;i<nlocal;++i) {
      const int ikind = getParticleKind(type[i]);
      if (ikind == BufferDirichlet || ikind == BufferNeumann) {
        phigrad[i][0] = -ae.e[0];
        phigrad[i][1] = -ae.e[1];
        phigrad[i][2] = -ae.e[2];
      }
    }
  }

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
// END::PBE NOX

// BEGIN::ST
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeSurfaceTension"
void PairISPH_Corrected::computeSurfaceTension() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_surface_tension);  

  switch (st.model) {
  case ContinuumSurfaceForce:
    computeSurfaceTension_ContinuumSurfaceForce();
    break;
  case PairwiseForce:
    computeSurfaceTension_PairwiseForce();
    break;
  default:
    error->all(FLERR, "Surface tension model is not supported: ContinuumSurfaceForce, PairwiseForce");    
  }

  FUNCT_EXIT(comm->me);
}

#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeSurfaceTension_ContinuumSurfaceForce"
void PairISPH_Corrected::computeSurfaceTension_ContinuumSurfaceForce() {
  FUNCT_ENTER(comm->me);

  FilterBinary filter;
  filter.setPairYes(Fluid, Fluid);

  clearCommArray(WorkVector);
  clearCommArray(WorkScalar);
  
  {
    // compute phase gradient and store them in the work vector
    FunctorPhaseGradient functor(this, atom->density, work3);
    functor.setFilter(&filter);

    PairFor(functor, functor.getNumberOfWork());

    // communicate the phase gradient
    comm_variable = WorkVector;
    comm_forward = 3;
    comm->forward_comm_pair(this);
  }
  
  {
    FixISPH_IgnorePhaseGradient *fix = NULL;
    for (int i=0;i<modify->nfix;++i) 
      if (strcmp(modify->fix[i]->style,"isph/ignore/phasegradient") == 0) {
        fix = dynamic_cast<FixISPH_IgnorePhaseGradient*>(modify->fix[i]);
        if (fix != NULL) 
          fix->ignorePhaseGradient(work3, sqrt(getDefaultCutSq()));
      }
  }

  {
    // normalize phase gradient to be phase normal and
    // store magnitude of phase gradients in work array
    FunctorNormalizeVector functor(this, work3, work);
    functor.setFilter(&filter);

    PairFor(functor, functor.getNumberOfWork());

    // communicate the phase gradient
    comm_variable = WorkScalar;
    comm_forward = 1;
    comm->forward_comm_pair(this);
  }

  {
    FunctorCorrectPhaseNormal functor(this, pnd, normal, work3);
    functor.setFilter(&filter);
    
    PairFor(functor, functor.getNumberOfWork());

    // communicate the phase gradient
    comm_variable = WorkVector;
    comm_forward = 3;
    comm->forward_comm_pair(this);
  }

  if (!pb.is_enabled)
    memcpy(&atom->psi[0], &work[0], sizeof(double)*atom->nlocal);
  
  // grad(c)
  if (!pb.is_enabled)
    memcpy(&atom->psigrad[0][0], &work3[0][0], sizeof(double)*3*atom->nlocal);

  // update atom->f with surface tension force  
  {
    FunctorContinuumSurfaceForce functor(this, work3, work, atom->f);
    functor.setFilter(&filter);
    
    PairFor(functor, functor.getNumberOfWork());
  }

  FUNCT_EXIT(comm->me);
}

#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeSurfaceTension_PairwiseForce"
void PairISPH_Corrected::computeSurfaceTension_PairwiseForce() {
  FUNCT_ENTER(comm->me);

  FilterBinary filter;
  filter.setPairYes(Fluid, Fluid);

  // update atom->f with surface tension force  
  {
    FunctorPairwiseForce functor(this, atom->f);
    functor.setFilter(&filter);
    functor._f_sum[0] = 0;
    functor._f_sum[1] = 0;
    functor._f_sum[2] = 0;
    PairFor(functor, functor.getNumberOfWork());
    
    std::cout << "sum of pairwise force = ";
    for (int k=0;k<3;++k)
      std::cout << functor._f_sum[k] << "  ";
    std::cout << "\n";
  }
  
  FUNCT_EXIT(comm->me);
}

// END::ST

// BEGIN::RS
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeForceFromRandomStress"
void PairISPH_Corrected::computeForceFromRandomStress()  {
  FUNCT_ENTER(comm->me);

  const int nmax = atom->nmax;
  const int dim = domain->dimension;
  //const int nlocal = atom->nlocal;

  // communication buffer is packed either 1 or 3; need to match
  double **rstress_x = NULL, **rstress_y = NULL, **rstress_z = NULL; 

  switch (dim) {
  case 3:
    memory->create(rstress_z, nmax, 3, "pair:computeForceFromRandomStress:rstress_z");
  case 2:
    memory->create(rstress_x, nmax, 3, "pair:computeForceFromRandomStress:rstress_x");
    memory->create(rstress_y, nmax, 3, "pair:computeForceFromRandomStress:rstress_y");
  }

  computeRandomStressTensor(rstress_x, rstress_y, rstress_z);

  {
    // generate force field based on the random tensor above
    FunctorRandomStress functor(this,
                                update->dt, 
                                atom->viscosity, atom->density, 
                                atom->f,
                                rstress_x, rstress_y, 
                                (dim > 2 ? rstress_z : NULL));
    
    FilterBinary filter;
    filter.setPairYes(Fluid, Fluid);
    
    functor.setFilter(&filter);
    PairFor(functor, functor.getNumberOfWork());
  }
  
  switch (dim) {
  case 3:
    memory->destroy(rstress_z);
  case 2:
    memory->destroy(rstress_x);
    memory->destroy(rstress_y);
  }
  
  FUNCT_EXIT(comm->me);
}
// END::RS

// BEGIN::RT
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeSoluteTransport"
void PairISPH_Corrected::computeSoluteTransportSpecies(const double dcoeff,
                                                         double *b) {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_solute_transport_species);

  clearCommArray(WorkScalar);
  {
    FunctorSoluteTransport functor(this, 
                                   atom->v, tr.is_particle_fixed,
                                   update->dt, tr.theta,
                                   dcoeff, b, work);
    
    PairFor(functor, functor.getNumberOfWork());
  }

  FUNCT_EXIT(comm->me);
}
// END::RT

// BEGIN::NS - Pressure correction scheme
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeHelmholtz"
void PairISPH_Corrected::computeHelmholtz(double *b) {
  FUNCT_ENTER(comm->me);
  
  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_helmholtz);

  clearCommArray(WorkVector);
  
  // construct matrix A and rhs 
  {
    switch (ns.boundary) {
    case NoBoundaryCond:     
    case HomogeneousNeumann:  
    case ConstExtension:
    case NavierSlip: {
#undef FunctorSymType
#define FunctorSymType(SymType)                                         \
      {                                                                 \
        FunctorIncompNavierStokesHelmholtz_##SymType functor(this,        \
                                                           update->dt, ns.theta, \
                                                           atom->viscosity, atom->density, atom->pressure, \
                                                           atom->f,     \
                                                           b, atom->nlocal, \
                                                           &work3[0][0]); \
        PairFor(functor, functor.getNumberOfWork());                    \
      }
      FunctorSelect;
      break;
    }
    case MorrisHolmes: {
#undef FunctorSymType
#define FunctorSymType(SymType)                                         \
      {                                                                 \
        FunctorIncompNavierStokesHelmholtz_MorrisHolmes_##SymType functor(this, \
                                                                        update->dt, ns.theta, \
                                                                        atom->viscosity, atom->density, atom->pressure, \
                                                                        atom->f, \
                                                                        b, atom->nlocal, \
                                                                        &work3[0][0]); \
        PairFor(functor, functor.getNumberOfWork());                    \
      }
      FunctorSelect;
      break;
    }
    case Dirichlet:
    case MorrisNormal: {
      break;
    }
    }  
  }
  
  // modification of matrix according to applied boundary conditions
  {
    switch (ns.boundary) {
    case NavierSlip: {
      if (ns.beta != 0.0) {
        FunctorBoundaryNavierSlip functor(this, ns.beta, normal, atom->density);
        PairFor(functor, functor.getNumberOfWork());
      }
      break;
    } 
    case Dirichlet: {
      FunctorBoundaryDirichlet functor(this, normal, b, atom->nlocal);
      PairFor(functor, functor.getNumberOfWork());
      break;
    }
    }
  }
  
  FUNCT_EXIT(comm->me);
}

#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeBlockHelmholtz"
void PairISPH_Corrected::computeBlockHelmholtz(double *b) {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_block_helmholtz);

  clearCommArray(WorkVector);
  
  // construct matrix A and rhs

#undef FunctorSymType
#define FunctorSymType(SymType)                                         \
  {                                                                     \
    FunctorIncompNavierStokesBlockHelmholtz_MorrisHolmes_##SymType functor(this, \
                                                                         update->dt, ns.theta, ns.beta, \
                                                                         atom->viscosity, atom->density, atom->pressure, \
                                                                         atom->f, \
                                                                         b, atom->nlocal, \
                                                                         &work3[0][0], normal); \
    PairFor(functor, functor.getNumberOfWork());                        \
  }
  FunctorSelect;

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computePoisson"
void PairISPH_Corrected::computePoisson(double *b) {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_poisson);

  clearCommArray(WorkScalar);
    
  {
    switch (ns.boundary) {
    case NoBoundaryCond:     
    case HomogeneousNeumann:  
    case ConstExtension:
    case NavierSlip: {
#undef FunctorSymType
#define FunctorSymType(SymType)                                         \
      {                                                                 \
        FunctorIncompNavierStokesPoisson_##SymType functor(this,          \
                                                         update->dt, normal, \
                                                         atom->density, \
                                                         vstar, atom->pressure, \
                                                         b, work);      \
        PairFor(functor, functor.getNumberOfWork());                    \
      }
      FunctorSelect;
      break;
    }
    case Dirichlet:
    case MorrisNormal: 
    case MorrisHolmes: {
#undef FunctorSymType
#define FunctorSymType(SymType)                                         \
      {                                                                 \
        FunctorIncompNavierStokesPoisson_MorrisHolmes_##SymType functor(this,     \
                                                                      update->dt, normal, \
                                                                      atom->density, \
                                                                      vstar, atom->pressure, \
                                                                      b, work); \
        PairFor(functor, functor.getNumberOfWork());                    \
      }
      FunctorSelect;
      break;
    }
    }  
  }

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::correctVelocity"
void PairISPH_Corrected::correctVelocity() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_correct_velocity);

  // note that we update on vstar not atom->v
#undef FunctorSymType
#define FunctorSymType(SymType)                 \
  {                                                                     \
    FunctorCorrectVelocity_##SymType functor(this, update->dt, atom->density, dp, vstar); \
    PairFor(functor, functor.getNumberOfWork());                        \
  }                                                                     
  FunctorSelect;
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::correctPressure"
void PairISPH_Corrected::correctPressure() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_correct_pressure);
#undef FunctorSymType
#define FunctorSymType(SymType)                                         \
  {                                                                     \
    FunctorCorrectPressure_##SymType functor(this, atom->pressure, dp, atom->nghost); \
    PairFor(functor, functor.getNumberOfWork());                        \
  }
  FunctorSelect;

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeVelocityDivergence"
void PairISPH_Corrected::computeVelocityDivergence(double **v, double *div) {
  FUNCT_ENTER(comm->me);

  int nlocal = atom->nlocal;

  FilterBinary filter;
  filter.setPairYes(Fluid, All);

  {
    memset(&div[0], 0, sizeof(double)*nlocal);
    switch (ns.boundary) {
    case HomogeneousNeumann:
    case NoBoundaryCond:
    case ConstExtension: {
      FunctorDivergence functor(this, v, 1.0, div);
      
      functor.setFilter(&filter);
      PairFor(functor, functor.getNumberOfWork());
      break;
    }
    case MorrisHolmes: {
      FunctorDivergence_MorrisHolmes functor(this, v, 1.0, div);
      
      functor.setFilter(&filter);
      PairFor(functor, functor.getNumberOfWork());
      break;
    }
    }
  }

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
// END::NS - Pressure correction scheme


// BEGIN::NS - Velocity correction scheme
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::predictAleVelocity"
void PairISPH_Corrected::predictAleVelocity() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_ale_predict_velocity);

  // // store BDF rhs contribution to vstar
  // ale.tdiff->diffVectorVariable(atom->nlocal, domain->dimension, 
  //                               atom->vprev, vstar);
  
  // switch (ns.boundary) {
  // case NoBoundaryCond:
  // case HomogeneousNeumann:
  // case ConstExtension:
  // case NavierSlip: {
  //   FunctorAlePredictVelocity functor(this, 
  //                                     update->dt, ale.tdiff->getGamma(),
  //                                     atom->viscosity, atom->density, 
  //                                     atom->v, atom->xdot, 
  //                                     atom->f,
  //                                     vstar);
  //   PairFor(functor, functor.getNumberOfWork());
  //   break;                                                                                           
  // }                                                                                                  
  // case Dirichlet:                                                                                    
  // case MorrisNormal: 
  // case MorrisHolmes: {                                                                               
  //   FunctorAlePredictVelocity_MorrisHolmes functor(this, 
  //                                                  update->dt, ale.tdiff->getGamma(),
  //                                                  atom->viscosity, atom->density, 
  //                                                  atom->v, atom->xdot, 
  //                                                  atom->f,
  //                                                  vstar);
  //   PairFor(functor, functor.getNumberOfWork());
  //   break;                                                                                           
  // }                                                                                                  
  // }
  
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeAlePoisson"
void PairISPH_Corrected::computeAlePoisson(double *b) {
  FUNCT_ENTER(comm->me);
  
  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_poisson);

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::correctAleVelocity"
void PairISPH_Corrected::correctAleVelocity() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_ale_correct_velocity);
  
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::computeAleHelmholtz"
void PairISPH_Corrected::computeAleHelmholtz(double *b) {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_helmholtz);
  
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
// END::NS - Velocity correction scheme

// BEGIN::Time/Shift
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::advanceTime"
void PairISPH_Corrected::advanceTime() {
  FUNCT_ENTER(comm->me);

  if (use_exact_solution) {
    PairISPH::advanceTime();
  } else {     
    if (ns.is_enabled) {
      if (ale.is_enabled) { 
        
      } else { 
        // particles are moved in final_integrate
        double dt = update->dt; 
#undef FunctorSymType
#define FunctorSymType(SymType)                 \
        {                                                               \
          FunctorAdvanceTimeBegin_##SymType functor(this, dt, atom->v, atom->pressure); \
          PairFor(functor, functor.getNumberOfWork());                  \
        }
        FunctorSelect;
        {
          FunctorAdvanceTimeEnd functor(this, dt, atom->v, atom->pressure, atom->nghost);
          PairFor(functor, functor.getNumberOfWork());
        }
      }
    }
  }
FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::shiftParticles"
void PairISPH_Corrected::shiftParticles(const double shift,
                                        const double shiftcut,
                                        const double solidweight) {
  FUNCT_ENTER(comm->me);
  // - re-neighboring and communication is done in fix/shift
  // - required precomputation, which depends on Corrected/MLS,
  //   is also done in fix/shift

  // compute maximum fluid particle velocity
  int *type  = atom->type;
  int nlocal = atom->nlocal;
  int dim    = domain->dimension;

  double vmax[4] = {}, vmax_all[4] = {}, vsum[2] = {}, vsum_all[2] = {};
  for (int i=0;i<nlocal;++i) {
    if (!(getParticleKind(type[i]) & Fluid))
      continue;

    double magsq = 0.0;
    for (int k=0;k<dim;++k) {
      double val = abs(atom->v[i][k]);
      vmax[k] = max(vmax[k], val);
      magsq += val*val;
    }
    const double mag = sqrt(magsq);
    vmax[3] = max(vmax[3], mag);
    ++vsum[0];
    vsum[1] += mag;
  }
  MPI_Allreduce(vmax, vmax_all, 4, MPI_DOUBLE, MPI_MAX, world);
  MPI_Allreduce(vsum, vsum_all, 2, MPI_DOUBLE, MPI_SUM, world);

  //const double vshift = (vsum_all[0] != 0.0 ? vsum_all[1]/vsum_all[0] : 0.0);
  const double vshift = vmax_all[3];
  
  // compute shift distance based on the current x
  double dt = update->dt;
  {
    // compute shift dr and store it on work3
    double coeff = shift*dt*vshift;
    FunctorComputeShift functor(this, coeff, shiftcut, solidweight, work3);
    PairFor(functor, functor.getNumberOfWork());
  }

  // apply shift
  if (ale.is_enabled) {
    
  } else {
    // shift is applied in final_integrate
    // shift particles and interpolate the velocity and pressure
#undef FunctorSymType
#define FunctorSymType(SymType)   \
    {                                                                   \
      FunctorApplyShift_##SymType functor(this, work3, atom->v, atom->pressure);  \
      PairFor(functor, functor.getNumberOfWork());                      \
    }
    FunctorSelect;
  }
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
// END::Time/Shift 


/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "PairISPH_Corrected::coeff"
void PairISPH_Corrected::coeff(int narg, char **arg) {
  FUNCT_ENTER(comm->me);
  
  if (narg < 4 || narg > 5)
    error->all(FLERR, "Incorrect args for pair_style isph coefficients");
  
  PairISPH::coeff(narg, arg);
  
  // ** kernel setup; Wendland is default
  {
    // type table size
    int ilo, ihi, jlo, jhi;
    force->bounds(arg[0], atom->ntypes, ilo, ihi);
    force->bounds(arg[1], atom->ntypes, jlo, jhi);

    // smoothing length is given from the input file
    double h_one = force->numeric(FLERR,arg[3]);
    double h_min = h_one;
    
    if (narg == 5)
      h_min = force->numeric(FLERR,arg[4]);
    
    // setup kernel
    string str("Kernel Function");
    double cut_one = 0.0, cut_one_sq = 0.0;

    if (g_params.isSublist(str)) {
      auto kernel_param = g_params.sublist(str);
      auto kernel_type = kernel_param.get("type", "Wendland");
      if      (kernel_type == "Wendland")
        kernel = new KernelFuncWendland(domain->dimension);
      else if (kernel_type == "Quintic")
        kernel = new KernelFuncQuintic(domain->dimension);
      else if (kernel_type == "Cubic")
        kernel = new KernelFuncCubic(domain->dimension);
      else
        error->all(FLERR, "Kernel is not in supported list: Wendland, Quintic, and Cubic");  

      cut_one = h_one*kernel_param.get("cut over h", 2.0);
      morris_safe_coeff = kernel_param.get("MorrisSafeCoeff", 0.43301);
    } else {
      kernel = new KernelFuncWendland(domain->dimension);
      cut_one = h_one*2.0;
      morris_safe_coeff = 0.43301; // suggested from morris sqrt(3)/4;
    }
    cut_one_sq = cut_one*cut_one;

    // cut can be different depending on materials;
    // here we use one cut
    int count = 0;
    for (int i=ilo;i<=ihi;++i) {
      for (int j=jlo;j<=jhi;++j) {
        // LAMMPS required
        cutsq[i][j] = cut_one_sq;
        setflag[i][j] = 1;

        // isph specific parameter
        // when particle kinds are different, reduce the smoothing length 
        // so make non-local operator behave like local operator.
        if (getParticleKind(i) == getParticleKind(j))
          h[i][j] = h_one;
        else
          h[i][j] = h_min;; 

        ++count;
      }
    }
    if (count == 0)
      error->all(FLERR,"Incorrect args for pair coefficients");
  }

  grow(1.2);

  FUNCT_EXIT(comm->me);
}

/* ----------------------------------------------------------------------
   communication
   ------------------------------------------------------------------------- */
int PairISPH_Corrected::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc) {  

  int m = PairISPH::pack_forward_comm(n, list, buf, pbc_flag, pbc);

  switch (comm_variable) {
  case NormalVector: {
    for (int ii=0;ii<n;++ii) {
      int i = list[ii];
      buf[m++] = pnd[i];
    }
    break;
  }
  }
  return m;
}

void PairISPH_Corrected::unpack_forward_comm(int n, int first, double *buf) { 

  PairISPH::unpack_forward_comm(n, first, buf);
  
  switch (comm_variable) {
  case NormalVector:
    memcpy(&pnd[first], &buf[n*3], sizeof(double)*n);
    break;
  }
}
