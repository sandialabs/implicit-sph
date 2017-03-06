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

// user
#include <iostream>
#include <fstream>
#include <string>
#include "utils.h"

#include "pair_isph_mls.h"

#include "pair_for.h"

using namespace std;
using namespace LAMMPS_NS;

typedef class PairISPH_MLS Self;  

// common operators
// ================
typedef class FunctorOuterGraph<Self> FunctorGraph;

// math-related operators for "mls" scheme
// ======================================= 
typedef class MLS::FunctorOuterMassMatrix<Self>               FunctorMassMatrix;
typedef class MLS::FunctorOuterMassMatrixCompactPoisson<Self> FunctorMassMatrixCompactPoisson;
typedef class MLS::FunctorOuterNormal<Self>                   FunctorNormal;
    
typedef class MLS::FunctorOuterHelper<Self>                        FunctorHelper;
typedef class MLS::FunctorOuterGradient<Self>                      FunctorGradient;
typedef class MLS::FunctorOuterGradientCompactPoisson<Self>        FunctorGradientCompactPoisson;
typedef class MLS::FunctorOuterDivergence<Self>                    FunctorDivergence;
typedef class MLS::FunctorOuterLaplacian<Self>                     FunctorLaplacian;
typedef class MLS::FunctorOuterLaplacianCompactPoisson<Self>       FunctorLaplacianCompactPoisson;
typedef class MLS::FunctorOuterLaplacianMatrix<Self>               FunctorLaplacianMatrix;
typedef class MLS::FunctorOuterLaplacianMatrixCompactPoisson<Self> FunctorLaplacianMatrixCompactPoisson;
typedef class MLS::FunctorOuterCurl<Self>                          FunctorCurl;
typedef class MLS::FunctorOuterCurlCurl<Self>                      FunctorCurlCurl;

typedef class MLS::FunctorOuterBoundaryNeumann<Self>               FunctorBoundaryNeumann;

// gradient base curl
typedef class FunctorOuterGradientDotOperatorMatrix<Self,
                                                    MLS::FunctorOuterGradientOperator
                                                    > FunctorGradientDotOperatorMatrix;

// no more use this as we have a 2nd order curlcurl; 
// the combination of curlcurl is not much better than 2nd order curlcurl 
// typedef class FunctorOuterCurlCurl<Self,MLS::FunctorOuterGradient> FunctorCurlCurl;
    
// Poisson-Boltzmann NOX-related operators for "mls" scheme
// ========================================================
typedef class FunctorOuterPoissonBoltzmannF<Self,
                                            MLS::FunctorOuterLaplacian
                                            > FunctorPoissonBoltzmannF;
typedef class FunctorOuterPoissonBoltzmannExtraF<Self> FunctorPoissonBoltzmannExtraF;
typedef class FunctorOuterPoissonBoltzmannJacobian<Self,
                                                   MLS::FunctorOuterLaplacianMatrix
                                                   > FunctorPoissonBoltzmannJacobian; 

// Navier-Stokes pressure correction scheme
// ========================================
typedef class FunctorOuterIncompNavierStokesHelmholtz<Self,
                                                      MLS::FunctorOuterLaplacianMatrix,
                                                      MLS::FunctorOuterGradient
                                                      > FunctorIncompNavierStokesHelmholtz;
typedef class FunctorOuterIncompNavierStokesPoisson<Self,
                                                    MLS::FunctorOuterLaplacianMatrix, 
                                                    MLS::FunctorOuterDivergence,
                                                    MLS::FunctorOuterGradientOperator
                                                    > FunctorIncompNavierStokesPoisson;

typedef class FunctorOuterCorrectVelocity<Self,
                                          MLS::FunctorOuterGradient
                                          > FunctorCorrectVelocity; 
typedef class FunctorOuterCorrectPressure<Self> FunctorCorrectPressure;

typedef class FunctorOuterAdvanceTimeBegin<Self,
                                           MLS::FunctorOuterGradient
                                           > FunctorAdvanceTimeBegin;
typedef class FunctorOuterAdvanceTimeEnd<Self> FunctorAdvanceTimeEnd;

typedef class FunctorOuterComputeShift<Self> FunctorComputeShift;
typedef class FunctorOuterApplyShift<Self,MLS::FunctorOuterGradient
                                     > FunctorApplyShift;

// Navier-Stokes velocity correction scheme: ALE
// =============================================
typedef class FunctorOuterAleAdvection<Self,
                                       MLS::FunctorOuterGradientOperator
                                       > FunctorAleAdvection;
typedef class FunctorOuterAleAdvectionMatrix<Self,
                                             MLS::FunctorOuterGradientOperator
                                             > FunctorAleAdvectionMatrix;

// Curcurl operator type
// =====================
//#define USE_GRADIENT_BASE_CURL
#undef  USE_GRADIENT_BASE_CURL

// Gradient based curlcurl
// =======================
#ifdef  USE_GRADIENT_BASE_CURL
typedef class FunctorOuterAlePredictVelocity<Self,
                                             MLS::FunctorOuterGradientOperator,
                                             MLS::FunctorOuterGradient,
                                             MLS::FunctorOuterGradient
                                             > FunctorAlePredictVelocity;
typedef class FunctorOuterAleIncompNavierStokesHelmholtz<Self,
                                                         MLS::FunctorOuterLaplacianMatrix,
                                                         MLS::FunctorOuterGradientOperator,
                                                         MLS::FunctorOuterGradient,
                                                         MLS::FunctorOuterGradient
                                                         > FunctorAleIncompNavierStokesHelmholtz;
#else
// MLS curlcurl
// ============
typedef class FunctorOuterAlePredictVelocityCurlCurl<Self,
                                                     MLS::FunctorOuterGradientOperator,
                                                     MLS::FunctorOuterCurlCurl
                                                     > FunctorAlePredictVelocity;
typedef class FunctorOuterAleIncompNavierStokesHelmholtzCurlCurl<Self,
                                                                 MLS::FunctorOuterLaplacianMatrix,
                                                                 MLS::FunctorOuterGradientOperator,
                                                                 MLS::FunctorOuterCurlCurl
                                                                 > FunctorAleIncompNavierStokesHelmholtz;
#endif

typedef class FunctorOuterAleCorrectVelocity<Self,
                                             MLS::FunctorOuterGradient
                                             > FunctorAleCorrectVelocity;

typedef class FunctorOuterAleIncompNavierStokesPoissonBoundary<Self,
                                                               MLS::FunctorOuterLaplacianMatrix, 
                                                               MLS::FunctorOuterDivergence,
                                                               MLS::FunctorOuterGradientOperator,
                                                               MLS::FunctorOuterBoundaryNeumann
                                                               > FunctorAleIncompNavierStokesPoisson;

typedef class FunctorOuterAleIncompNavierStokesCompactPoissonBoundary<Self,
                                                                      MLS::FunctorOuterLaplacianMatrixCompactPoisson
                                                                      > FunctorAleIncompNavierStokesCompactPoisson;

typedef class FunctorOuterAleTrackParticles<Self> FunctorAleTrackParticles;
typedef class FunctorOuterAleApplyShift<Self> FunctorAleApplyShift;

// NS - Post processing
typedef class FunctorOuterTractionVector<Self,
                                         MLS::FunctorOuterGradient
                                         > FunctorTractionVector;

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::PairISPH_MLS"
PairISPH_MLS::PairISPH_MLS(LAMMPS *lmp) 
  : PairISPH(lmp), 
    M(NULL), Ma(NULL), Mb(NULL), np(0) {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor initalizationTimerMonitor(*g_timer_init);

  discretization = MLS;

  nmax = 0;

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::~PairISPH_MLS"
PairISPH_MLS::~PairISPH_MLS() {
  FUNCT_ENTER(comm->me);

  delete kernel;
  delete basis;

  memory->destroy(Ma);
  memory->destroy(Mb);

  FUNCT_EXIT(comm->me);
}

// BEGIN::COMMON
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::grow"
void PairISPH_MLS::grow(double extend) {
  FUNCT_ENTER(comm->me);

  if (atom->nmax > nmax) {
    // nmax is updated here
    PairISPH::grow(extend);

    const int ndof   = basis->ndof(ISPH_MLS_MAX_ORDER);
    const int ndofpl = ndof + 1;  // lagrange multiplier

    memory->grow(Ma, nmax, ndof*ndof,     "pair:Ma");
    memory->grow(Mb, nmax, ndofpl*ndofpl, "pair:Mb");
    
  }
  M = Ma;

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computePre"
void PairISPH_MLS::computePre(const int mask) {
  FUNCT_ENTER(comm->me);

  computeMassMatrix();

  if (boundary_particle != NoParticle)
   computeNormals();
  
  if (cp.is_enabled) 
    computeMassMatrixCompactPoisson();

  FUNCT_EXIT(comm->me);
}  

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computeMassMatrix"
void PairISPH_MLS::computeMassMatrix() {
  FUNCT_ENTER(comm->me);
  
  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_mass_matrix);

  // compute for all interacting with all neighbors
  FunctorMassMatrix functor(this);
  functor.usePseudoInverse(cp.is_enabled);

  PairFor(functor, functor.getNumberOfWork());

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computeMassMatrixCompactPoisson"
void PairISPH_MLS::computeMassMatrixCompactPoisson() {
  FUNCT_ENTER(comm->me);
  
  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_mass_matrix_compact_poisson);

  enterCompactPoisson();

  // compute for all interacting with all neighbors
  FunctorMassMatrixCompactPoisson functor(this, normal);
  functor.usePseudoInverse(true);

  PairFor(functor, functor.getNumberOfWork());

  exitCompactPoisson();

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computeNormals"
void PairISPH_MLS::computeNormals() {
  FUNCT_ENTER(comm->me);
  
  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_normals);

  // orientaiton is set from solid to fluid particles
  double orient[MaxParticleKind] = {};
  orient[Fluid]    = -1.0;
  orient[Solid]    =  1.0;
  orient[Boundary] =  1.0;

  // conventional MLS boundary particles
  switch (boundary_particle) {
  case Boundary: {
    // compute normal for boundary particles
    FilterExclusiveOr filter;
    filter.setPairYes(All);

    FunctorNormal functor(this, &orient[0], normal);

    functor.setFilter(&filter);
    PairFor(functor, functor.getNumberOfWork());
    break;
  }
  case Solid: {
    // compute normal for solid and fluid; here boundary particle does not exist
    FilterExclusiveOr filter;
    filter.setPairYes(All);

    FunctorNormal functor(this, &orient[0], normal);

    functor.setFilter(&filter);
    PairFor(functor, functor.getNumberOfWork());
    break;
  }
  }

  comm_variable = NormalVector;
  comm_forward = 3;
  comm->forward_comm_pair(this);

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::enterCompactPoisson"
void PairISPH_MLS::enterCompactPoisson() {
  FUNCT_ENTER(comm->me);
  M = Mb;
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::exitCompactPoisson"
void PairISPH_MLS::exitCompactPoisson() {
  FUNCT_ENTER(comm->me);
  M = Ma; 
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::isCompactPoisson"
bool PairISPH_MLS::isCompactPoisson() {
  return (M == Mb);
}

/* ---------------------------------------------------------------------- */
// END::COMMON

// BEGIN::PBE NOX
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computeF"
bool PairISPH_MLS::computeF(const Epetra_Vector &psitmp,
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
#define __FUNCT__ "PairISPH_MLS::computeJacobian"
bool PairISPH_MLS::computeJacobian(const Epetra_Vector &psitmp,
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
  }

  FUNCT_EXIT(comm->me);
  return true; 
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computePreconditioner"
bool PairISPH_MLS::computePreconditioner(const Epetra_Vector &x,
                                         Epetra_Operator &M,
                                         Teuchos::ParameterList *param) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "We use Jacobian only at this moment: isph");  
  FUNCT_EXIT(comm->me);
  return false;
}  

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computePsiGradient"
void PairISPH_MLS::computePsiGradient() {
  FUNCT_ENTER(comm->me);

  // cleanup
  memset(&atom->psigrad[0][0], 0, sizeof(double)*3*atom->nlocal);

  FunctorGradient functor(this, atom->psi, 1.0, atom->psigrad);
  
  FilterBinary filter;
  filter.setPairYes(Fluid | Boundary, 
                    Fluid | Boundary);
  functor.setFilter(&filter);
  
  PairFor(functor, functor.getNumberOfWork());
  
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
// END::PBE NOX

// BEGIN::NS - Pressure correction scheme
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computeHelmholtz"
void PairISPH_MLS::computeHelmholtz(double *b) {
  FUNCT_ENTER(comm->me);
  
  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_helmholtz);
  
  // construct matrix A and rhs 
  {
    // no morris type boundary condition
    switch (ns.boundary) {
    case NoBoundaryCond:     
    case HomogeneousNeumann:  
    case ConstExtension: {
      FunctorIncompNavierStokesHelmholtz functor(this, 
                                                 update->dt, ns.theta, 
                                                 atom->viscosity, atom->density, atom->pressure, 
                                                 atom->f,
                                                 b, atom->nlocal, 
                                                 &work3[0][0]);
      PairFor(functor, functor.getNumberOfWork());
      break;
    }
    default:
      error->all(FLERR, "MLS does not support NavierSlip and Morris type boundary conditions, use ConstExtension");
      break;
    }  
  }

  // boundary condition
  // no boundary now 

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computePoisson"
void PairISPH_MLS::computePoisson(double *b) {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_poisson);

  FunctorIncompNavierStokesPoisson functor(this,
                                           update->dt, normal,
                                           atom->density, 
                                           vstar, atom->pressure,
                                           b, work);
  PairFor(functor, functor.getNumberOfWork());

  //A.crs->Print(cout);

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::correctVelocity"
void PairISPH_MLS::correctVelocity() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_correct_velocity);

  // note that we update on vstar not atom->v
  FunctorCorrectVelocity functor(this, update->dt, atom->density, dp, vstar);
  PairFor(functor, functor.getNumberOfWork());

  FUNCT_EXIT(comm->me);
}
   
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::correctPressure"
void PairISPH_MLS::correctPressure() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_correct_pressure);
  
  FunctorCorrectPressure functor(this, atom->pressure, dp, atom->nghost);
  PairFor(functor, functor.getNumberOfWork());

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
// END::NS - Pressure correction scheme

// BEGIN::NS - Velocity correction scheme
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::predictAleVelocity"
void PairISPH_MLS::predictAleVelocity() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_ale_predict_velocity);

  // store BDF rhs contribution to vstar
  ale.tdiff->diffVectorVariable(atom->nlocal, domain->dimension, 
                                atom->vprev, vstar);

  // update vstar
  FunctorAlePredictVelocity functor(this, 
                                    update->dt, ale.tdiff->getGamma(),
                                    atom->viscosity, atom->density, 
                                    atom->v, atom->xdot, 
                                    atom->f,
                                    vstar);
  PairFor(functor, functor.getNumberOfWork());

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computeAlePoisson"
void PairISPH_MLS::computeAlePoisson(double *b) {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_poisson);

  ale.tdiff->diffVectorVariable(atom->nlocal, domain->dimension,
                                atom->vprev, work3);

  comm_variable = WorkVector;
  comm_forward = 3;
  comm->forward_comm_pair(this);

  if (cp.is_enabled) {
    const int nmax = atom->nmax;

    double *f = NULL, *g = NULL;
    memory->grow(f, nmax, "pair:mls:ale:cp:f");
    memory->grow(g, nmax, "pair:mls:ale:cp:g");

    memset(f, 0, sizeof(double)*nmax);
    memset(g, 0, sizeof(double)*nmax);

    const double gamma_over_dt = ale.tdiff->getGamma()/update->dt;

    { // compute f = \gamma_over_dt \div v
      FilterBinary filter;
      filter.setPairYes(All, All);

      FunctorDivergence functor(this, vstar, -gamma_over_dt, f); 
      functor.setFilter(&filter);
      
      PairFor(functor, functor.getNumberOfWork());

      temp = f;
      comm_variable = TempScalar;
      comm_forward = 1;
      comm->forward_comm_pair(this);
    }
    { // compute g = \gamma_over_dt ((w[i] - v[i]) \cdot n[i])
      // - w is given velocity in the boundary region
      // - v is predicted velocity
      // - n normal vector at i
      FilterBinary filter;
      filter.setPairYes(Boundary);

      FunctorBoundaryNeumann functor(this, normal, work3, vstar, gamma_over_dt, g);
      functor.setFilter(&filter);
      
      PairFor(functor, functor.getNumberOfWork());

      temp = g;
      comm_variable = TempScalar;
      comm_forward = 1;
      comm->forward_comm_pair(this);
    }
    { // compute poisson matrix and modify the entrees according to the boundary condition 
      enterCompactPoisson();
      FunctorAleIncompNavierStokesCompactPoisson functor(this,
                                                         atom->density,
                                                         f, normal, g,
                                                         b);

      PairFor(functor, functor.getNumberOfWork());
      exitCompactPoisson();
    }

    memory->destroy(f);
    memory->destroy(g);
  } else {
    FunctorAleIncompNavierStokesPoisson functor(this,
                                                update->dt, ale.tdiff->getGamma(),
                                                normal,
                                                atom->density, 
                                                work3, vstar,
                                                b);
    PairFor(functor, functor.getNumberOfWork());
  }

  //A.crs->Print(cout);

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::correctAleVelocity"
void PairISPH_MLS::correctAleVelocity() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_ale_correct_velocity);

  FunctorAleCorrectVelocity functor(this, 
                                    update->dt, ale.tdiff->getGamma(), 
                                    atom->density, 
                                    atom->pressure, 
                                    vstar);
  PairFor(functor, functor.getNumberOfWork());
  
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computeAleHelmholtz"
void PairISPH_MLS::computeAleHelmholtz(double *b) {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_helmholtz);
  
  // TODO:: boundary condition
  {
    FunctorAleIncompNavierStokesHelmholtz functor(this, 
                                                  update->dt, ale.tdiff->getGamma(), 
                                                  atom->viscosity, atom->density, 
                                                  vstar, 
                                                  atom->v, atom->xdot,
                                                  b, atom->nlocal);
    PairFor(functor, functor.getNumberOfWork());
  }

  //A.crs->Print(cout);

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computeVelocityCurl"
void PairISPH_MLS::computeVelocityCurl(double **v, double **curl) {
  FUNCT_ENTER(comm->me);
  int nlocal = atom->nlocal;
  
  FilterBinary filter;
  filter.setPairYes(Fluid | Boundary,
                    Fluid | Boundary);
  
  memset(&curl[0][0], 0, sizeof(double)*3*nlocal);
  FunctorCurl functor(this, v, 1.0, curl);
  
  functor.setFilter(&filter);
  PairFor(functor, functor.getNumberOfWork());

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computeVelocityDivergence"
void PairISPH_MLS::computeVelocityDivergence(double **v, double *div) {
  int nlocal = atom->nlocal;
  
  FilterBinary filter;
  filter.setPairYes(Fluid | Boundary,
                    Fluid | Boundary);
  
  memset(&div[0], 0, sizeof(double)*nlocal);
  FunctorDivergence functor(this, v, 1.0, div);
  
  functor.setFilter(&filter);
  PairFor(functor, functor.getNumberOfWork());
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computeTractionVector"
void PairISPH_MLS::computeTractionVector(double *p, double **v, double **traction) {
  int nlocal = atom->nlocal;
  
  FilterBinary filter;
  filter.setPairYes(Boundary,
                    Fluid | Boundary);
  
  memset(&traction[0][0], 0, sizeof(double)*3*nlocal);
  FunctorTractionVector functor(this, 
                                atom->density,
                                atom->viscosity,
                                p, v, 
                                normal, traction);
  
  functor.setFilter(&filter);
  PairFor(functor, functor.getNumberOfWork());
}

/* ---------------------------------------------------------------------- */
// END::NS - Velocity correction scheme

// BEGIN::Time/Shift
/* ---------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PairISPH_MLS::computeUnitTests"
void PairISPH_MLS::computeUnitTests() {
  FUNCT_ENTER(comm->me);

  MLS::TestSuite ts(this);
  
  ts.doUnitTests();

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
// END::NS - Velocity correction scheme

// BEGIN::Time/Shift
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::advanceTime"
void PairISPH_MLS::advanceTime() {
  FUNCT_ENTER(comm->me);

  if (use_exact_solution) {
    PairISPH::advanceTime();
  } else {
    if (ale.is_enabled) { 
    
      // particles are moved in initial_integrate;
      // note that Gc, Lc are not computed yet at this moment
      // which also means that we do not use any mathematical operators here
      double dt     = update->dt; 
      int ntimestep = update->ntimestep;
      
      int nlocal = atom->nlocal;
      int dim    = domain->dimension;
      
      // update BDF vectors
      ale.tdiff->updateTimestepBegin(ntimestep, dt);
      ale.tdiff->setCurrentVectorVariable(nlocal, dim, atom->v, atom->vprev);
      ale.tdiff->setCurrentVectorVariable(nlocal, dim, atom->x, atom->xprev);
      ale.tdiff->updateTimestepEnd();
      
      // update v, xdot and x 
      
      // atom->v := extrpolated v
      ale.tdiff->extrapolateVectorVariable(nlocal, dim, atom->vprev, atom->v);
      
      // xdot := atom->v
      memcpy(&atom->xdot[0][0], &atom->v[0][0], sizeof(double)*3*nlocal);
      
      // copy the current position vector to workspace (backup)
      memcpy(&work3[0][0], &atom->x[0][0], sizeof(double)*3*nlocal);    
      
      // recover relative varialbes; the current x becomes the reference
      ale.tdiff->recoverRelativeVectorVariable(nlocal, dim, atom->xprev);
      
      // store BDF rhs contribution in atom->x
      ale.tdiff->diffVectorVariable(nlocal, dim, atom->xprev, atom->x);
      
      // update x := x + dt * xdot
      FunctorAleTrackParticles functor(this,                                                    
                                       dt, ale.tdiff->getGamma(), 
                                       atom->xdot);                                             
      PairFor(functor, functor.getNumberOfWork());                   
      
      // update the current relative distance in xprev
      ale.tdiff->updateRelativeVectorVariable(nlocal, dim, atom->x, work3, atom->xprev);
      
    } else { 
      // particles are moved in final_integrate
      double dt = update->dt; 
      {
        FunctorAdvanceTimeBegin functor(this, dt, atom->v, atom->pressure);
        PairFor(functor, functor.getNumberOfWork());           
      }
      {
        FunctorAdvanceTimeEnd functor(this, dt, atom->v, atom->pressure, atom->nghost);
        PairFor(functor, functor.getNumberOfWork());
      }
    }
  }

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH_MLS::shiftParticles"
void PairISPH_MLS::shiftParticles(const double shift,
                                  const double shiftcut,
                                  const double boundaryweight) {
  // - re-neighboring and communication is done in fix/shift
  // - required precomputation, which depends on Corrected/MLS,
  //   is also done in fix/shift

  // compute maximum fluid particle velocity
  int *type  = atom->type;
  int nlocal = atom->nlocal;
  int dim    = domain->dimension;

  double vmax[4] = {}, vmax_all[4] = {}, vsum[2] = {}, vsum_all[2] = {};
  for (int i=0;i<nlocal;++i) {
    if (getParticleKind(type[i]) != Fluid)
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
    FunctorComputeShift functor(this, coeff, shiftcut, boundaryweight, work3);
    PairFor(functor, functor.getNumberOfWork());
  }

  // apply shift
  if (ale.is_enabled) {
    // add the relative distance on the current xdr
    ale.tdiff->modifyCurrentVectorVariable('-', nlocal, dim, work3, atom->xprev);

    // shift is applied in initial_integrate
    // velocity and pressure are already either updated or extrapolated
    FunctorAleApplyShift functor(this, 
                                 dt, ale.tdiff->getGamma(),
                                 work3,
                                 atom->xdot);
    PairFor(functor, functor.getNumberOfWork());
  } else {
    // shift is applied in final_integrate
    // shift particles and interpolate the velocity and pressure
    FunctorApplyShift functor(this, work3, atom->v, atom->pressure);
    PairFor(functor, functor.getNumberOfWork());
  }

}
/* ---------------------------------------------------------------------- */
// END::Time/Shift 


/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "PairISPH_MLS::coeff"
void PairISPH_MLS::coeff(int narg, char **arg) {
  FUNCT_ENTER(comm->me);

  if (narg < 4 || narg > 5)
    error->all(FLERR, "Incorrect args for pair_style isph coefficients");

  PairISPH::coeff(narg, arg);

  // ** kernel cut length setup
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
    
    string str("Kernel Function");
    double cut_one = 0.0, cut_one_sq = 0.0;

    if (g_params.isSublist(str)) {
      auto kernel_param = g_params.sublist(str);
      cut_one = h_one*kernel_param.get("cut over h", 1.0);
      np = kernel_param.get("MLS Order", 2);
      is_interpolation_enabled = (kernel_param.get("Use Interpolation", "Enabled") == "Enabled");
    } else {
      cut_one = h_one*1.0;
      np = 2;
      is_interpolation_enabled = true;
    }
    cut_one_sq = cut_one*cut_one;

    // create basis and kernel functions; this shoud go coeff
    kernel = new KernelFuncMLS();
    basis  = new ScaledTaylorMonomial(domain->dimension, is_interpolation_enabled);
    
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

  // ** ns setup
  if (ns.is_enabled || ut.is_enabled) {
    string str("Incompressible Navier Stokes");
    if (g_params.isSublist(str)) {
      auto ns_param = g_params.sublist(str);

      string substr("Compact Poisson");
      if (ns_param.isSublist(substr)) {
        cp.is_enabled = true;
        
        auto cp_param = ns_param.sublist(substr);
        
        // actual control parameters are computed 
        // tau*h^dim for interiorand tau*h^{dim-1} for boundary
        cp.tau_interior = cp_param.get("tau:interior", 0.01);
        cp.tau_boundary = cp_param.get("tau:boundary",  0.01);
      } else {
        cp.is_enabled = false;
        cp.tau_interior = 0.0;
        cp.tau_boundary = 0.0;
      }
    }
  }
  
  grow(1.2);

  FUNCT_EXIT(comm->me);
}

