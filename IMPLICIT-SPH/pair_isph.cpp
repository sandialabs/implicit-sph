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
#include "bond.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"

#include "modify.h"
#include "compute.h"
#include "random_mars.h"

// user
#include <iostream>
#include <fstream>
#include <string>
#include "utils.h"

#include "bond_isph.h"

#include "time_bdf.h"

#include "precond.h"
#include "precond_ifpack.h"
#include "precond_ml.h"

#include "solver_nox.h"
#include "solver_nox_impl.h"
#include "solver_nox_aztecOO.h"
#include "solver_nox_stratimikos.h"

#include "solver_lin.h"
#include "solver_lin_belos.h"

#include "pair_isph.h"

#include "pair_for.h"

using namespace std;
using namespace LAMMPS_NS;

typedef class PairISPH Self;
typedef class FunctorOuterGraph<Self> FunctorGraph;
typedef class FunctorOuterGraphBoundary<Self> FunctorGraphBoundary;
typedef class FunctorOuterElectrostaticForce<Self> FunctorElectrostaticForce;

typedef class FunctorOuterAdvanceTimeExactSolutionBegin<Self> FunctorAdvanceTimeExactSolutionBegin;
typedef class FunctorOuterAdvanceTimeExactSolutionEnd<Self>   FunctorAdvanceTimeExactSolutionEnd;

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::PairISPH"
PairISPH::PairISPH(LAMMPS *lmp)
  : Pair(lmp),
    normal(NULL),
    vstar(NULL), dp(NULL),
    work(NULL), work3(NULL),
    random(NULL) {
  FUNCT_ENTER(comm->me);

  restartinfo = 0;

  // coeff set these parameters from xml file
  use_exact_solution = false;

  pb.is_enabled = false;
  st.is_enabled = false;
  ns.is_enabled = false;
  ale.is_enabled = false;
  ut.is_enabled = false;

  use_part = false;

  is_once = 0;
  nmax = 0;

  resetParticleKindMap();

  // do not grow here

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::~PairISPH"
PairISPH::~PairISPH() {
  FUNCT_ENTER(comm->me);

  if (allocated) {
    memory->destroy(pinfo);
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(h);

    memory->destroy(st.pf.s);
  }

  delete nl_solver;
  delete li_solver;
  delete prec;

  if (ale.is_enabled)
    delete ale.tdiff;

  if (st.is_enabled) {
    switch (st.model) {
    case ContinuumSurfaceForce: delete st.csf.color; break;
    case PairwiseForce: delete st.pf.force; break;
    }
  }

  memory->destroy(vstar);
  memory->destroy(dp);

  memory->destroy(normal);

  memory->destroy(work3);
  memory->destroy(work);

  if (random != NULL)
    delete random;

  FUNCT_EXIT(comm->me);
}

// BEGIN::COMMON
/* ---------------------------------------------------------------------- */
int PairISPH::getDiscretizationInfo() const {
  return discretization;
}
/* ---------------------------------------------------------------------- */
void PairISPH::resetParticleKindMap() {
  for (int i=0;i<MaxParticleKind;++i)
    particle_kind_map[i] = i;
}
/* ---------------------------------------------------------------------- */
void PairISPH::mapParticleKind(const int from, const int to) {
  particle_kind_map[from] = to;
}
/* ---------------------------------------------------------------------- */
int PairISPH::getParticleKind(const int itype) const { //, const int i) const {
  return particle_kind_map[pinfo[0][itype]];
}
/* ---------------------------------------------------------------------- */
bool PairISPH::isParticleFixed(const int itype) const {
  return pinfo[1][itype];
}
/* ---------------------------------------------------------------------- */
int PairISPH::getParticlePhase(const int itype) const {
  return pinfo[2][itype];
}
/* ---------------------------------------------------------------------- */
bool PairISPH::isAleEnabled() const {
  return ale.is_enabled;
}
/* ---------------------------------------------------------------------- */
double PairISPH::getDefaultCutSq() const {
  return cutsq[1][1];
}
/* ---------------------------------------------------------------------- */
double** PairISPH::getVstar() const {
  return vstar;
}
/* ---------------------------------------------------------------------- */
void PairISPH::setUsePartToComputeNormal(const bool flag) {
  use_part = flag;
}
/* ---------------------------------------------------------------------- */
void PairISPH::clearCommArray(const int comm_variable) {
  const int nmax = atom->nmax;

  double *ptr = NULL;
  int stride = 0;

  switch (comm_variable) {
  case Vfrac: {
    ptr = &atom->vfrac[0];
    stride = 1;
    break;
  }
  case NormalVector: {
    ptr = &normal[0][0];
    stride = 3;
    break;
  }
  case Velocity: {
    ptr = &atom->v[0][0];
    stride = 3;
    break;
  }
  case Pressure: {
    ptr = &atom->pressure[0];
    stride = 1;
    break;
  }
  case Vstar: {
    ptr = &vstar[0][0];
    stride = 3;
    break;
  }
  case DeltaP: {
    ptr = &dp[0];
    stride = 1;
    break;
  }
  case Psi: {
    ptr = &atom->psi[0];
    stride = 1;
    break;
  }
  case WorkScalar: {
    ptr = &work[0];
    stride = 1;
    break;
  }
  case WorkVector: {
    ptr = &work3[0][0];
    stride = 3;
    break;
  }
  case TempScalar: {
    ptr = &temp[0];
    stride = 1;
    break;
  }
  case TempVector: {
    ptr = &temp3[0][0];
    stride = 3;
    break;
  }
  default:
    error->all(FLERR, "Comm: packing empty buffer");
    break;
  }
  memset(ptr, 0, sizeof(double)*nmax*stride);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::printParticles"
void PairISPH::printParticles(ostream &s, string name) {
  FUNCT_ENTER(comm->me);

  double *vfrac = atom->vfrac;
  double **x = atom->x;
  double **v = atom->v;
  double  *p = atom->pressure;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  if (comm->me == 0) {
    s << scientific << name << endl;
    for (int i=0;i<(nlocal+nghost);++i) {
      if (i == nlocal)
        s << " ------------ ghost particles ----------------" << endl;

      if (i<4 || (i > nlocal && i < nlocal +4)) {
        s << "i = " << i
          <<"   tag = " << tag[i]
          << "  vfrac = " << vfrac[i]
          << "  x = " << x[i][0] << "  " << x[i][1]
          << "  v = " << v[i][0] << "  " << v[i][1]
          << "  vstar = " << vstar[i][0] << "  " << vstar[i][1]
          << "  p = " << p[i]
          << "  dp = " << dp[i]
          << endl;
      }
    }
  }

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::grow"
void PairISPH::grow(double extend) {
  FUNCT_ENTER(comm->me);

  if (atom->nmax > nmax) {
    nmax = static_cast<int>(atom->nmax * extend);

    memory->grow(vstar,    nmax,       3, "pair:vstar");
    memory->grow(dp,       nmax,          "pair:dp");

    memory->grow(normal,   nmax,       3, "pair:normal");

    memory->grow(work3,    nmax,       3, "pair:work3");
    memory->grow(work,     nmax,          "pair:work");

  }
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::initializeSolvers"
void PairISPH::initializeSolvers() {
  FUNCT_ENTER(comm->me);

  {
    string str("Preconditioner");
    if (g_params.isSublist(str)) {
      auto prec_param = g_params.sublist(str);

      // create and load the default paramlist
      auto prec_package = prec_param.get("Precond Package", "ML");
      if      (prec_package == "Ifpack")
        prec = new PrecondWrapper_Ifpack(world);
      else if (prec_package == "ML")
        prec = new PrecondWrapper_ML(world);
      else
        error->all(FLERR, "Preconditioner is not in supported list: Ifpack, ML");
      auto param = prec->setParameters();
      //writeParameterListToXmlFile(*param, "prec.xml");

      // update the param if file exists
      string prec_file = prec_param.get("XML File","none");
      if (prec_file != "none") {
        ifstream file(prec_file);
        if (file.good()) {
          Teuchos::updateParametersFromXmlFile(prec_file, Teuchos::inoutArg(*param));
        } else {
          error->all(FLERR, "Preconditioner XML file does not exist");
        }
      }

      // customize the param with remaind params in g_params
      for (auto iter=prec_param.begin();iter!=prec_param.end();++iter) {
        string name = prec_param.name(iter);
        if (name != "Precond Package" &&
            name != "XML File")
          param->setEntry(name, prec_param.getEntry(name));
      }

      // print the updated parameter set
      if (comm->me == 0)
        cout << " - Preconditioner parameters for linear solver - \n\n" << *param << endl;

    } else {
      prec = new PrecondWrapper_ML(world);
      prec->setParameters();
    }
  }

  {
    // create a linear solver
    li_solver = new SolverLinear(world);

    // load default parameters
    li_solver->setParameters();
  }
  {
    //nl_solver = new SolverNonlinear(SolverNonlinear::MatrixFree, *this, world);
    nl_solver = new SolverNonlinear(SolverNonlinear::Analytic, *this, world);

    // load default parameters
    nl_solver->setParameters();

    // convergence test
    nl_solver->setConvergenceTests();
  }

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeGraph"
void PairISPH::computeGraph() {
  FUNCT_ENTER(comm->me);
  {
    FunctorGraph functor(this, nodalmap, &tags_in_cut);
    PairFor(functor, functor.getNumberOfWork());
  }
  if (ns.is_block_helmholtz_enabled) {
    FunctorGraphBoundary functor(this, nodalmap, normal, &tags_in_boundary);
    PairFor(functor, functor.getNumberOfWork());
  }
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computePre"
void PairISPH::computePre(const int mask) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed");
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "PairISPH::computeNormals"
void PairISPH::computeNormals() {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed");
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "PairISPH::computeZeroMeanPressure"
double PairISPH::computeZeroMeanPressure(double *f, bool update) {
  FUNCT_ENTER(comm->me);

  int *ilist = list->ilist;
  int *type = atom->type;

  const int inum = list->inum;

  // if singular poisson pressure field exist for a fluid particle only.
  int nlocal = 0, nglobal = 0;
  double mysum = 0.0;

  const bool is_singular_poisson = (ns.singular_poisson != NotSingularPoisson);
  for (int ii=0;ii<inum;++ii) {
    const int i = ilist[ii];
    const int itype = type[i];
    const int ikind = getParticleKind(itype);

    if (ikind == Solid) {
      // clean noisy pressure on solid particles
      f[i] = 0.0;
    } else {
      // fluid pressure is only considered to have zero mean
      mysum += f[i];
      ++nlocal;
    }
  }

  double mean_val = 0.0;
  MPI_Allreduce(&mysum,  &mean_val, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&nlocal, &nglobal,  1, MPI_INT,    MPI_SUM, world);

  mean_val = mean_val/nglobal;

  const int nsize = (atom->nlocal + atom->nghost);
  if (update)
    for (int i=0;i<nsize;++i)
      f[i] -= mean_val*(getParticleKind(type[i]) != Solid);

  FUNCT_EXIT(comm->me);

  return mean_val;
}

/* ---------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "PairISPH::refreshParticles"
void PairISPH::refreshParticles() {
  FUNCT_ENTER(comm->me);

  if (!neighbor->dist_check && (neighbor->delay <= 1.0)) {
    domain->pbc();
    domain->reset_box();
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();
    comm->exchange();
    if (atom->sortfreq > 0) atom->sort();
    comm->borders();
    neighbor->build();
  } else {
    comm->forward_comm();
  }

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "PairISPH::modifySingularMatrix"
void PairISPH::modifySingularMatrix(const int row, double &diag, double &b) {
  FUNCT_ENTER(comm->me);

  // this should be called only once
  switch (ns.singular_poisson) {
  case NotSingularPoisson:
  case NullSpace:
    // do nothing
    break;
  case PinZero: {
    int ncol;
    double *values;

    A.crs->ExtractGlobalRowView(row, ncol, values);
    memset(values, 0, sizeof(double)*ncol);

    diag = -1.0;
    b = 0.0;
    break;
  }
  case DoubleDiag: {
    diag *= 1.5;
    break;
  }
  }

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
// END::COMMON

// BEGIN::PBE NOX
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeF"
bool PairISPH::computeF(const Epetra_Vector &psitmp,
                        Epetra_Vector &f,
                        const FillType flag) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
  return false;
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeJacobian"
bool PairISPH::computeJacobian(const Epetra_Vector &psitmp,
                               Epetra_Operator &J) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
  return false;
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computePreconditioner"
bool PairISPH::computePreconditioner(const Epetra_Vector &x,
                                     Epetra_Operator &M,
                                     Teuchos::ParameterList *param) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
  return false;
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computePsiGradient"
void PairISPH::computePsiGradient() {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PairISPH::computePoissonBoltzmann"
void PairISPH::computePoissonBoltzmann() {

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_poisson_boltzmann);

  // set nodal map on NOX
  nl_solver->setNodalMap(nodalmap);

  // wrap the potential field
  nl_solver->createSolutionVector(atom->psi);

  // solution fields are initialized
  // use the previous solutions as initial guess
  if (is_once == 0)
    nl_solver->setInitialSolution(SolverNonlinear::Random);

  // set jacobian matrix
  nl_solver->setJacobianMatrix(A.crs);

  // run newton solver;
  // at the beginning of each newton iteration, Psi should be communicated
  nl_solver->solveProblem("PoissonBoltzmann");

  // communicate psi
  comm_variable = Psi;
  comm_forward = 1;
  comm->forward_comm_pair(this);

  // compute electric field
  computePsiGradient();

  // invalidate A
  A.is_filled = 0;
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeAppliedElectricPotential"
void PairISPH::computeAppliedElectricPotential(double *b) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computePhiGradient"
void PairISPH::computePhiGradient() {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeAppliedElectricField"
void PairISPH::computeAppliedElectricField() {
  FUNCT_ENTER(comm->me);

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_applied_electric_field);

  const int nlocal = atom->nlocal;

  li_solver->setNodalMap(nodalmap);
  li_solver->setMatrix(A.crs);
  prec->setMatrix(A.crs);

  // Laplace equation
  {
    // wrap/create [nlocal x 1] multivector for x and b
    li_solver->createSolutionMultiVector(atom->phi, nlocal, 1);
    li_solver->createLoadMultiVector(NULL, nlocal, 1);

    // b = x
    auto x = li_solver->getSolutionMultiVector();
    auto b = li_solver->getLoadMultiVector();
    *b = *x;

    // build system of equations
    computeAppliedElectricPotential(b->Values());
  }

  // solve the Laplace equation
  {
    Teuchos::TimeMonitor evalTimerMonitor(*g_timer_solve_applied_electric_potential);
    li_solver->solveProblem(prec, "AppliedElectricPotential");
  }

  temp = atom->phi;
  comm_variable = TempScalar;
  comm_forward = 1;
  comm->forward_comm_pair(this);
  temp = NULL;

  // invalidate A
  A.is_filled = 0;

  // compute the applied electric field
  computePhiGradient();

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */

#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeElectrostaticForce"
void PairISPH::computeElectrostaticForce() {
  FUNCT_ENTER(comm->me);

  FunctorElectrostaticForce functor(this,
                                    (ae.is_enabled ? atom->phigrad : NULL),
                                    atom->psi,
                                    atom->psigrad,
                                    atom->f);
  PairFor(functor, functor.getNumberOfWork());

  FUNCT_EXIT(comm->me);
}

// END::PBE NOX

// BEGIN::ST - CSF model for multiphase flow
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeSurfaceTension"
void PairISPH::computeSurfaceTension() {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}

// END::ST - CSF model for multiphase flow

// BEGIN::RS - Random stress
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeRandomStressTensor"
void PairISPH::computeRandomStressTensor(double **rstress_x,
                                         double **rstress_y,
                                         double **rstress_z) {
  FUNCT_ENTER(comm->me);

  const int dim = domain->dimension;
  const int nlocal = atom->nlocal;

  {
    // this local filter only yes for fluid particles
    FilterBinary filter;
    filter.setPairYes(Fluid);

    // memory allocation for random stress tensor
    // should be done (memory->grow) before this
    int *type = atom->type;
    int *ilist = list->ilist;
    double rs_tensor_temp[3][3],rs_tensor[3][3];  // temporary work space for tensor
    double rs_tensor_trace;

    // let's use ilist consistently
    for (int ii=0;ii<nlocal;++ii) {
      const int i = ilist[ii];
      if (filter.yes(getParticleKind(type[i]))) {
        // generate a random tensor
        rs_tensor_trace = 0.0;   
        for (int k2=0;k2<dim;++k2) 
          for (int k1=0;k1<dim;++k1) 
            rs_tensor_temp[k2][k1] = random->gaussian(); 
          
        for (int k2=0;k2<dim;++k2) 
          for (int k1=0;k1<dim;++k1) 
            rs_tensor[k2][k1] = 0.5*(rs_tensor_temp[k2][k1]+rs_tensor_temp[k1][k2]); //it has to be symmetric
          
        for (int k2=0;k2<dim;++k2)
          rs_tensor_trace += rs_tensor[k2][k2];

        for (int k2=0;k2<dim;++k2) 
          rs_tensor[k2][k2] -= rs_tensor_trace/dim; //make it traceless

        // assign random tensor for each particle
        for (int k=0;k<dim;++k) {
          rstress_x[i][k] = rs_tensor[k][0];
          rstress_y[i][k] = rs_tensor[k][1];
        }
        if (dim > 2)
          for (int k=0;k<dim;++k)
            rstress_z[i][k] = rs_tensor[k][2];
      }
    }
    // communicate for ghost particles; we only setup for nlocal (particles owned by this processor)
    temp3 = rstress_x;
    comm_variable = TempVector;
    comm_forward = 3;
    comm->forward_comm_pair(this);

    temp3 = rstress_y;
    comm_variable = TempVector;
    comm_forward = 3;
    comm->forward_comm_pair(this);

    if (dim > 2) {
      temp3 = rstress_z;
      comm_variable = TempVector;
      comm_forward = 3;
      comm->forward_comm_pair(this);
    }
    temp3 = NULL; // to prevent potential misuse of this temporary variable
  }

  FUNCT_EXIT(comm->me);
}

#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeRandomStress"
void PairISPH::computeForceFromRandomStress() {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");  
  FUNCT_EXIT(comm->me);
}

// END::RS - Random stress

// BEGIN::RT - Transport equation
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeSoluteTransport"
void PairISPH::computeSoluteTransport() {
  FUNCT_ENTER(comm->me);

  if (update->ntimestep == 0) {
    FUNCT_EXIT(comm->me);
    return;
  }

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_solute_transport);

  for (unsigned int cid=0;cid<ISPH_MAX_CONCENTRATION;++cid) {
    if (tr.mask[cid]) {
      const int nlocal = atom->nlocal;

      li_solver->setNodalMap(nodalmap);
      li_solver->setMatrix(A.crs);
      prec->setMatrix(A.crs);

      // Solute transport equation
      {
        // wrap/create [nlocal x 1] multivector for x and b
        li_solver->createSolutionMultiVector(atom->concentration[cid], nlocal, 1);
        li_solver->createLoadMultiVector(NULL, nlocal, 1);

        // b = x (=concentration);
        auto x = li_solver->getSolutionMultiVector();
        auto b = li_solver->getLoadMultiVector();
        *b = *x;

        // build system of equations
        computeSoluteTransportSpecies(tr.dcoeff[cid], b->Values());
      }

      // solve solute transport
      {
        Teuchos::TimeMonitor evalTimerMonitor(*g_timer_solve_solute_transport_species);

        string name = "Solute Transport: Species " + to_string(cid);
        li_solver->solveProblem(prec, name.c_str());
      }

      temp = atom->concentration[cid];
      comm_variable = TempScalar;
      comm_forward = 1;
      comm->forward_comm_pair(this);
      temp = NULL;

      // invalidate A
      A.is_filled = 0;
    }
  }

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeSoluteTransport"
void PairISPH::computeSoluteTransportSpecies(const double dcoeff, double *b) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}
// END::RT - Transport equation


// BEGIN::NS - Pressure correction scheme
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeHelmholtz"
void PairISPH::computeHelmholtz(double *b) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeBlockHelmholtz"
void PairISPH::computeBlockHelmholtz(double *b) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computePoisson"
void PairISPH::computePoisson(double *b) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::correctVelocity"
void PairISPH::correctVelocity() {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::correctPressure"
void PairISPH::correctPressure() {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PairISPH::computeIncompressibleNavieStokes"
void PairISPH::computeIncompressibleNavierStokes() {
  FUNCT_ENTER(comm->me);

  if (update->ntimestep == 0) {
    FUNCT_EXIT(comm->me);
    return;
  }

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_imcompressible_navier_stokes);

  const int nlocal = atom->nlocal;
  const int dim = domain->dimension;
  const int *type = atom->type;

  li_solver->setNodalMap(nodalmap);
  li_solver->setMatrix(A.crs);
  prec->setMatrix(A.crs);

  // Helmholtz:
  {
    // copy v -> vstar with tranpose; used as initial guess
    // all atom associate varialbes has dimension 3 for vectors
    util.copyDenseMatrix(3, nlocal, &atom->v[0][0],
                         true, &vstar[0][0]);

    // wrap/create [nlocal x dim] multivector for x and b
    li_solver->createSolutionMultiVector(&vstar[0][0], nlocal, dim);
    li_solver->createLoadMultiVector(NULL, nlocal, dim);

    // b = x (=v);
    auto x = li_solver->getSolutionMultiVector();
    auto b = li_solver->getLoadMultiVector();
    *b = *x;

    // build blocks
    if (ns.is_block_helmholtz_enabled) {
      computeBlockHelmholtz(b->Values());
      li_solver->setBlockBegin();
      for (int k2=0;k2<dim;++k2)
        for (int k1=0;k1<dim;++k1)
        {
          // for now complete fill for all blocks
          li_solver->setBlock(k1, k2, A_blk(k1,k2));
        }
      li_solver->setBlockEnd();
      li_solver->setMatrixIsBlocked(true);
    } else {
      computeHelmholtz(b->Values());
    }

    // solve helmholtz
    {
      Teuchos::TimeMonitor evalTimerMonitor(*g_timer_solve_helmholtz);

      const double theta = (ns.theta > 0 ? ns.theta : -ns.theta);
      if (theta < ISPH_EPSILON) 
        *x = *b;        
      else 
        if (ns.is_block_helmholtz_enabled) 
          li_solver->solveBlockProblem(prec, "Block 3x3 Helmholtz");
        else 
          li_solver->solveProblem(prec, "Helmholtz");
    }
    
    // recover vstar by transposing vstar[nlocal x dim] -> vstar[dim x nlocal]
    util.transposeDenseMatrix(nlocal, 3, &vstar[0][0]);

    comm_variable = Vstar;
    comm_forward = 3;
    comm->forward_comm_pair(this);

    // invalidate A
    A.is_filled = 0;
  }

  // Poisson
  {
    // wrap/create [nlocal x 1] vector for x and b
    li_solver->createSolutionMultiVector(dp, nlocal, 1);
    li_solver->createLoadMultiVector(NULL, nlocal, 1);

    // build system of equations
    auto b = li_solver->getLoadMultiVector();
    computePoisson(b->Values());

    // consider the pure neumann case; i.e., no colloids
    Epetra_IntSerialDenseVector null_mask(nlocal);
    if (ns.singular_poisson == NullSpace) {
      FilterBinary filter;
      filter.setPairYes(Solid);
      for (int i=0;i<nlocal;++i)
        null_mask[i] = !filter.yes(type[i]);

      li_solver->setNullVectorMask(&null_mask);
      li_solver->setMatrixIsSingular(true);
    }

    // solve poisson
    {
      Teuchos::TimeMonitor evalTimerMonitor(*g_timer_solve_poisson);
      li_solver->setInitialSolution(SolverLinear::Zero);
      li_solver->solveProblem(prec, "Poisson");
    }

    // remove singular condition
    li_solver->setMatrixIsSingular(false);

    comm_variable = DeltaP;
    comm_forward = 1;
    comm->forward_comm_pair(this);

    // zero mean pressure
    if (ns.is_incremental_pressure_used)
      computeZeroMeanPressure(dp);

    // invalidate A
    A.is_filled = 0;
  }

  // correction
  correctVelocity();
  correctPressure();

  FUNCT_EXIT(comm->me);
}
// END::NS - Pressure correction scheme

// BEGIN::NS - Velocity correction scheme
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::predictAleVelocity"
void PairISPH::predictAleVelocity() {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeAlePoisson"
void PairISPH::computeAlePoisson(double *b) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::correctAleVelocity"
void PairISPH::correctAleVelocity() {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeAleHelmholtz"
void PairISPH::computeAleHelmholtz(double *b) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PairISPH::computeAleIncompressibleNavieStokes"
void PairISPH::computeAleIncompressibleNavierStokes() {
  FUNCT_ENTER(comm->me);

  if (update->ntimestep == 0) {
    FUNCT_EXIT(comm->me);
    return;
  }

  Teuchos::TimeMonitor evalTimerMonitor(*g_timer_compute_imcompressible_navier_stokes);

  const int nlocal = atom->nlocal;
  const int dim = domain->dimension;

  // Step 1: predict velocity explicitly
  {
    predictAleVelocity(); // predicted velocity is stored in vstar

    comm_variable = Vstar;
    comm_forward = 3;
    comm->forward_comm_pair(this);
  }

  // matrix setup for solving implicitly
  li_solver->setNodalMap(nodalmap);
  li_solver->setMatrix(A.crs);
  prec->setMatrix(A.crs);

  // Step 2: solve Poisson problem
  {
    // wrap/create [nlocal x 1] vector for x and b
    li_solver->createSolutionMultiVector(&atom->pressure[0], nlocal, 1);
    li_solver->createLoadMultiVector(NULL, nlocal, 1);

    // build system of equations
    auto b = li_solver->getLoadMultiVector();
    computeAlePoisson(b->Values());

    // consider the pure neumann case; i.e., no colloids
    li_solver->setMatrixIsSingular(ns.singular_poisson == NullSpace);

    // solve poisson
    {
      Teuchos::TimeMonitor evalTimerMonitor(*g_timer_solve_poisson);
      li_solver->setInitialSolution(SolverLinear::Zero);
      li_solver->solveProblem(prec, "AlePoisson");
    }

    // remove singular condition
    li_solver->setMatrixIsSingular(false);

    comm_variable = Pressure;
    comm_forward = 1;
    comm->forward_comm_pair(this);

    // zero mean pressure
    computeZeroMeanPressure(atom->pressure);

    // invalidate A
    A.is_filled = 0;
  }

  // Step 3: correct velocity; vstar will have the updated velocity field
  {
    correctAleVelocity();
  }

  // Step 4: solve Helmholtz problem
  {
    // so far, inputs are vstar = \hat\hat u^{n+1}; atom->v := u_h^{n+1}

    // build system of equatios
    li_solver->createLoadMultiVector(NULL, nlocal, dim);
    auto b = li_solver->getLoadMultiVector();
    computeAleHelmholtz(b->Values());

    // set initial guess as v; transpose v
    util.transposeDenseMatrix(3, nlocal, &atom->v[0][0]);
    li_solver->createSolutionMultiVector(&atom->v[0][0], nlocal, dim);

    // solve helmholtz
    {
      Teuchos::TimeMonitor evalTimerMonitor(*g_timer_solve_helmholtz);
      li_solver->solveProblem(prec, "AleHelmholtz");
    }

    // recover atom->v
    util.transposeDenseMatrix(nlocal, 3, &atom->v[0][0]);

    comm_variable = Velocity;
    comm_forward = 3;
    comm->forward_comm_pair(this);

    // invalidate A
    A.is_filled = 0;
  }

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PairISPH::computeUnitTests"
void PairISPH::computeUnitTests() {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::advanceTime"
void PairISPH::advanceTime() {
  FUNCT_ENTER(comm->me);

  if (!use_exact_solution)
    error->all(FLERR, "PairISPH:: Not allowed in PairISPH unless use_exact_solution flag on");

  double dt = update->dt;
  double **v = atom->v;
  double **vnp1 = vstar;
  double *p = atom->pressure;
  double **xtmp = work3;
  {
    FunctorAdvanceTimeExactSolutionBegin functor(this, dt, v, vnp1, p, xtmp);
    PairFor(functor, functor.getNumberOfWork());
  }
  {
    FunctorAdvanceTimeExactSolutionEnd functor(this, dt, v, vnp1, p);
    PairFor(functor, functor.getNumberOfWork());
  }

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::shiftParticles"
void PairISPH::shiftParticles(const double shift,
                              const double shiftcut,
                              const double nonfluidweight) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeVelocityCurl"
void PairISPH::computeVelocityCurl(double **v, double **curl) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeVelocityDivergence"
void PairISPH::computeVelocityDivergence(double **v, double *div) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "PairISPH::computeTractionVector"
void PairISPH::computeTractionVector(double *p, double **v, double **traction) {
  FUNCT_ENTER(comm->me);
  error->all(FLERR, "PairISPH:: Not allowed in PairISPH");
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PairISPH::compute"
void PairISPH::compute(int eflag, int vflag) {
  FUNCT_ENTER(comm->me);

  grow(); // sync pair data structure to atom
  int nlocal = atom->nlocal;
  int dim = domain->dimension;

  // pre-computation and normal
  computePre();

  if (boundary_particle != NoParticle)
    computeNormals();

  if (use_exact_solution) {
    // do nothing
  } else {
    // create a nodal map
    Epetra_IntSerialDenseVector idx(Copy, atom->tag, nlocal);
    nodalmap = new Epetra_Map(-1, nlocal, idx.Values(), 1, Epetra_MpiComm(world));

    // graph is constructed: tags_in_cut, tags_in_boundary
    computeGraph();

    // create A
    {
      A.crs = new Epetra_CrsMatrix(Copy, *tags_in_cut);
      A.crs->FillComplete();
      A.diagonal = new Epetra_Vector(*nodalmap);
      A.scaled_laplace_diagonal = new Epetra_Vector(*nodalmap);
      A.is_filled = 0;

      // create block matrix
      if (ns.is_block_helmholtz_enabled) {
        for (int k2=0;k2<dim;++k2) {
          A_blk.diagonal[k2] = new Epetra_Vector(*nodalmap);
          A_blk.scaled_laplace_diagonal[k2] = new Epetra_Vector(*nodalmap);
          for (int k1=0;k1<dim;++k1) {
            if (k1==k2)
              A_blk.set(k1,k2, new Epetra_CrsMatrix(Copy, *tags_in_cut));
            else
             A_blk.set(k1,k2, new Epetra_CrsMatrix(Copy, *tags_in_boundary));
            A_blk(k1,k2)->FillComplete();
          }
        }
        li_solver->createBlockMatrix(dim, "Block 3x3");
      }
    }

    // if ML is used, send coordinates
    PrecondWrapper_ML *precond_ml = dynamic_cast<PrecondWrapper_ML*>(prec);
    Epetra_MultiVector *xt = NULL;

    if (precond_ml != NULL) {
      xt = new Epetra_MultiVector(*nodalmap, 3, false);
      double **xt_ptr = xt->Pointers();

      util.copyDenseMatrix(3, nlocal, &atom->x[0][0],
                           true, &xt_ptr[0][0]);
      precond_ml->setCoordinates(dim,
                                 xt_ptr[0],
                                 xt_ptr[1],
                                 xt_ptr[2]);
    }

    // body; this is not quite necessary but let's do ths to make sure f is zeros
    // force_clear in verlet.cpp does the same job
    memset(&atom->f[0][0], 0, sizeof(double)*nlocal*3);

    if (ae.is_enabled)
      computeAppliedElectricField();

    if (pb.is_enabled) {
      computePoissonBoltzmann();
      computeElectrostaticForce();
    }

    if (st.is_enabled)
      computeSurfaceTension();

    if (ns.is_enabled) {
      if (ns.is_bond_interaction_enabled) {
        // in verlet, it computes sequentially pair and bond contribution
        // since we compute everything in pair, bond should be computed here 
        // and later should be disabled. we do this by testing the conversion. 
        auto bond = dynamic_cast<BondISPH*>(force->bond);
        if (bond != NULL) {
          bond->setComputeEnabled(true);
          force->bond->compute(eflag, vflag);
          bond->setComputeEnabled(false);
        }
      }
      if (ns.is_fluctuation_enabled)
        computeForceFromRandomStress();
        
      if (ale.is_enabled)
        computeAleIncompressibleNavierStokes();
      else
        computeIncompressibleNavierStokes();
    }

    if (tr.is_enabled)
      computeSoluteTransport();

    if (ut.is_enabled)
      computeUnitTests();

    // clear
    if (xt != NULL)
      delete xt;

    {
      if (ns.is_block_helmholtz_enabled) {
        for (int k2=0;k2<dim;++k2) {
          delete A_blk.diagonal[k2];
          delete A_blk.scaled_laplace_diagonal[k2];
          for (int k1=0;k1<dim;++k1)
            if (A_blk(k1,k2) != NULL)
              delete A_blk(k1,k2);
        }
        li_solver->freeBlockMatrix();
      }

      delete A.crs;
      delete A.diagonal;
      delete A.scaled_laplace_diagonal;
      A.is_filled = 0;
    }
    delete tags_in_cut;
    if (ns.is_block_helmholtz_enabled)
      delete tags_in_boundary;

    delete nodalmap;
  }

  ++is_once;

  Teuchos::TimeMonitor::summarize(cout);

  FUNCT_EXIT(comm->me);
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PairISPH::allocate"
void PairISPH::allocate() {
  FUNCT_ENTER(comm->me);

  allocated = 1;
  const int n = atom->ntypes;

  memory->create(pinfo, 3, n + 1, "pair:pinfo");

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(h, n + 1, n + 1, "pair:h");

  memory->create(st.pf.s, n + 1, n + 1, "pair:st:pf:s");

  FUNCT_EXIT(comm->me);
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairISPH::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
               "Illegal number of setting arguments for pair_style isph");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "PairISPH::coeff"
void PairISPH::coeff(int narg, char **arg) {
  FUNCT_ENTER(comm->me);

  if (narg < 3)
    error->all(FLERR, "Incorrect args for pair_style isph coefficients");

  if (!allocated)
    allocate();

  // ** load top-level parameterlist; check the file existence
  {
    ifstream file(arg[2]);
    if (file.good())
      Teuchos::updateParametersFromXmlFile(arg[2], Teuchos::inoutArg(g_params));
    else
      error->all(FLERR, "Pair coefficient isph (global) parameter file does not exist");
  }

  // ** debugging flag
  {
    use_exact_solution = (g_params.get("Use Exact Solution","Disabled") == "Enabled");
  }

  {
    string str("Random Number Generator");
    if (g_params.isSublist(str)) {
      auto random_conf = g_params.sublist(str);
      const int seed = random_conf.get("seed", 37482);
      random = new RanMars(lmp, seed + comm->me);
    } else {
      random = new RanMars(lmp, 37482 + comm->me);
    }
  }

  const int ntypes = atom->ntypes;
  const int max_ntypes = ntypes+1;

  // ** particle information
  {
    const string str("Particle Information");
    if (g_params.isSublist(str)) {
      auto particle_info = g_params.sublist(str);
      const string t_prefix("type:");
      for (int i=0;i<max_ntypes;++i) {
        const string itype = t_prefix+to_string(i);
        const string p = particle_info.get(itype, "none");
        if (p == "none") {
          // particle type information should be set appropriately
          //stringstream ss;
          //ss << "Particle information is not found for " << itype;
          //error->all(FLERR, ss.str().c_str());
          // if nothing there, copy type and set it movable
          pinfo[0][i] = i;
          pinfo[1][i] = 0;
          pinfo[2][i] = 0;
        } else {
          // detect type attributes
          // ----------------------

          // particle kind
          {
            if      (p.find("fluid")   !=string::npos)          pinfo[0][i] = Fluid;
            else if (p.find("solid")   !=string::npos)          pinfo[0][i] = Solid;
            else if (p.find("boundary")!=string::npos)          pinfo[0][i] = Boundary;
            else if (p.find("buffer-dirichlet") !=string::npos) pinfo[0][i] = BufferDirichlet;
            else if (p.find("buffer-neumann")   !=string::npos) pinfo[0][i] = BufferNeumann;
            else {
              stringstream ss;
              ss << "Incorrect particle kind for " << itype << " : " << p << endl
                 << "Available kind attributes are fluid, solid, boundary (if mls is used), buffer-dirichlet, and buffer-neumann";
              error->all(FLERR, ss.str().c_str());
            }
          }

          // particle mobility
          {
            pinfo[1][i] = (p.find("fixed")!=string::npos);
          }

          // particle phase
          {
            const string p_prefix("phase:");
            for (int j=0;j<max_ntypes;++j) {
              const string jphase = p_prefix+to_string(j);
              if (p.find(jphase) != string::npos) {
                pinfo[2][i] = j;
                break;
              }
            }
          }

        }
      }
    } else {
      // parameterlist for particle information is necessary now
      error->all(FLERR, "Parameterlist of <Particle Information> does not exist");
      // if parameterlist is not found, copy type and set it movable
      // for (int i=0;i<max_ntypes;++i) {
      //   pinfo[0][i] = i;
      //   pinfo[1][i] = 0;
      // }
    }

    bool solid_exist = false, boundary_exist = false;
    for (int i=0;i<max_ntypes;++i) {
      if (pinfo[0][i] == Solid)
        solid_exist = true;
      if (pinfo[0][i] == Boundary)
        boundary_exist = true;
    }

    if (solid_exist && boundary_exist)
      error->all(FLERR, "Incorrect use of particle kind: solid and boundary particles should not be mixed");
    else
      boundary_particle = (solid_exist ? Solid : (boundary_exist ? Boundary : NoParticle));
  }

  {
    string str("Physics Configuration");
    bool dummy_is_enabled;
    if (g_params.isSublist(str)) {
      auto physics_conf = g_params.sublist(str);
      ae.is_enabled = (physics_conf.get("Applied Electric Field","Disabled") == "Enabled");
      pb.is_enabled = (physics_conf.get("Poisson Boltzmann","Disabled") == "Enabled");
      st.is_enabled = (physics_conf.get("Surface Tension","Disabled") == "Enabled");
      tr.is_enabled = (physics_conf.get("Solute Transport","Disabled") == "Enabled");
      ns.is_enabled = (physics_conf.get("Incompressible Navier Stokes","Disabled") == "Enabled");
      ut.is_enabled = (physics_conf.get("Unit Tests","Disabled") == "Enabled");
      dummy_is_enabled = (physics_conf.get("Dummy","Disabled") == "Enabled");
    }
    if (!ae.is_enabled && 
        !pb.is_enabled && 
        !ns.is_enabled && 
        !ut.is_enabled && 
        !tr.is_enabled && 
        !dummy_is_enabled)
      error->all(FLERR, "No physics is configured; set Physics Configuration");
  }

  // ** st setup
  if (st.is_enabled) {
    string str("Surface Tension");
    if (g_params.isSublist(str)) {
      auto st_param = g_params.sublist(str);
      auto model = st_param.get("Modeling Method","ContinuumSurfaceForce");

      if (model == "ContinuumSurfaceForce") {
        st.model = ContinuumSurfaceForce;

        if (st_param.isSublist(model)) {
          auto csf_param = st_param.sublist(model);
          auto color = csf_param.get("Color Gradient","Corrected");
          if      (color == "Corrected")
            st.csf.color = new ColorFunctionCorrected();
          else if (color == "Adami")
            st.csf.color = new ColorFunctionAdami();
          else
            error->all(FLERR, "Color Gradient is not in supported list: Corrected and Adami");

          st.csf.alpha = csf_param.get("alpha", 1.0);
          st.csf.theta = csf_param.get("theta", 0.0);
          st.csf.epsilon = csf_param.get("epsilon", 0.01);
          st.csf.kappa = csf_param.get("kappa", 100.0);
        } else {
          error->all(FLERR, "ParameterList is not found: ContinuumSurfaceForce");
        }
      } else if (model == "PairwiseForce") {
        st.model = PairwiseForce;

        if (st_param.isSublist(model)) {
          auto pf_param = st_param.sublist(model);
          auto pf_model = pf_param.get("Force Model","TartakovskyMeakin");

          if      (pf_model == "TartakovskyMeakin")
            st.pf.force = new PairwiseForceFunction_TartakovskyMeakin(domain->dimension);
          else if (pf_model == "TartakovskyPanchenkoVar1")
            st.pf.force = new PairwiseForceFunction_TartakovskyPanchenkoVar1(domain->dimension);
          else if (pf_model == "TartakovskyPanchenkoVar2")
            st.pf.force = new PairwiseForceFunction_TartakovskyPanchenkoVar2(domain->dimension);
          else
            error->all(FLERR, "Force Model is not in supported list: TartakovskyMeakin, TartakovskyPanchenkoVar1, TartakovskyPanchenkoVar2");

          const string s_prefix("s:");
          for (int j=0;j<max_ntypes;++j) {
            for (int i=0;i<max_ntypes;++i) {
              const string sij = s_prefix+to_string(i)+":"+to_string(j);
              const double val = pf_param.get(sij, 0.0);
              const int iphase = getParticlePhase(i);
              const int jphase = getParticlePhase(j);

              st.pf.s[iphase][jphase] = val;
            }
          }
        } else {
          error->all(FLERR, "ParameterList is not found: PairwiseForce");
        }
      } else {
        error->all(FLERR, "Parameter is not supported: ContinuumSurfaceForce, PairwiseForce");
      }
    } else {
      error->all(FLERR, "Parameter is not supported: Surface Tension");
    }
  }

  // ** rt setup
  if (tr.is_enabled) {
    string str("Solute Transport");
    if (g_params.isSublist(str)) {
      auto rt_param = g_params.sublist(str);

      // multiple species
      const string d_prefix("d:");
      for (int i=0;i<ISPH_MAX_CONCENTRATION;++i) {
        const string di = d_prefix+to_string(i);
        tr.mask[i] = false;
        if (rt_param.isParameter(di)) {
          tr.dcoeff[i] = rt_param.get(di, 0.0);
          tr.mask[i] = true;
        }
      }

      // time stepping
      tr.theta = rt_param.get("theta", 0.5);

      // eulerian transport
      tr.is_particle_fixed = (rt_param.get("Use Fixed Particles","Disabled") == "Enabled");

    } else {
      error->all(FLERR, "ParameterList is not found: Solute Transport");
    }
  }
  // ** ae setup
  if (ae.is_enabled) {
    string str("Applied Electric Field");
    if (g_params.isSublist(str)) {
      auto ae_param = g_params.sublist(str);
      ae.is_solid2fluid = (ae_param.get("Solid2Fluid", "Disabled") == "Enabled");
      ae.smooth_phi = ae_param.get("SmoothPhi", 0);

      ae.e[0] = ae_param.get("e.x", 0.0);
      ae.e[1] = ae_param.get("e.y", 0.0);
      ae.e[2] = ae_param.get("e.z", 0.0);

    } else {
      error->all(FLERR, "Fail to find out param sublist: Poisson Boltzmann");
    }
  } else {
    ae.is_solid2fluid = false;
    ae.smooth_phi = 0;

    ae.e[0] = 0.0;
    ae.e[1] = 0.0;
    ae.e[2] = 0.0;
  }
  
  // ** pb setup
  if (pb.is_enabled) {
    string str("Poisson Boltzmann");
    if (g_params.isSublist(str)) {
      auto pb_param = g_params.sublist(str);
      pb.psiref = pb_param.get("psiref", 1.0);
      pb.is_linearized = pb_param.get("linearized", 0);

      string extra_f("Extra F");
      pb.has_extra_f = pb_param.isSublist(extra_f);
      if (pb.has_extra_f)
        pb.extra_f = pb_param.sublist(extra_f);

      pb.ezcb = pb_param.get("ezcb", 0.0);

      pb.e[0] = pb_param.get("e.x", 0.0);
      pb.e[1] = pb_param.get("e.y", 0.0);
      pb.e[2] = pb_param.get("e.z", 0.0);

      pb.gamma = pb_param.get("gamma", 0.0);

      switch (boundary_particle) {
      case NoParticle: {
        pb.boundary = NoBoundaryCond;
        break;
      }
      case Solid: {
        auto boundary = pb_param.get("Boundary","ConstExtension");

        if      (boundary == "ConstExtension")       pb.boundary = ConstExtension;
        else if (boundary == "HomogeneousNeumann")   pb.boundary = HomogeneousNeumann;
        else if (boundary == "MorrisHolmes")         pb.boundary = MorrisHolmes;
        else if (boundary == "MorrisNormal")         pb.boundary = MorrisNormal;
        else
          error->all(FLERR, "Parameter is not supported: ConstExtension, HomogeneousNeumann, MorrisHolmes and MorrisNormal");
        break;
      }
      case Boundary: {
        pb.boundary = ConstExtension;
        break;
      }
      }

    } else {
      error->all(FLERR, "Fail to find out param sublist: Poisson Boltzmann");
    }
  } else {
    pb.psiref  = 1.0;
    pb.is_linearized = 0;
    pb.has_extra_f = false;

    pb.ezcb = 0.0;
    pb.e[0] = 0.0;
    pb.e[1] = 0.0;
    pb.e[2] = 0.0;

    pb.gamma = 0.0;

    if (boundary_particle)
      pb.boundary = ConstExtension; // default
    else
      pb.boundary = NoBoundaryCond;
  }

  // ** ns setup
  if (ns.is_enabled) {
    string str("Incompressible Navier Stokes");
    if (g_params.isSublist(str)) {
      auto ns_param = g_params.sublist(str);

      // thyra block matrix solver interface
      ns.is_block_helmholtz_enabled = (ns_param.get("Block Helmholtz","Disabled") == "Enabled");

      // theta time stepping scheme
      ns.theta = ns_param.get("theta", 0.5);

      // Wenxiao robin boundary condition coefficient
      ns.beta  = ns_param.get("beta", 1.0);

      // constant force field
      ns.g[0]  = ns_param.get("g.x", 0.0);
      ns.g[1]  = ns_param.get("g.y", 0.0);
      ns.g[2]  = ns_param.get("g.z", 0.0);

      // Wenxiao random stress for fluctuating hydrodynamics
      ns.is_fluctuation_enabled = (ns_param.get("Apply Thermal Fluctuation","Disabled") == "Enabled");
      ns.kBT = ns_param.get("kBT", 0.0);

      // molecule bond interaction setup
      ns.is_bond_interaction_enabled = (ns_param.get("Apply Bond Interaction","Disabled") == "Enabled");

      // ALE flags
      ale.is_enabled = (ns_param.get("Arbitrary Lagrangian Eulerian","Disabled") == "Enabled");
      if (ale.is_enabled)
        ale.tdiff = new TimeBDF(ns_param.get("order", 2));

      // Incremental pressue
      ns.is_incremental_pressure_used = (ns_param.get("Use Incremental Pressure","Enabled") == "Enabled");

      // Momentum preserve form
      ns.is_momentum_preserve_operator_used = (ns_param.get("Use Momentum Preserve Operator","Enabled") == "Enabled");

      // Singular problem technique
      auto singular_poisson = ns_param.get("Singular Poisson","NotSingular");

      if      (singular_poisson.find("NotSingular") != string::npos) ns.singular_poisson = NotSingularPoisson;
      else if (singular_poisson.find("NullSpace")   != string::npos) ns.singular_poisson = NullSpace;
      else if (singular_poisson.find("PinZero")     != string::npos) ns.singular_poisson = PinZero;
      else if (singular_poisson.find("DoubleDiag")  != string::npos) ns.singular_poisson = DoubleDiag;
      else
        error->all(FLERR, "Parameter is not supported: NullSpace, PinZero, DoubleDiag");

      //const bool is_singular_poisson_parameter_forced = (singular_poisson.find("Forced") != string::npos);

      //if (ns.singular_poisson == NotSingularPoisson && !ns.is_incremental_pressure_used)
      //  error->all(FLERR, "Non-incremental pressure formulat cannot be combined with NotSingularPoisson");

      // boundary particle and over-riding necessary parameters
      switch (boundary_particle) {
      case NoParticle: {
        ns.boundary = NoBoundaryCond;

        // over-ride if the parameter is setup wrong
        if (ns.singular_poisson == NotSingularPoisson)
          ns.singular_poisson = NullSpace;
        break;
      }
      case Solid: {
        auto boundary = ns_param.get("Boundary","MorrisHolmes");

        if      (boundary == "NavierSlip")         ns.boundary = NavierSlip;
        else if (boundary == "HomogeneousNeumann") ns.boundary = HomogeneousNeumann;
        else if (boundary == "Dirichlet")          ns.boundary = Dirichlet;
        else if (boundary == "MorrisHolmes")       ns.boundary = MorrisHolmes;
        else if (boundary == "MorrisNormal")       ns.boundary = MorrisNormal;
        else
          error->all(FLERR, "Parameter is not supported: NavierSlip, HomogeneousNeumann, Dirichlet, MorrisHolmes and MorrisNormal");

        // over-ride if the parameter is setup wrong
        //if (!is_singular_poisson_parameter_forced)
        //  ns.singular_poisson = NotSingularPoisson;
        break;
      }
      case Boundary: {
        ns.boundary = ConstExtension;

        if (!ale.is_enabled)
          ns.singular_poisson = NotSingularPoisson;
        break;
      }
      }

    } else {
      error->all(FLERR, "Fail to find out param sublist: Incompressible Navier Stokes");
    }
  } else {

    ns.is_block_helmholtz_enabled = false;

    ns.theta = 0.5;
    ns.beta  = 1000.0;
    ns.g[0]  = 0.0;
    ns.g[1]  = 0.0;
    ns.g[2]  = 0.0;
    ns.is_fluctuation_enabled = false;
    ns.kBT   = 0.0;

    switch (boundary_particle) {
    case NoParticle: // SPH or MLS -> set default SPH
      ns.boundary = NoBoundaryCond;
      ns.singular_poisson = NullSpace;
      ale.is_enabled = false;
      ale.tdiff = NULL;
      break;
    case Solid:      // SPH
      ns.boundary = MorrisHolmes;
      ns.singular_poisson = NotSingularPoisson;
      ale.is_enabled = false;
      ale.tdiff = NULL;
      break;
    case Boundary:   // MLS
      ns.boundary = ConstExtension;
      ns.singular_poisson = NullSpace;
      ale.is_enabled = true;
      ale.tdiff = new TimeBDF(2);
      break;
    }

  }

  // over-ride beta to be zero if ntypes = 1 (homogeneous material)
  if (atom->ntypes == 1)
    ns.beta = 0.0;

  // show the parameterlist
  if (comm->me == 0)
    cout << g_params << endl;

  // create solvers and customize them with given xml params
  initializeSolvers();

  FUNCT_EXIT(comm->me);
}

/* ----------------------------------------------------------------------
   init pair style
   ------------------------------------------------------------------------- */

void PairISPH::init_style() {
  // sanity check

  // request a full neighbor list; not enough documentation about this
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

double PairISPH::init_one(int i, int j) {

  if (setflag[i][j] == 0)
    error->all(FLERR,"All pair isph coeffs are not set");

  cutsq[j][i] = cutsq[i][j];
  h[j][i] = h[i][j];

  return sqrt(cutsq[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairISPH::single(int i, int j, int itype, int jtype,
                        double rsq, double factor_coul, double factor_lj,
                        double &fforce) {
  // not sure if this needs to be zero or not.
  //fforce = 0.0;
  return 0.0;
}

/* ----------------------------------------------------------------------
   communication
   ------------------------------------------------------------------------- */
int PairISPH::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc) {
  int m = 0;
  switch (comm_variable) {
  case Vfrac: {
    double *vfrac = atom->vfrac;
    for (int ii=0;ii<n;++ii) {
      int i = list[ii];
      buf[m++] = vfrac[i];
    }
    break;
  }
  case NormalVector: {
    for (int ii=0;ii<n;++ii) {
      int i = list[ii];
      buf[m++] = normal[i][0];
      buf[m++] = normal[i][1];
      buf[m++] = normal[i][2];
    }
    break;
  }
  case Velocity: {
    double **v = atom->v;
    for (int ii=0;ii<n;++ii) {
      int i = list[ii];
      buf[m++] = v[i][0];
      buf[m++] = v[i][1];
      buf[m++] = v[i][2];
    }
    break;
  }
  case Pressure: {
    double *p = atom->pressure;
    for (int ii=0;ii<n;++ii) {
      int i = list[ii];
      buf[m++] = p[i];
    }
    break;
  }
  case Vstar: {
    for (int ii=0;ii<n;++ii) {
      int i = list[ii];
      buf[m++] = vstar[i][0];
      buf[m++] = vstar[i][1];
      buf[m++] = vstar[i][2];
    }
    break;
  }
  case DeltaP: {
    for (int ii=0;ii<n;++ii) {
      int i = list[ii];
      buf[m++] = dp[i];
    }
    break;
  }
  case Psi: {
    double *psi = atom->psi;
    for (int ii=0;ii<n;++ii) {
      int i = list[ii];
      buf[m++] = psi[i];
    }
    break;
  }
  case WorkScalar: {
    for (int ii=0;ii<n;++ii) {
      int i = list[ii];
      buf[m++] = work[i];
    }
    break;
  }
  case WorkVector: {
    for (int ii=0;ii<n;++ii) {
      int i = list[ii];
      buf[m++] = work3[i][0];
      buf[m++] = work3[i][1];
      buf[m++] = work3[i][2];
    }
    break;
  }
  case TempScalar: {
    for (int ii=0;ii<n;++ii) {
      int i = list[ii];
      buf[m++] = temp[i];
    }
    break;
  }
  case TempVector: {
    for (int ii=0;ii<n;++ii) {
      int i = list[ii];
      buf[m++] = temp3[i][0];
      buf[m++] = temp3[i][1];
      buf[m++] = temp3[i][2];
    }
    break;
  }
  default:
    error->all(FLERR, "Comm: packing empty buffer");
    break;
  }
  return m;
}

void PairISPH::unpack_forward_comm(int n, int first, double *buf) {
  switch (comm_variable) {
  case Vfrac: {
    memcpy(&atom->vfrac[first], &buf[0], sizeof(double)*n);
    break;
  }
  case NormalVector: {
    memcpy(&normal[first][0], &buf[0], sizeof(double)*n*3);
    break;
  }
  case Velocity: {
    memcpy(&atom->v[first][0], &buf[0], sizeof(double)*n*3);
    break;
  }
  case Pressure: {
    memcpy(&atom->pressure[first], &buf[0], sizeof(double)*n);
    break;
  }
  case Vstar: {
    memcpy(&vstar[first][0], &buf[0], sizeof(double)*n*3);
    break;
  }
  case DeltaP: {
    memcpy(&dp[first], &buf[0], sizeof(double)*n);
    break;
  }
  case Psi: {
    memcpy(&atom->psi[first], &buf[0], sizeof(double)*n);
    break;
  }
  case WorkScalar: {
    memcpy(&work[first], &buf[0], sizeof(double)*n);
    break;
  }
  case WorkVector: {
    memcpy(&work3[first][0], &buf[0], sizeof(double)*n*3);
    break;
  }
  case TempScalar: {
    memcpy(&temp[first], &buf[0], sizeof(double)*n);
    break;
  }
  case TempVector: {
    memcpy(&temp3[first][0], &buf[0], sizeof(double)*n*3);
    break;
  }
  default:
    error->all(FLERR, "Comm: unpacking empty buffer");
  }
}
