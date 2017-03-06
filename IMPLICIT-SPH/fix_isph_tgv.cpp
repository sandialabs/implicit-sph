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

#include "stdio.h"
#include "string.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "modify.h"
#include "compute.h"
#include "domain.h"
#include "math.h"

#include <iostream>
#include "utils.h"

#include "pair_isph_corrected.h"
//#include "pair_isph_mls.h"

#include "fix_isph_tgv.h"

using namespace LAMMPS_NS;
using namespace FixConst;

template<class ArgPair> 
inline void 
FixISPH_TGV::compute_error(ArgPair *pair) {
  
  int nlocal = atom->nlocal;
  int ntotal = atom->natoms;
  
  double *viscosity = atom->viscosity;
  double *density = atom->density;
  double **x = atom->x;

  double **vnp1;
  if (pair->isAleEnabled())
    vnp1 = atom->v;
  else
    vnp1 = pair->getVstar();

  double *pressure = atom->pressure;
  double dt = update->dt;
  int    ntimestep = update->ntimestep;
  double t = dt * ntimestep;
  if (_stat != NULL)
    t = _stat->vector[0];

  double errorp = 0, erroru = 0, normp = 0, normu=0, p;
  double Umax = 0.1;
  double vex[3] = { 0.0, 0.0, 0.0 };
  
  // for calculating error
  int myid = comm->me;
  double *sendbuff = new double[4], *recvbuff = new double[4];
  double pres_average = 0.0, exact_pres_average = 0.0;
  
  double partial_sum = 0, total_sum;
  for (int i = 0; i < atom->nlocal; ++i)
    partial_sum += pressure[i];
  
  MPI_Allreduce(&partial_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, world);
  pres_average = total_sum / atom->natoms;
  
  double pres_avg_diff = pres_average - exact_pres_average;
  
  for (int i = 0; i < nlocal; i++) {
    
    //error calculation (in discrete l2 norm) for Taylor-Green vortex
    p = 0.25 * density[i] * (cos(2 * x[i][0]) + cos(2 * x[i][1])) * Umax * Umax * exp(-4.0 * viscosity[i] * t);
    vex[0] = Umax * exp(-2.0 * viscosity[i] * t) * sin(x[i][0]) * cos(x[i][1]);
    vex[1] = -Umax * exp(-2.0 * viscosity[i] * t) * cos(x[i][0]) * sin(x[i][1]);
    vex[2] = 0;

    if (pair->isAleEnabled() && ntimestep < ISPH_BDF_MAX_ORDER) {
      for (int k=0;k<3;++k)
        atom->v[i][k] = vex[k];
      atom->pressure[i] = p;
    }
    
    normp += (p - exact_pres_average) * (p - exact_pres_average);
    normu += std::pow(vex[0], 2) + std::pow(vex[1], 2) + std::pow(vex[2],2);
    
    // we let pressure have average zero
    errorp += std::pow(pressure[i] - p - pres_avg_diff, 2);
    erroru += std::pow(vnp1[i][0] - vex[0], 2) + std::pow(vnp1[i][1] - vex[1], 2) + std::pow(vnp1[i][2]-vex[2],2);
  }
  
  sendbuff[0] = errorp;
  sendbuff[1] = normp;
  sendbuff[2] = erroru;
  sendbuff[3] = normu;
  
  MPI_Reduce(sendbuff, recvbuff, 4, MPI_DOUBLE, MPI_SUM, 0, world);
  if (myid == 0) {
    errorp = sqrt(recvbuff[0] / ntotal);
    normp = sqrt(recvbuff[1] / ntotal);
    erroru = sqrt(recvbuff[2] / ntotal);
    normu = sqrt(recvbuff[3] / ntotal);

    if (pair->isAleEnabled() && ntimestep < ISPH_BDF_MAX_ORDER) 
      std::cout << std::endl 
                << "FixISPH_TGV::compute_error - Store exact solution to compare time error"
                << std::endl;
    
    printf("time step:%d, time:%e \n\t pressure l2 error (norm): %.12e (%.8e)\n\t velocity l2 error (norm): %.12e (%.8e)\n", ntimestep, t, errorp, normp, erroru, normu);
  }
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_TGV::FixISPH_TGV"  
FixISPH_TGV::FixISPH_TGV(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);
  
  if (atom->isph_flag != 1)
    error->all(FLERR, "Fix ipsh command requires a matching atom_style, i.e., isph");

  if (narg < 3)
    error->all(FLERR,"Illegal number of arguments for fix isph command");

  time_integrate = 0;

  _stat = NULL;
  if (narg > 3 && strcmp(arg[3], "null") != 0) {
    int istat = modify->find_compute(arg[3]);
    if (istat >= 0) {
      _stat = modify->compute[istat];
      if (strcmp(_stat->style,"isph/status") != 0)
        error->all(FLERR,"Input argument for compute style does not match to isph/status");
    }
  }

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_TGV::~FixISPH_TGV"
FixISPH_TGV::~FixISPH_TGV() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);
}  

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_TGV::setmask"  
int FixISPH_TGV::setmask() {
  FUNCT_ENTER(comm->me);

  int mask = 0;
  mask |= FINAL_INTEGRATE;

  FUNCT_EXIT(comm->me);
  return mask;
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_TGV::final_integrate"  
void FixISPH_TGV::final_integrate() {
  FUNCT_ENTER(comm->me);

  int is_once = 0;

  if (!is_once) {
    PairISPH_Corrected *pair = dynamic_cast<PairISPH_Corrected*>(force->pair_match("isph/corrected", 0));
    if (pair != NULL) {
      compute_error<PairISPH_Corrected>(pair);
      ++is_once;
    }
  }
  if (!is_once) {
    // PairISPH_MLS *pair = dynamic_cast<PairISPH_MLS*>(force->pair_match("isph/mls", 0));
    // if (pair != NULL) {
    //   compute_error<PairISPH_MLS>(pair);
    //   ++is_once;
    // }
  }
  if (is_once != 1)
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH_Corrected, PairISPH_MLS");

  FUNCT_EXIT(comm->me);
}

