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

#include "string.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

#include <iostream>
#include <fstream>
#include "utils.h"

#include "pair_isph_corrected.h"
//#include "pair_isph_mls.h"

#include "compute_isph_status.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_Status::ComputeISPH_Status"
ComputeISPH_Status::ComputeISPH_Status(LAMMPS *lmp, int narg, char **arg) 
  : Compute(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Compute ipsh/status command requires a matching atom_style, i.e., isph");

  if (narg < 3)
    error->all(FLERR,"Illegal number of arguments for compute isph command");

  _os = &cout;
  if (narg >= 3 && strcmp(arg[3],"null") != 0) 
    _os = new ofstream(arg[3]);
  
  _vflag = false;
  if (narg >= 4 && strcmp(arg[4],"null") != 0) 
    if (strcmp(arg[4],"vstar") == 0) _vflag = true;

  _c[0] = 0.0; _c[1] = 0.0; _c[2] = 0.0;
  _r = -1.0;


  if (narg >= 8) {
    _c[0] = force->numeric(FLERR, arg[5]);
    _c[1] = force->numeric(FLERR, arg[6]);
    _c[2] = force->numeric(FLERR, arg[7]);
    _r = force->numeric(FLERR, arg[8]);
  }

  _mask = ~0;
  if (narg > 9) {
    _mask = 0;
    for (int i=9;i<narg;++i)
      _mask |= (1<<(atoi(arg[i])-1));
  }

  for (int i=0;i<ISPH_COMPUTE_STATUS_MAX_SIZE;++i)
    _stat[i] = 0.0;
  
  vector_flag = 1;
  size_vector = ISPH_COMPUTE_STATUS_MAX_SIZE;

  // register data field
  vector = _stat;

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_Status::~ComputeISPH_Status"
ComputeISPH_Status::~ComputeISPH_Status() {
  FUNCT_ENTER(comm->me);

  if (_os != &cout)
    delete _os;

  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_Status::init"
void ComputeISPH_Status::init() {
  FUNCT_ENTER(comm->me);

  int count = 0;
  for (int i=0;i<modify->ncompute;++i)
    count += (strcmp(modify->compute[i]->style,"isph/status") == 0);

  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute isph/status");

  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_Status::compute_vector"
void ComputeISPH_Status::compute_vector() {
  FUNCT_ENTER(comm->me);

  invoked_vector = update->ntimestep;

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair == NULL)
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");

  int *type = atom->type;

  int nlocal = atom->nlocal;

  double *density = atom->density;
  double *vfrac = atom->vfrac;
  double **v = (_vflag ? pair->getVstar() : atom->v); 
  double **x = atom->x;

  double stat[ISPH_COMPUTE_STATUS_MAX_SIZE] = {};

  // note: lammps masking groupbit is not properly set
  //       if types are modified after they are established.
  //       i.e., inline generation and porous modification 
  
  unsigned int cnt = 0;
  for (int i=0;i<nlocal;++i) {
    if (_mask & (1<<(type[i]-1))) {
      const double vx = v[i][0], vy = v[i][1], vz = v[i][2];
      const double vmagsq = vx*vx + vy*vy + vz*vz;
      const double mass = density[i]*vfrac[i];

      double xx[3];
      xx[0] = x[i][0] - _c[0];
      xx[1] = x[i][1] - _c[1];
      xx[2] = x[i][2] - _c[2];

      const double r = sqrt(xx[0]*xx[0] + xx[1]*xx[1] + xx[2]*xx[2]);

      stat[1] += vx;   // sum vx
      stat[2] += vy;   // sum vy
      stat[3] += vz;   // sum vz

      stat[4] += vfrac[i];  // volume
      stat[5] += mass;      // mass

      stat[6] += 0.5*mass*vmagsq;  // kinetic energy

      // a particular interest 
      if (_r > 0 && r <= _r)
        stat[7]  = max(stat[7], sqrt(vmagsq));  // max velocity

      ++cnt;
    }
  }
  stat[0] = cnt;;

  memset(&_stat[1], 0, sizeof(double)*(ISPH_COMPUTE_STATUS_MAX_SIZE-1));

  MPI_Allreduce(&stat[0], &_stat[1], 7, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&stat[7], &_stat[8], 1, MPI_DOUBLE, MPI_MAX, world);

  // accumulate time
  if (update->ntimestep == 0)
    _stat[0] = 0.0;
  else
    _stat[0] += update->dt;
  
  const int w = 15;

  // write root processor only
  if (comm->me == 0) {
    (*_os) << scientific
           << setw(w) << _stat[0] << "   "    // time
           << setw(w) << _stat[1] << "   "    // # of particles
           << setw(w) << _stat[2] << "   "    // vx
           << setw(w) << _stat[3] << "   "    // vy
           << setw(w) << _stat[4] << "   "    // vz
           << setw(w) << _stat[5] << "   "    // volume
           << setw(w) << _stat[6] << "   "    // mass
           << setw(w) << _stat[7] << "   "    // kinetic energy
           << setw(w) << _stat[8] << "   "    // max velocity
           << fixed << endl;
  }

  FUNCT_EXIT(comm->me);  
}

