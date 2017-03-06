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
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

#include <iostream>
#include <fstream>

#include "utils.h"

// not used for sph
//#include "pair_isph_corrected.h"
//#include "pair_isph_mls.h"

#include "compute_isph_status_flow_past_cylinder.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_StatusFlowPastCylinder::ComputeISPH_StatusFlowPastCylinder"
ComputeISPH_StatusFlowPastCylinder::ComputeISPH_StatusFlowPastCylinder(LAMMPS *lmp, int narg, char **arg) 
  : Compute(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Compute ipsh/status command requires a matching atom_style, i.e., isph");

  if (narg < 16)
    error->all(FLERR,"Illegal number of arguments for compute isph command");

  // arg[3] : set output file 
  _os = &cout;
  if (strcmp(arg[3], "null") != 0)
    _os = new ofstream(arg[3]);

  // arg[4,5,6] cylinder radius length 
  _cyl_radius = force->numeric(FLERR, arg[4]);
  _cyl_lo     = force->numeric(FLERR, arg[5]);
  _cyl_hi     = force->numeric(FLERR, arg[6]);

  // arg[7,8,9] cylinder center
  _cyl_center[0] = force->numeric(FLERR, arg[7]);
  _cyl_center[1] = force->numeric(FLERR, arg[8]);
  _cyl_center[2] = force->numeric(FLERR, arg[9]);

  // arg[10] cylinder axis component
  if      (strcmp(arg[10],"x") == 0) _axis = 0;
  else if (strcmp(arg[10],"y") == 0) _axis = 1;
  else if (strcmp(arg[10],"z") == 0) _axis = 2;
  else
    error->all(FLERR,"Not supported keyword: x, y, z");

  // arg[11,12,13] drag vector
  _drag[0] = force->numeric(FLERR, arg[11]);
  _drag[1] = force->numeric(FLERR, arg[12]);
  _drag[2] = force->numeric(FLERR, arg[13]);

  {
    const double norm = sqrt(util.dotVectors(3, _drag, _drag));
    for (int i=0;i<3;++i)
      _drag[i] /= norm;
  }
  
  // arg[14,15,16] lift vector
  _lift[0] = force->numeric(FLERR, arg[14]);
  _lift[1] = force->numeric(FLERR, arg[15]);
  _lift[2] = force->numeric(FLERR, arg[16]);

  {
    const double norm = sqrt(util.dotVectors(3, _lift, _lift));
    for (int i=0;i<3;++i)
      _lift[i] /= norm;
  }

  // type masks
  _mask = ~0;
  if (narg > 17) {
    _mask = 0;
    for (int i=17;i<narg;++i)
      _mask |= (1<<(atoi(arg[i])-1));
  }

  vector_flag = 1;
  size_vector = ISPH_COMPUTE_STATUS_MAX_SIZE;

  for (int i=0;i<ISPH_COMPUTE_STATUS_MAX_SIZE;++i)
    _stat[i] = 0.0;

  vector = _stat;

  _nmax = 0;
  _traction = NULL;

  if (comm->me == 0) {
    cout << " - ComputeISPH_StatusFlowPastCylinder - " << endl
         << "   output data file = " << arg[3] << endl
         << "   cylinder radius = " << _cyl_radius << endl
         << "   cylinder low and high = " << _cyl_lo << ", " <<  _cyl_hi << endl
         << "   cylinder center = "
         << _cyl_center[0] << ", " << _cyl_center[1] << ", " << _cyl_center[2] << endl
         << "   cylinder axis = " << arg[15] << ", " << _axis << endl
         << endl
         << "   drag = "
         << _drag[0] << ", " << _drag[1] << ", " << _drag[2] << endl
         << "   lift = "
         << _lift[0] << ", " << _lift[1] << ", " << _lift[2] << endl
         << endl;
  }
  
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_StatusFlowPastCylinder::~ComputeISPH_StatusFlowPastCylinder"
ComputeISPH_StatusFlowPastCylinder::~ComputeISPH_StatusFlowPastCylinder() {
  FUNCT_ENTER(comm->me);

  if (_os != &cout)
    delete _os;
  
  if (_traction != NULL)
    memory->destroy(_traction);

  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_StatusFlowPastCylinder::init"
void ComputeISPH_StatusFlowPastCylinder::init() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);  
}

#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_StatusFlowPastCylinder::compute_vector"
void ComputeISPH_StatusFlowPastCylinder::compute_vector() {

  FUNCT_ENTER(comm->me);

  invoked_vector = update->ntimestep;

  const int nlocal = atom->nlocal;
  
  // grow tranction array if necessary
  if (nlocal > _nmax) {
    _nmax = atom->nmax;
    memory->grow(_traction, _nmax, 3, "compute:status-flow-pst-cylinder:traction");
    memset(&_traction[0][0], 0, sizeof(double)*_nmax*3);
  }

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair == NULL)
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");

  double stat[ISPH_COMPUTE_STATUS_MAX_SIZE] = {};

  // compute traction vector
  pair->computeTractionVector(atom->pressure, atom->v, _traction);

  int *type = atom->type; 
  double **x = atom->x;
  const int dim = domain->dimension;

  unsigned int cnt = 0;
  for (int i=0;i<nlocal;++i) {
    if (_mask & (1<<(type[i]-1))) {
      const double *t = &_traction[i][0];
      stat[1] += util.dotVectors(dim, t, _drag);  // Cd
      stat[2] += util.dotVectors(dim, t, _lift);  // Cl

      double s[3], r[3];
      for (int k=0;k<dim;++k)
        r[k] = x[i][k] - _cyl_center[k];

      s[0] = t[1]*r[2] - t[2]*r[1];
      s[1] = t[2]*r[0] - t[0]*r[2];
      s[2] = t[0]*r[1] - t[1]*r[0];

      stat[3] += s[_axis];                        // Cm
      ++cnt;
    }
  }
  stat[0] = cnt;

  memset(&_stat[1], 0, sizeof(double)*(ISPH_COMPUTE_STATUS_MAX_SIZE-1));

  MPI_Allreduce(&stat[0], &_stat[1], 5, MPI_DOUBLE, MPI_SUM, world);

  const double area = 2*M_PI*_cyl_radius*(_cyl_hi - _cyl_lo);
  const double da = area/_stat[1];
  
  for (int i=2;i<5;++i)
    _stat[i] *= da;

  const int w = 15;

  // write root processor only
  if (comm->me == 0) {
    (*_os) << scientific
           << setw(w) << _stat[0] << "   "    // time
           << setw(w) << _stat[1] << "   "    // # of particles
           << setw(w) << _stat[2] << "   "    // Cd
           << setw(w) << _stat[3] << "   "    // Cl
           << setw(w) << _stat[4] << "   "    // Cm
           << fixed << endl;
  }

  FUNCT_EXIT(comm->me);

}
