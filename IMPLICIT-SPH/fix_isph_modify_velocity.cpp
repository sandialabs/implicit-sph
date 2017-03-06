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

#include "fix_isph_modify_velocity.h"

using namespace LAMMPS_NS;
using namespace FixConst;
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyVelocity::FixISPH_ModifyVelocity"
FixISPH_ModifyVelocity::FixISPH_ModifyVelocity(LAMMPS *lmp, int narg, char **arg) 
  : Fix(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Fix ipsh/modify/velocity command requires a matching atom_style, i.e., isph");

  if (narg < 11)
    error->all(FLERR,"Illegal number of arguments for fix isph/modify/velocity command");
  
  // default values
  _vx = 0.0; _vy = 0.0; _vz = 0.0;

  _xmin = 0.0; _xmax = 0.0;
  _ymin = 0.0; _ymax = 0.0;
  _zmin = 0.0; _zmax = 0.0;
  
  _vx = force->numeric(FLERR, arg[3]);
  _vy = force->numeric(FLERR, arg[4]);
  _vz = force->numeric(FLERR, arg[5]);
 
  _xmin  = force->numeric(FLERR, arg[6]);
  _xmax  = force->numeric(FLERR, arg[7]);

  _ymin  = force->numeric(FLERR, arg[8]);
  _ymax  = force->numeric(FLERR, arg[9]);
 
  _zmin  = force->numeric(FLERR, arg[10]);
  _zmax  = force->numeric(FLERR, arg[11]);

  _mask = ~0;
  if (narg >= 12) {
    _mask = 0;
    for (int i=12;i<narg;++i) 
      _mask |= (1<<(atoi(arg[i])-1));
  }
  
  if (comm->me == 0) {
    const int w = 15; 
    cout << "FixISPH_ModifyVelocity:: velocity " << _vx << ", " << _vy << ", " << _vz << endl;
    cout << "                         mask " << _mask << endl; 
    cout << scientific
         << setw(w) << _xmin << "   " << setw(w) << _xmax << "   "
         << setw(w) << _ymin << "   " << setw(w) << _ymax << "   "
         << setw(w) << _zmin << "   " << setw(w) << _zmax << "   "
         << fixed << endl;
  }

  time_integrate = 0;
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyVelocity::~FixISPH_ModifyVelocity"
FixISPH_ModifyVelocity::~FixISPH_ModifyVelocity() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyVelocity::setmask"
int FixISPH_ModifyVelocity::setmask() {
  FUNCT_ENTER(comm->me);

  int mask = 0;

  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;

  FUNCT_EXIT(comm->me);

  return mask;
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyVelocity::initial_integrate"
void FixISPH_ModifyVelocity::initial_integrate(int dummy) {
  FUNCT_ENTER(comm->me);
  modify_velocity();
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyVelocity::final_integrate"
void FixISPH_ModifyVelocity::final_integrate() {
  FUNCT_ENTER(comm->me);
  modify_velocity();
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyVelocity::modify_velocity"
void FixISPH_ModifyVelocity::modify_velocity() {
  FUNCT_ENTER(comm->me);

  // get accesses to data structures
  const int nsize = (atom->nlocal + atom->nghost);
  
  double **x = atom->x;
  int *type = atom->type;
  
  const bool check_x = (_xmin != _xmax);
  const bool check_y = (_ymin != _ymax);
  const bool check_z = (_zmin != _zmax);

  double **v = atom->v;
  
  for (int i=0;i<nsize;++i) 
    if (_mask & (1<<(type[i]-1))) {
      double *x_at_i = x[i];
      bool is_modify = true;
      
      is_modify &= (check_x ? (_xmin < x_at_i[0] && x_at_i[0] < _xmax) : true);
      is_modify &= (check_y ? (_ymin < x_at_i[1] && x_at_i[1] < _ymax) : true);
      is_modify &= (check_z ? (_zmin < x_at_i[2] && x_at_i[2] < _zmax) : true);
      
      if (is_modify) {
        v[i][0] = _vx;
        v[i][1] = _vy;
        v[i][2] = _vz;
      }
    }
  
  FUNCT_EXIT(comm->me);
}


