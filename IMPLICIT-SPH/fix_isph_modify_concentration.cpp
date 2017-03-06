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

#include "fix_isph_modify_concentration.h"

using namespace LAMMPS_NS;
using namespace FixConst;
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyConcentration::FixISPH_ModifyConcentration"
FixISPH_ModifyConcentration::FixISPH_ModifyConcentration(LAMMPS *lmp, int narg, char **arg) 
  : Fix(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Fix ipsh/modify/concentration command requires a matching atom_style, i.e., isph");

  if (narg < 11)
    error->all(FLERR,"Illegal number of arguments for fix isph/modify/concentration command");
  
  // default values
  _stat = NULL;
  _concentration = NULL;
  _value = 0.0;

  _xmin = 0.0; _xmax = 0.0;
  _ymin = 0.0; _ymax = 0.0;
  _zmin = 0.0; _zmax = 0.0;
  
  for (int i=0;i<modify->ncompute;++i) 
    if (strcmp(modify->compute[i]->style,"isph/status") == 0)
      _stat = modify->compute[i];
  
  if (strcmp(arg[3], "null") != 0) {
    const int iconcentration = modify->find_compute(arg[3]);
    if (iconcentration >= 0) {
      _concentration = modify->compute[iconcentration];
      if (strcmp(_concentration->style,"isph/concentration") != 0)
        error->all(FLERR,"Input argument for compute style does not match to isph/concentration");
    }
  }
  
  _t_begin = force->numeric(FLERR, arg[4]);
  _t_end   = force->numeric(FLERR, arg[5]);

  _value = force->numeric(FLERR, arg[6]);
 
  _xmin  = force->numeric(FLERR, arg[7]);
  _xmax  = force->numeric(FLERR, arg[8]);

  _ymin  = force->numeric(FLERR, arg[9]);
  _ymax  = force->numeric(FLERR, arg[10]);
 
  _zmin  = force->numeric(FLERR, arg[11]);
  _zmax  = force->numeric(FLERR, arg[12]);

  _mask = ~0;
  if (narg >= 13) {
    _mask = 0;
    for (int i=13;i<narg;++i)
      _mask |= (1<<(atoi(arg[i])-1));
  }
  
  if (comm->me == 0) {
    const int w = 15; 
    cout << "FixISPH_ModifyConcentration:: compute field " << arg[3] << ", " << _value << endl;
    cout << "                              mask " << _mask << endl; 
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
#define __FUNCT__ "FixISPH_ModifyConcentration::~FixISPH_ModifyConcentration"
FixISPH_ModifyConcentration::~FixISPH_ModifyConcentration() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyConcentration::setmask"
int FixISPH_ModifyConcentration::setmask() {
  FUNCT_ENTER(comm->me);

  int mask = 0;

  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;

  FUNCT_EXIT(comm->me);

  return mask;
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyConcentration::initial_integrate"
void FixISPH_ModifyConcentration::initial_integrate(int dummy) {
  FUNCT_ENTER(comm->me);
  modify_concentration();
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyConcentration::final_integrate"
void FixISPH_ModifyConcentration::final_integrate() {
  FUNCT_ENTER(comm->me);
  modify_concentration();
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyConcentration::modify_concentration"
void FixISPH_ModifyConcentration::modify_concentration() {
  FUNCT_ENTER(comm->me);

  // current time 
  const int ntimestep = update->ntimestep;
  const double dt = update->dt;
  const double t = (_stat == NULL ? ntimestep*dt : _stat->vector[0]);

  // modify time for a specific range
  if (t > _t_begin && (_t_end < 0 || t < _t_end)) {
    // get accesses to data structures
    const int nsize = (atom->nlocal + atom->nghost);
    
    double **x = atom->x;
    int *type = atom->type;
    
    const bool check_x = (_xmin != _xmax);
    const bool check_y = (_ymin != _ymax);
    const bool check_z = (_zmin != _zmax);
    
    double *concentration;
    if (_concentration != NULL) {
      _concentration->compute_peratom(); // assignment (cost nothing).
      concentration = _concentration->vector_atom;
    } else {
      // default 
      concentration = atom->concentration[0]; 
    }
    
    for (int i=0;i<nsize;++i) 
      if (_mask & (1<<(type[i]-1))) {
        double *x_at_i = x[i];
        bool is_modify = true;
        
        is_modify &= (check_x ? (_xmin < x_at_i[0] && x_at_i[0] < _xmax) : true);
        is_modify &= (check_y ? (_ymin < x_at_i[1] && x_at_i[1] < _ymax) : true);
        is_modify &= (check_z ? (_zmin < x_at_i[2] && x_at_i[2] < _zmax) : true);
        
        if (is_modify)
          concentration[i] = _value;
      }
    
    // PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
    // if (pair != NULL) {
    //   pair->refreshParticles();
    // } else {
    //   error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");
    // }
  }

  FUNCT_EXIT(comm->me);
}


