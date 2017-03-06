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

#include "fix_isph.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */  
#undef  __FUNCT__
#define __FUNCT__ "FixISPH::FixISPH"
FixISPH::FixISPH(LAMMPS *lmp, int narg, char **arg) 
  : Fix(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1) 
    error->all(FLERR, "Fix ipsh command requires a matching atom_style, i.e., isph");
  if (narg < 3)
    error->all(FLERR,"Illegal number of arguments for fix isph command");
  
  time_integrate = 1;

  _stat = NULL;
  if (narg > 3 && strcmp(arg[3], "null") != 0) {
    const int istat = modify->find_compute(arg[3]);
    if (istat >= 0) {
      _stat = modify->compute[istat];
      if (strcmp(_stat->style,"isph/status") != 0)
        error->all(FLERR,"Input argument for compute style does not match to isph/status");
    }
  }

  _use_variable_timestep = false;
  if (_stat != NULL && narg == 7) {
    _use_variable_timestep = true;
    _cfl  = force->numeric(FLERR, arg[4]);
    _dx   = force->numeric(FLERR, arg[5]);
    _umin = force->numeric(FLERR, arg[6]);
  }
  
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH::~FixISPH"
FixISPH::~FixISPH() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH::setmask"
int FixISPH::setmask() {
  FUNCT_ENTER(comm->me);

  int mask = 0;

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair != NULL) {
    // ALE (aka velocity correction) move/shift before pair->compute()
    mask |= INITIAL_INTEGRATE; 

    // pressure correction move/shift it pair->compute()
    mask |= FINAL_INTEGRATE;   

  } else {
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");
  }
  
  FUNCT_EXIT(comm->me);

  return mask;
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH::initial_integrate"
void FixISPH::initial_integrate(int dummy) {
  FUNCT_ENTER(comm->me);

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair != NULL) {
    // ALE (velocity correction) move particles
    if (pair->isAleEnabled()) {
      pair->advanceTime();
      pair->refreshParticles();
    } 
    show_current_timestep(__FUNCT__);
  } else {
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");
  }

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH::final_integrate"
void FixISPH::final_integrate() {
  FUNCT_ENTER(comm->me);

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair != NULL) {
    // pressure correction move particles
    if (!pair->isAleEnabled()) {
      pair->advanceTime();
    } 
    show_current_timestep(__FUNCT__);
  } else {
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");
  }

  if (_stat != NULL) {
    _stat->compute_vector();

    if (_use_variable_timestep) {
      const double *result = _stat->vector;
      const double umax  = max(result[8], _umin);
      update->dt = _cfl*_dx/umax;
    }
  }

  string str = "isph/status";
  for (int i=0;i<modify->ncompute;++i) {
    if (strncmp(modify->compute[i]->style,str.c_str(), str.size()) == 0) {
      auto stat = modify->compute[i];
      if (stat != _stat)
        stat->compute_vector();
    }
  }
  
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
void FixISPH::show_current_timestep(const char *description) {
  if (comm->me == 0) {
    int ntimestep = update->ntimestep;
    double dt = update->dt;
    
    if (_use_variable_timestep)
      cout << endl 
           << description << endl
           << "         Timestep     = " << ntimestep 
           << "         Current time = " << _stat->vector[0]
           << endl;
    else
      cout << endl 
           << description << endl
           << "         Timestep     = " << ntimestep 
           << "         Current time = " << ntimestep*dt
           << endl;
  }
}
