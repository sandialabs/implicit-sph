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

#include "fix_isph_quit.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */  
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Quit::FixISPH_Quit"
FixISPH_Quit::FixISPH_Quit(LAMMPS *lmp, int narg, char **arg) 
  : Fix(lmp, narg, arg) {
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

  _tfinal = -1.0;
  if (_stat != NULL && narg == 5) 
    _tfinal = force->numeric(FLERR, arg[4]);
  
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Quit::~FixISPH_Quit"
FixISPH_Quit::~FixISPH_Quit() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Quit::setmask"
int FixISPH_Quit::setmask() {
  FUNCT_ENTER(comm->me);

  int mask = 0;
  mask |= FINAL_INTEGRATE;   
  
  FUNCT_EXIT(comm->me);

  return mask;
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Quit::final_integrate"
void FixISPH_Quit::final_integrate() {
  FUNCT_ENTER(comm->me);

  if (_stat != NULL && _tfinal > 0.0) {
    const double tnow = _stat->vector[0];
    if (tnow > _tfinal) {
      FUNCT_EXIT(comm->me);
      MPI_Finalize();
      exit(0);
    }
  }
  
  FUNCT_EXIT(comm->me);
}
