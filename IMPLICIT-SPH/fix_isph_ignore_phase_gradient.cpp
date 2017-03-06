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
#include "domain.h"
#include "math.h"

#include <iostream>
#include "utils.h"

#include "pair_isph_corrected.h"
//#include "pair_isph_mls.h"

#include "fix_isph_ignore_phase_gradient.h"

using namespace LAMMPS_NS;
using namespace FixConst;
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_IgnorePhaseGradient::FixISPH_IgnorePhaseGradient"
FixISPH_IgnorePhaseGradient::FixISPH_IgnorePhaseGradient(LAMMPS *lmp, int narg, char **arg) 
  : Fix(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Fix ipsh/ignore/phasegradient command requires a matching atom_style, i.e., isph");

  if (narg < 6)
    error->all(FLERR,"Illegal number of arguments for fix isph/modify/type command");
  
  if      (strcmp(arg[3],"x") == 0) _axis = 0;
  else if (strcmp(arg[3],"y") == 0) _axis = 1;
  else if (strcmp(arg[3],"z") == 0) _axis = 2;

  _pt_on_axis = force->numeric(FLERR, arg[4]);
  _thres_over_cut = force->numeric(FLERR, arg[5]);

  if (comm->me == 0) {
    cout << "FixISPH_IgnorePhaseGradient:: " << endl
         << "    axis               = " << _axis << endl
         << "    point on axis      = " << _pt_on_axis << endl 
         << "    threshold over cut = " << _thres_over_cut << endl;
  }

  time_integrate = 0;
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_IgnorePhaseGradient::~FixISPH_IgnorePhaseGradient"
FixISPH_IgnorePhaseGradient::~FixISPH_IgnorePhaseGradient() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_IgnorePhaseGradient::setmask"
int FixISPH_IgnorePhaseGradient::setmask() {
  FUNCT_ENTER(comm->me);

  int mask = 0;

  FUNCT_EXIT(comm->me);

  return mask;
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_IgnorePhaseGradient::ignorePhaseGradient"
void FixISPH_IgnorePhaseGradient::ignorePhaseGradient(double **grad, const double cut) const {
  FUNCT_ENTER(comm->me);

  const double cut_applied = cut*_thres_over_cut;

  const double val_min = _pt_on_axis - cut_applied;
  const double val_max = _pt_on_axis + cut_applied;

  const int nsize = (atom->nlocal + atom->nghost);

  double **x = atom->x;
  
  for (int i=0;i<nsize;++i) {
    double *x_at_i = x[i];
    if (val_min < x_at_i[_axis] && x_at_i[_axis] < val_max)
      memset(grad[i], 0, sizeof(double)*3);
  }
  
  FUNCT_EXIT(comm->me);
}


