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

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

#include "pair_isph_corrected.h"
#include "pair_isph_mls.h"

#include "fix_isph_analytic.h"

using namespace LAMMPS_NS;
using namespace FixConst;
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Analytic::FixISPH_Analytic"
FixISPH_Analytic::FixISPH_Analytic(LAMMPS *lmp, int narg, char **arg)
  : Fix(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Fix isph/analytic command requires a matching atom_style, i.e., isph");

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair == NULL) {
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");
  }

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Analytic::~FixISPH_Analytic"
FixISPH_Analytic::~FixISPH_Analytic() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Shift::setmask"
int FixISPH_Analytic::setmask() {
  FUNCT_ENTER(comm->me);

  int mask = 0;

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair != NULL) {
    if (pair->isAleEnabled()) {
      // ALE (aka velocity correction) move/shift before pair->compute()
      mask |= INITIAL_INTEGRATE;
    } else {
      // pressure correction move/shift it pair->compute()
      mask |= FINAL_INTEGRATE;
    }
  } else {
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");
  }

  FUNCT_EXIT(comm->me);

  return mask;
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Analytic::initial_integrate"
void FixISPH_Analytic::initial_integrate(int dummy) {
  FUNCT_ENTER(comm->me);

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair != NULL) {
    
    if (pair->isAleEnabled()) {

      int nlocal = atom->nlocal;
      double **_X = atom->x;
      double **v;
      double t = update->dt*update->ntimestep;

      v = atom->v;

      for (int i=0;i<nlocal;++i) {

        if (pair->getParticleKind(type[i]) == PairISPH_MLS::Fluid)
          continue;

        double x(_X[i][0]), y(_X[i][1]), z(_X[i][2]);

        v[i][0] = sin(t)*sin(x)*(cos(y)-cos(z));
        v[i][1] = -sin(t)*sin(y)*(cos(x)-cos(z));
        v[i][2] = sin(t)*sin(z)*(cos(x)-cos(y));
      }
    } else {
      error->all(FLERR, "Attempt to shift particles in initial_integrate while ALE is disabled");
    }
  } else {
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");
  }

  FUNCT_EXIT(comm->me);
}


/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Analytic::final_integrate"
void FixISPH_Analytic::final_integrate() {
  FUNCT_ENTER(comm->me);

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair != NULL) {
    if (!pair->isAleEnabled()) {

      int nlocal = atom->nlocal;
      int *type = atom->type;
      double **_X = atom->x;

      double **v;
      double t = update->dt*update->ntimestep;

      v = pair->getVstar();

      for (int i=0;i<nlocal;++i) {
        if (pair->getParticleKind(type[i]) == PairISPH_MLS::Fluid)
          continue;

        double x(_X[i][0]), y(_X[i][1]), z(_X[i][2]);
        v[i][0] = sin(t)*sin(x)*(cos(y)-cos(z));
        v[i][1] = -sin(t)*sin(y)*(cos(x)-cos(z));
        v[i][2] = sin(t)*sin(z)*(cos(x)-cos(y));
      }

    } else {
      error->all(FLERR, "Attempt to shift particles in final_integrate while ALE is enabled");
    }
  } else {
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");
  }

  FUNCT_EXIT(comm->me);
}


