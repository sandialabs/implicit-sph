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
//#include "pair_isph_mls.h"

#include "fix_isph_shift.h"

using namespace LAMMPS_NS;
using namespace FixConst;
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Shift::FixISPH_Shift"
FixISPH_Shift::FixISPH_Shift(LAMMPS *lmp, int narg, char **arg) 
  : Fix(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Fix ipsh/shift command requires a matching atom_style, i.e., isph");

  // default values
  _shift = 0.005;
  
  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair != NULL) 
    _shiftcut = sqrt(pair->getDefaultCutSq()); // default is fluid/fluid
  else 
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");

  _nonfluidweight = 0.25;
  
  // over-ride shifting parameters
  switch (narg) {
  case 6: {
    const double tmp = force->numeric(FLERR, arg[5]);
    _shiftcut = (tmp < 0.0 ? _shiftcut : tmp);
  }
  case 5: {
    const double tmp = force->numeric(FLERR, arg[4]);
    _nonfluidweight = (tmp < 0.0 ? _nonfluidweight : tmp);
  }
  case 4: {
    const double tmp = force->numeric(FLERR, arg[3]);
    _shift = (tmp < 0.0 ? _shift : tmp);
    break;
  }
  }

  time_integrate = 0;
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Shift::~FixISPH_Shift"
FixISPH_Shift::~FixISPH_Shift() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Shift::setmask"
int FixISPH_Shift::setmask() {
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
#define __FUNCT__ "FixISPH_Shift::initial_integrate"
void FixISPH_Shift::initial_integrate(int dummy) {
  FUNCT_ENTER(comm->me);

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair != NULL) {
    
    if (pair->isAleEnabled()) { 
      pair->computePre();
      pair->shiftParticles(_shift,_shiftcut,_nonfluidweight);
      pair->refreshParticles();
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
#define __FUNCT__ "FixISPH_Shift::final_integrate"
void FixISPH_Shift::final_integrate() {
  FUNCT_ENTER(comm->me);

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair != NULL) {
    if (!pair->isAleEnabled()) { 
      pair->refreshParticles();   
      pair->computePre();
      pair->shiftParticles(_shift,_shiftcut,_nonfluidweight);
    } else {
      error->all(FLERR, "Attempt to shift particles in final_integrate while ALE is enabled");
    }
  } else {
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");
  }

  FUNCT_EXIT(comm->me);
}


