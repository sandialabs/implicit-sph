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
#include "utils.h"

#include "pair_isph_corrected.h"
//#include "pair_isph_mls.h"

#include "compute_isph_concentration.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_Concentration::ComputeISPH_Concentration"
ComputeISPH_Concentration::ComputeISPH_Concentration(LAMMPS *lmp, int narg, char **arg)
  : Compute(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Compute ipsh/concentration command requires a matching atom_style, i.e., isph");

  if (narg < 3)
    error->all(FLERR, "Illegal number of arguments for compute isph command");

  _cid = (narg == 3 ? atoi(arg[3]) : 0);

  if (_cid >= ISPH_MAX_CONCENTRATION)
    error->all(FLERR, "Concentration ID is bigger than its max ID");

  peratom_flag = 1;
  size_peratom_cols = 0; // 0 = single vector

  vector_atom = atom->concentration[_cid];

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_Concentration::~ComputeISPH_Concentration"
ComputeISPH_Concentration::~ComputeISPH_Concentration() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_Concentration::init"
void ComputeISPH_Concentration::init() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_Concentration::compute_peratom"
void ComputeISPH_Concentration::compute_peratom() {
  FUNCT_ENTER(comm->me);

  invoked_peratom = update->ntimestep;

  // this array can grow
  vector_atom = atom->concentration[_cid];

  FUNCT_EXIT(comm->me);
}
