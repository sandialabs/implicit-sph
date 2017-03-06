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

#include "compute_isph_velocity_divergence.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_VelocityDivergence::ComputeISPH_VelocityDivergence"
ComputeISPH_VelocityDivergence::ComputeISPH_VelocityDivergence(LAMMPS *lmp, int narg, char **arg) 
  : Compute(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Compute ipsh/velocity-curl command requires a matching atom_style, i.e., isph");

  if (narg < 3)
    error->all(FLERR,"Illegal number of arguments for compute isph command");

  _vflag = 0;
  for (int i=3;i<narg;++i)
    if (strcmp(arg[i],"vstar") == 0) _vflag = 1;
  
  peratom_flag = 1;
  size_peratom_cols = 0; // 0 = vector, N = columns in peratom array

  _nmax = 0;
  _div = NULL;

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_VelocityDivergence::~ComputeISPH_VelocityDivergence"
ComputeISPH_VelocityDivergence::~ComputeISPH_VelocityDivergence() {
  FUNCT_ENTER(comm->me);

  memory->destroy(_div);

  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_VelocityDivergence::init"
void ComputeISPH_VelocityDivergence::init() {
  FUNCT_ENTER(comm->me);

  int count = 0;
  for (int i=0;i<modify->ncompute;++i)
    count += (strcmp(modify->compute[i]->style,"isph/velocity-curl") == 0);

  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute isph/velocity-curl");

  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_VelocityDivergence::compute_peratom"
void ComputeISPH_VelocityDivergence::compute_peratom() {
  FUNCT_ENTER(comm->me);

  if (invoked_peratom != update->ntimestep) {
    invoked_peratom = update->ntimestep;
    
    // grow curl array if necessary
    if (atom->nlocal > _nmax) {
      _nmax = atom->nmax;
      memory->grow(_div, _nmax, "compute:velocity-div:div");
      
      vector_atom = _div;
    }
    
    PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
    if (pair == NULL)
      error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");
    
    double **v = (_vflag ? pair->getVstar() : atom->v); 
    
    pair->computeVelocityDivergence(v, _div);
  }

  FUNCT_EXIT(comm->me);  
}

