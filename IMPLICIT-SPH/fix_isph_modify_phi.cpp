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

#include "fix_isph_modify_phi.h"

using namespace LAMMPS_NS;
using namespace FixConst;
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyPhi::FixISPH_ModifyPhi"
FixISPH_ModifyPhi::FixISPH_ModifyPhi(LAMMPS *lmp, int narg, char **arg) 
  : Fix(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Fix ipsh/modify/phi command requires a matching atom_style, i.e., isph");

  if (narg < 3)
    error->all(FLERR,"Illegal number of arguments for fix isph/modify/phi command");

  _phi = NULL;  
  _eapp[0] = 0.0; _eapp[1] = 0.0; _eapp[2] = 0.0;

  int cnt = 0;
  if (strcmp(arg[3], "null") == 0) {
    if (narg < 6)
      error->all(FLERR,"Illegal number of arguments for fix isph/modify/phi command");

    _eapp[0] = force->numeric(FLERR, arg[4]);
    _eapp[1] = force->numeric(FLERR, arg[5]);
    _eapp[2] = force->numeric(FLERR, arg[6]);
    cnt = 7;
  } else {
    const int iphi = modify->find_compute(arg[3]);
    if (iphi >= 0) {
      _phi = modify->compute[iphi];
      string str = "isph/applied-electric-potential";
      if (strncmp(_phi->style, str.c_str(), str.size()) != 0)
        error->all(FLERR,"Input argument for compute style does not match to isph/applied-electric-potential");
    }
    cnt = 4;
  }

  _mask = ~0;
  if (narg > cnt) {
    _mask = 0;
    for (int i=cnt;i<narg;++i)
      _mask |= (1<<(atoi(arg[i])-1));
  }
  
  if (comm->me == 0) {
    if (_phi == NULL) {
      cout << "FixISPH_ModifyPhi:: applied efiled " << _eapp[0] << ", " << _eapp[1] << ", " << _eapp[2] << endl;
    } else {
      cout << "FixISPH_ModifyPhi:: compute phi " << arg[3] << ", " << _phi->style << endl;
    }
    cout << "                    mask " << _mask << endl; 
  }

  time_integrate = 0;
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyPhi::~FixISPH_ModifyPhi"
FixISPH_ModifyPhi::~FixISPH_ModifyPhi() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyPhi::setmask"
int FixISPH_ModifyPhi::setmask() {
  FUNCT_ENTER(comm->me);

  int mask = 0;

  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;

  FUNCT_EXIT(comm->me);

  return mask;
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyPhi::initial_integrate"
void FixISPH_ModifyPhi::initial_integrate(int dummy) {
  FUNCT_ENTER(comm->me);
  
  // get accesses to data structures
  const int nsize = (atom->nlocal + atom->nghost);
  
  if (_phi == NULL) {
    double **x = atom->x;
    int *type = atom->type;
    double *phi = atom->phi;
    
    for (int i=0;i<nsize;++i) 
      if (_mask & (1<<(type[i]-1))) {
        const double *x_at_i = x[i];
        phi[i] = -1.0*(_eapp[0]*x_at_i[0] + 
                       _eapp[1]*x_at_i[1] + 
                       _eapp[2]*x_at_i[2]);
      }
  } else {
    _phi->compute_peratom();
  }
  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_ModifyPhi::final_integrate"
void FixISPH_ModifyPhi::final_integrate() {
  FUNCT_ENTER(comm->me);

  if (_phi != NULL) 
    _phi->compute_scalar();

  FUNCT_EXIT(comm->me);
}

