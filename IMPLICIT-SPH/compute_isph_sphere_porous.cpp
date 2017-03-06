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
#include <fstream>

#include "math.h"
#include "utils.h"

#include "pair_isph_corrected.h"
//#include "pair_isph_mls.h"

#include "compute_isph_sphere_porous.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_SpherePorous::ComputeISPH_SpherePorous"
ComputeISPH_SpherePorous::ComputeISPH_SpherePorous(LAMMPS *lmp, int narg, char **arg) 
  : Compute(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Compute ipsh/status command requires a matching atom_style, i.e., isph");

  if (narg < 5)
    error->all(FLERR,"Illegal number of arguments for compute isph command");

  // arg[3] sphere type id
  _stype = atoi(arg[3]);

  // arg[4] porous file : coordinates and radii 
  _sphere_coords = NULL;
  _sphere_radius = NULL;

  { // set the cut length for deleting inside particles
    PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
    _sphere_hollow_gap = sqrt(pair->getDefaultCutSq());
  }
  
  if (comm->me == 0) {
    cout << " - ComputeISPH_SpherePorous - " << endl
         << "   sphere data file = " << arg[4] << endl;
  }

  readSphereCoordinates(arg[4]);
  setParticleGeometry();
  
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_SpherePorous::~ComputeISPH_SpherePorous"
ComputeISPH_SpherePorous::~ComputeISPH_SpherePorous() {
  FUNCT_ENTER(comm->me);
  
  memory->destroy(_sphere_coords);
  memory->destroy(_sphere_radius);
  
  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_SpherePorous::init"
void ComputeISPH_SpherePorous::init() {
  FUNCT_ENTER(comm->me);

  int count = 0;
  for (int i=0;i<modify->ncompute;++i)
    count += (strcmp(modify->compute[i]->style,"isph/sphere/porous") == 0);

  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute isph/sphere/porous");

  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_SpherePorous::readSphereCoordinates"
void ComputeISPH_SpherePorous::readSphereCoordinates(const char *filename) {
  FUNCT_ENTER(comm->me);

  double x[4];
  
  // count spheres    
  _n_spheres = 0;
  for (ifstream ifs(filename);!ifs.eof();++_n_spheres) 
    ifs >> x[0] >> x[1] >> x[2] >> x[3];

  // allocate
  if (_n_spheres) {
    memory->grow(_sphere_coords, _n_spheres, 3, "compute:sphere:sphere_coords");
    memory->grow(_sphere_radius, _n_spheres,    "compute:sphere:sphere_radius");
  }

  // read sphere coordinates
  ifstream ifs(filename);
  for (int i=0;i<_n_spheres;++i) {
    ifs >> x[0] >> x[1] >> x[2] >> x[3];
    
    _sphere_coords[i][0] = x[0];
    _sphere_coords[i][1] = x[1];
    _sphere_coords[i][2] = x[2];
    _sphere_radius[i]    = x[3];
  }
  
  if (comm->me == 0) {
    cout << " - ComputeISPH_SpherePorous - " << endl
         << "   # of spheres in file = " << _n_spheres << endl;
  }

  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_SpherePorous::setParticleGeometry"
void ComputeISPH_SpherePorous::setParticleGeometry() {
  FUNCT_ENTER(comm->me);

  {
    // get accesses to data structures
    const int nsize = (atom->nlocal + atom->nghost);

    double **x = atom->x;
    int *type = atom->type;

    PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
    if (pair != NULL) {
      pair->setUsePartToComputeNormal(true);
      int *part = atom->part;

      // check all particles including ghost particles
      for (int i=0;i<nsize;++i) {
        int pid = 0;
        double *x_at_i = x[i];
        if (is_coords_in_spheres(x_at_i, pid)) { 
          //type[i] = (pid ? _stype : 0); part[i] = pid; // inside spheres
          type[i] = _stype; 
        } 
        // else {
        //   type[i] = 1; // set with dummy for now
        // }
      }

      // if type == 0 they should be erased
    } else {
      error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH"); 
    }
  }
  
  FUNCT_EXIT(comm->me);  
}

