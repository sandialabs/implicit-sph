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

/* ----------------------------------------------------------------------
 Contributing author: Mike Parks (SNL), Nathaniel Trask (Brown)
 ------------------------------------------------------------------------- */

#include "float.h"
#include "stdlib.h"
#include "atom_vec_isph.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "math.h"

#include "macrodef.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecISPH::AtomVecISPH(LAMMPS *lmp) : AtomVec(lmp) {
  molecular = 1;
  bonds_allow = 1;
  mass_type = 0;

  comm_x_only = 0;
  comm_f_only = 0;

  // number of doubles to be communicated
  // x[3], v[3], vfrac, pressure, rmass, psi, phi
  // xdot[3],
  // vprev[ISPH_BDF_MAX_ORDER][3], xprev[ISPH_BDF_MAX_ORDER][3]
  // concentration[ISPH_MAX_CONCENTRATION]
  size_forward = 11 + 3 + 6*ISPH_BDF_MAX_ORDER + ISPH_MAX_CONCENTRATION;
  size_reverse =  3; // atom bond need this

  // tag, type, mask + part, viscosity, density, eps, sigma + (size_forward) + molecule
  size_border = size_forward + 8 + 1;

  size_velocity = 0; // we don't need this

  size_data_atom = 11;  // 8 + 3
  size_data_vel  = 4; // I do not know this.
  xcol_data      = 9;  // location starting coordinates in input file (last 3 are always coordinates)

  // style flag; used in sanity check
  atom->isph_flag = 1;
  atom->molecule_flag = 1;

  // additional variables that need to be handled in LAMMPS
  // note that velocity and position are handled without flag by default in LAMMPS
  atom->pressure_flag = 1;

  atom->eps_flag      = 1;
  atom->sigma_flag    = 1;

  atom->psi_flag      = 1;
  atom->psi0_flag     = 1;
  atom->psigrad_flag  = 1;

  atom->phi_flag      = 1;
  atom->phigrad_flag  = 1;

  atom->concentration_flag = 1;

  // material parameters
  atom->vfrac_flag     = 1;
  atom->rmass_flag     = 1;
  atom->viscosity_flag = 1;

  atom->part_flag      = 1;
}

/* ----------------------------------------------------------------------
 grow atom arrays
 n = 0 grows arrays by DELTA
 n > 0 allocates arrays to size n
 ------------------------------------------------------------------------- */

void AtomVecISPH::grow(int n) {
  if (n == 0)
    nmax += DELTA;
  else
    nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR, "Per-processor system is too big");

  tag   = memory->grow(atom->tag,   nmax, "atom:tag");
  type  = memory->grow(atom->type,  nmax, "atom:type");
  mask  = memory->grow(atom->mask,  nmax, "atom:mask");
  image = memory->grow(atom->image, nmax, "atom:image");

  part  = memory->grow(atom->part,  nmax, "atom:part");

  x       = memory->grow(atom->x,       nmax, 3, "atom:x");
  v       = memory->grow(atom->v,       nmax, 3, "atom:v");

  xdot    = memory->grow(atom->xdot,    nmax, 3, "atom:xdot");

  vprev   = memory->grow(atom->vprev,   nmax, ISPH_BDF_MAX_ORDER, 3, "atom:vprev");
  xprev   = memory->grow(atom->xprev,   nmax, ISPH_BDF_MAX_ORDER, 3, "atom:xprev");

  f       = memory->grow(atom->f,       nmax, 3, "atom:f");

  rmass     = memory->grow(atom->rmass,     nmax, "atom:rmass");
  viscosity = memory->grow(atom->viscosity, nmax, "atom:viscosity");
  vfrac     = memory->grow(atom->vfrac,     nmax, "atom:vfrac");
  pressure  = memory->grow(atom->pressure,  nmax, "atom:pressure");
  density   = memory->grow(atom->density,   nmax, "atom:density");

  eps       = memory->grow(atom->eps,       nmax, "atom:eps");
  sigma     = memory->grow(atom->sigma,     nmax, "atom:sigma");

  psi       = memory->grow(atom->psi,       nmax, "atom:psi");
  psi0      = memory->grow(atom->psi0,      nmax, "atom:psi0");
  psigrad   = memory->grow(atom->psigrad,   nmax, 3, "atom:psigrad");

  phi       = memory->grow(atom->phi,       nmax, "atom:phi");
  phigrad   = memory->grow(atom->phigrad,   nmax, 3, "atom:phigrad");

  // multivector
  concentration = memory->grow(atom->concentration,   ISPH_MAX_CONCENTRATION, nmax, "atom:concentration");

  // molecule 
  molecule = memory->grow(atom->molecule,nmax,"atom:molecule");

  nspecial = memory->grow(atom->nspecial,nmax,3,"atom:nspecial");
  special = memory->grow(atom->special,nmax,atom->maxspecial,"atom:special");

  num_bond = memory->grow(atom->num_bond,nmax,"atom:num_bond");
  bond_type = memory->grow(atom->bond_type,nmax,atom->bond_per_atom,
                           "atom:bond_type");
  bond_atom = memory->grow(atom->bond_atom,nmax,atom->bond_per_atom,
                           "atom:bond_atom");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
 reset local array ptrs
 ------------------------------------------------------------------------- */

void AtomVecISPH::grow_reset() {
  tag   = atom->tag;
  type  = atom->type;
  mask  = atom->mask;
  image = atom->image;

  part  = atom->part;

  x    = atom->x;
  v    = atom->v;

  xdot  = atom->xdot;

  vprev = atom->vprev;
  xprev = atom->xprev;

  f    = atom->f;

  vfrac     = atom->vfrac;
  rmass     = atom->rmass;
  viscosity = atom->viscosity;
  pressure  = atom->pressure;
  density   = atom->density;

  eps       = atom->eps;
  sigma     = atom->sigma;

  psi       = atom->psi;
  psi0      = atom->psi0;
  psigrad   = atom->psigrad;

  phi       = atom->phi;
  phigrad   = atom->phigrad;

  concentration = atom->concentration;

  molecule = atom->molecule;
  nspecial = atom->nspecial; 
  special = atom->special;
  num_bond = atom->num_bond; 
  bond_type = atom->bond_type;
  bond_atom = atom->bond_atom;
}

/* ----------------------------------------------------------------------
 copy atom I info to atom J
 ------------------------------------------------------------------------- */

void AtomVecISPH::copy(int i, int j, int delflag) {
  tag[j]   = tag[i];
  type[j]  = type[i];
  mask[j]  = mask[i];
  image[j] = image[i];

  part[j]  = part[i];

  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];

  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  xdot[j][0] = xdot[i][0];
  xdot[j][1] = xdot[i][1];
  xdot[j][2] = xdot[i][2];

  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    xprev[j][k][0] = xprev[i][k][0];
    xprev[j][k][1] = xprev[i][k][1];
    xprev[j][k][2] = xprev[i][k][2];
  }

  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    vprev[j][k][0] = vprev[i][k][0];
    vprev[j][k][1] = vprev[i][k][1];
    vprev[j][k][2] = vprev[i][k][2];
  }

  vfrac[j]     = vfrac[i];
  rmass[j]     = rmass[i];
  pressure[j]  = pressure[i];
  viscosity[j] = viscosity[i];
  density[j]   = density[i];

  eps[j]   = eps[i];
  sigma[j] = sigma[i];

  psi0[j] = psi0[i];
  psi[j]  = psi[i];

  psigrad[j][0] = psigrad[i][0];
  psigrad[j][1] = psigrad[i][1];
  psigrad[j][2] = psigrad[i][2];

  phi[j] = phi[i];

  phigrad[j][0] = phigrad[i][0];
  phigrad[j][1] = phigrad[i][1];
  phigrad[j][2] = phigrad[i][2];

  for (int k=0;k<ISPH_MAX_CONCENTRATION;++k) 
    concentration[k][j] = concentration[k][i];
  
  molecule[j] = molecule[i];
  
  num_bond[j] = num_bond[i];
  for (int k = 0; k < num_bond[j]; k++) {
    bond_type[j][k] = bond_type[i][k];
    bond_atom[j][k] = bond_atom[i][k];
  }
  
  nspecial[j][0] = nspecial[i][0];
  nspecial[j][1] = nspecial[i][1];
  nspecial[j][2] = nspecial[i][2];
  for (int k = 0; k < nspecial[j][2]; k++) 
    special[j][k] = special[i][k];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i, j, delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecISPH::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, m;
  double dx, dy, dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];

      // coords
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];

      // velocity
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];

      // ale vector
      buf[m++] = xdot[j][0];
      buf[m++] = xdot[j][1];
      buf[m++] = xdot[j][2];

      // position history
      for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
        buf[m++] = xprev[j][k][0];
        buf[m++] = xprev[j][k][1];
        buf[m++] = xprev[j][k][2];
      }

      // velocity history
      for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
        buf[m++] = vprev[j][k][0];
        buf[m++] = vprev[j][k][1];
        buf[m++] = vprev[j][k][2];
      }

      // pressure
      buf[m++] = pressure[j];

      // volume
      buf[m++] = vfrac[j];

      // mass
      buf[m++] = rmass[j];

      // psi
      buf[m++] = psi[j];

      // phi
      buf[m++] = phi[j];

      // concentration
      for (int k=0;k<ISPH_MAX_CONCENTRATION;++k) 
        buf[m++] = concentration[k][j];

    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
      dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
      dz = pbc[2] * domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];

      // coords
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;

      // velocity
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];

      // ale vecvtor
      buf[m++] = xdot[j][0];
      buf[m++] = xdot[j][1];
      buf[m++] = xdot[j][2];

      // position history
      for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
        buf[m++] = xprev[j][k][0];
        buf[m++] = xprev[j][k][1];
        buf[m++] = xprev[j][k][2];
      }

      // velocity history
      for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
        buf[m++] = vprev[j][k][0];
        buf[m++] = vprev[j][k][1];
        buf[m++] = vprev[j][k][2];
      }

      // pressure
      buf[m++] = pressure[j];

      // volume
      buf[m++] = vfrac[j];

      // mass
      buf[m++] = rmass[j];

      // psi
      buf[m++] = psi[j];

      // phi
      buf[m++] = phi[j];

      for (int k=0;k<ISPH_MAX_CONCENTRATION;++k) 
        buf[m++] = concentration[k][j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecISPH::pack_comm_vel(int n, int *list, double *buf, int pbc_flag, int *pbc)

{
  error->all(FLERR, "AtomVecISPH::pack_comm_vel - not necessary");
  return -1;
}

void AtomVecISPH::unpack_comm(int n, int first, double *buf) {
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax)
      grow(0);

    // coords
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];

    //velocity
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];

    // ale vector
    xdot[i][0] = buf[m++];
    xdot[i][1] = buf[m++];
    xdot[i][2] = buf[m++];

    // position history
    for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
      xprev[i][k][0] = buf[m++];
      xprev[i][k][1] = buf[m++];
      xprev[i][k][2] = buf[m++];
    }

    // velocity history
    for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
      vprev[i][k][0] = buf[m++];
      vprev[i][k][1] = buf[m++];
      vprev[i][k][2] = buf[m++];
    }

    //pressure
    pressure[i] = buf[m++];

    // volume
    vfrac[i] = buf[m++];

    // rmass (used in lammps core)
    rmass[i] = buf[m++];

    // psi; electric potential
    psi[i] = buf[m++];

    // psi; electric potential
    phi[i] = buf[m++];

    for (int k=0;k<ISPH_MAX_CONCENTRATION;++k)
      concentration[k][i] = buf[m++];

  }

}

void AtomVecISPH::unpack_comm_vel(int n, int first, double *buf) {
  error->all(FLERR, "AtomVecISPH::unpack_comm_vel - not necessary");
}

/* ---------------------------------------------------------------------- */

int AtomVecISPH::pack_reverse(int n, int first, double *buf) {
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecISPH::unpack_reverse(int n, int *list, double *buf) {
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecISPH::pack_border(int n, int *list, double *buf, int pbc_flag, int *pbc) {
  int i, j, m;
  double dx, dy, dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];

      // position
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];

      // velocity
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];

      // ale vector
      buf[m++] = xdot[j][0];
      buf[m++] = xdot[j][1];
      buf[m++] = xdot[j][2];

      // position history
      for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
        buf[m++] = xprev[j][k][0];
        buf[m++] = xprev[j][k][1];
        buf[m++] = xprev[j][k][2];
      }

      // velocity history
      for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
        buf[m++] = vprev[j][k][0];
        buf[m++] = vprev[j][k][1];
        buf[m++] = vprev[j][k][2];
      }

      // pressure
      buf[m++] = pressure[j];

      // psi; electric potential
      buf[m++] = psi[j];

      // phi; electric potential for background field
      buf[m++] = phi[j];

      // concentration
      for (int k=0;k<ISPH_MAX_CONCENTRATION;++k) 
        buf[m++] = concentration[k][j];

      // eps; dielectric constant
      buf[m++] = eps[j];

      // sigma; conductivity constant
      buf[m++] = sigma[j];

      // volume
      buf[m++] = vfrac[j];

      //viscosity
      buf[m++] = viscosity[j];

      //density
      buf[m++] = density[j];

      //mass
      buf[m++] = rmass[j];

      //part
      buf[m++] = ubuf(part[j]).d;

      //tag
      buf[m++] = ubuf(tag[j]).d;

      //type
      buf[m++] = ubuf(type[j]).d;

      //mask
      buf[m++] = ubuf(mask[j]).d;

      // molecule
      buf[m++] = ubuf(molecule[j]).d;
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      printf("we cannot deal with triclinic case yet!!!\n");
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];

      // position
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;

      //velocity
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];

      // ale vector
      buf[m++] = xdot[j][0];
      buf[m++] = xdot[j][1];
      buf[m++] = xdot[j][2];

      // position history
      for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
        buf[m++] = xprev[j][k][0];
        buf[m++] = xprev[j][k][1];
        buf[m++] = xprev[j][k][2];
      }

      // velocity history
      for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
        buf[m++] = vprev[j][k][0];
        buf[m++] = vprev[j][k][1];
        buf[m++] = vprev[j][k][2];
      }

      //pressure
      buf[m++] = pressure[j];

      // psi; electric potential
      buf[m++] = psi[j];

      // phi; electric potential
      buf[m++] = phi[j];

      // concentration
      for (int k=0;k<ISPH_MAX_CONCENTRATION;++k)
        buf[m++] = concentration[k][j];

      // eps; electric potential
      buf[m++] = eps[j];

      // sigma; conductivity
      buf[m++] = sigma[j];

      // volume
      buf[m++] = vfrac[j];

      //viscosity
      buf[m++] = viscosity[j];

      //density
      buf[m++] = density[j];

      //mass
      buf[m++] = rmass[j];

      //part
      buf[m++] = ubuf(part[j]).d;

      //tag
      buf[m++] = ubuf(tag[j]).d;

      //type
      buf[m++] = ubuf(type[j]).d;

      //mask
      buf[m++] = ubuf(mask[j]).d;

      // molecule
      buf[m++] = ubuf(molecule[j]).d;
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecISPH::pack_border_vel(int n, int *list, double *buf, int pbc_flag, int *pbc) {
  error->all(FLERR, "AtomVecISPH::pack_border_vel - not necessary");
  return -1;
}

/* ---------------------------------------------------------------------- */

void AtomVecISPH::unpack_border(int n, int first, double *buf) {
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax)
      grow(0);

    // position
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];

    // velocity
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];

    // ale vector
    xdot[i][0] = buf[m++];
    xdot[i][1] = buf[m++];
    xdot[i][2] = buf[m++];

    // position history
    for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
      xprev[i][k][0] = buf[m++];
      xprev[i][k][1] = buf[m++];
      xprev[i][k][2] = buf[m++];
    }

    // velocity history
    for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
      vprev[i][k][0] = buf[m++];
      vprev[i][k][1] = buf[m++];
      vprev[i][k][2] = buf[m++];
    }

    // pressure
    pressure[i] = buf[m++];

    // electric potential
    psi[i] = buf[m++];

    // electric potential
    phi[i] = buf[m++];

    // concentration
    for (int k=0;k<ISPH_MAX_CONCENTRATION;++k)
      concentration[k][i] = buf[m++];
    
    // eps; electric potential
    eps[i] = buf[m++];

    // sigma; electric potential
    sigma[i] = buf[m++];

    // volume
    vfrac[i] = buf[m++];

    // viscosity
    viscosity[i] = buf[m++];

    // density
    density[i] = buf[m++];

    // mass
    rmass[i] = buf[m++];

    //part
    part[i] = (int) ubuf(buf[m++]).i;

    // tag
    tag[i] = (tagint) ubuf(buf[m++]).i;

    //type
    type[i] = (int) ubuf(buf[m++]).i;

    // mask
    mask[i] = (int) ubuf(buf[m++]).i;

    // molecule
    molecule[i] = (tagint) ubuf(buf[m++]).i;
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecISPH::unpack_border_vel(int n, int first, double *buf) {
  error->all(FLERR, "AtomVecISPH::unpack_border_vel - not necessary");
}

/* ----------------------------------------------------------------------
 pack data for atom I for sending to another proc
 xyz must be 1st 3 values, so comm::exchange() can test on them
 ------------------------------------------------------------------------- */

int AtomVecISPH::pack_exchange(int i, double *buf) {
  int m = 1;

  // position
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];

  //velocity
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  // ale vector
  buf[m++] = xdot[i][0];
  buf[m++] = xdot[i][1];
  buf[m++] = xdot[i][2];

  // position history
  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    buf[m++] = xprev[i][k][0];
    buf[m++] = xprev[i][k][1];
    buf[m++] = xprev[i][k][2];
  }

  // velocity history
  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    buf[m++] = vprev[i][k][0];
    buf[m++] = vprev[i][k][1];
    buf[m++] = vprev[i][k][2];
  }

  //force
  // buf[m++] = f[i][0];
  // buf[m++] = f[i][1];
  // buf[m++] = f[i][2];

  // pressure
  buf[m++] = pressure[i];

  // electric potential
  buf[m++] = psi[i];

  // concentration
  for (int k=0;k<ISPH_MAX_CONCENTRATION;++k)
    buf[m++] = concentration[k][i];

  // mass
  buf[m++] = rmass[i];

  // viscosity
  buf[m++] = viscosity[i];

  // volume
  buf[m++] = vfrac[i];

  // density
  buf[m++] = density[i];

  //part
  buf[m++] = ubuf(part[i]).d;

  //tag
  buf[m++] = ubuf(tag[i]).d;

  //type
  buf[m++] = ubuf(type[i]).d;

  //mask
  buf[m++] = ubuf(mask[i]).d;

  //image
  buf[m++] = ubuf(image[i]).d;

  // molecule
  buf[m++] = ubuf(molecule[i]).d;
  buf[m++] = ubuf(num_bond[i]).d;
  for (int k = 0; k < num_bond[i]; k++) {
    buf[m++] = ubuf(bond_type[i][k]).d;
    buf[m++] = ubuf(bond_atom[i][k]).d;
  }

  buf[m++] = ubuf(nspecial[i][0]).d;
  buf[m++] = ubuf(nspecial[i][1]).d;
  buf[m++] = ubuf(nspecial[i][2]).d;
  for (int k = 0; k < nspecial[i][2]; k++) 
    buf[m++] = ubuf(special[i][k]).d;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i, &buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecISPH::unpack_exchange(double *buf) {
  int nlocal = atom->nlocal;
  if (nlocal == nmax)
    grow(0);

  int m = 1;

  // coord
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];

  // velocity
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  // ale vector
  xdot[nlocal][0] = buf[m++];
  xdot[nlocal][1] = buf[m++];
  xdot[nlocal][2] = buf[m++];

  // position history
  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    xprev[nlocal][k][0] = buf[m++];
    xprev[nlocal][k][1] = buf[m++];
    xprev[nlocal][k][2] = buf[m++];
  }

  // velocity history
  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    vprev[nlocal][k][0] = buf[m++];
    vprev[nlocal][k][1] = buf[m++];
    vprev[nlocal][k][2] = buf[m++];
  }

  // f[nlocal][0] = buf[m++];
  // f[nlocal][1] = buf[m++];
  // f[nlocal][2] = buf[m++];

  pressure[nlocal] = buf[m++];
  psi[nlocal] = buf[m++];

  for (int k=0;k<ISPH_MAX_CONCENTRATION;++k)
    concentration[k][nlocal] = buf[m++];

  rmass[nlocal] = buf[m++];
  viscosity[nlocal] = buf[m++];
  vfrac[nlocal] = buf[m++];

  density[nlocal] = buf[m++];

  part[nlocal] = (int) ubuf(buf[m++]).i;

  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  molecule[nlocal] = (tagint) ubuf(buf[m++]).i;

  num_bond[nlocal] = (int) ubuf(buf[m++]).i;
  for (int k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    bond_atom[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  nspecial[nlocal][0] = (int) ubuf(buf[m++]).i;
  nspecial[nlocal][1] = (int) ubuf(buf[m++]).i;
  nspecial[nlocal][2] = (int) ubuf(buf[m++]).i;
  for (int k = 0; k < nspecial[nlocal][2]; k++)
    special[nlocal][k] = (tagint) ubuf(buf[m++]).i;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->unpack_exchange(nlocal, &buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
 size of restart data for all atoms owned by this proc
 include extra data stored by fixes
 ------------------------------------------------------------------------- */

int AtomVecISPH::size_restart() {
  int i;

  int nlocal = atom->nlocal;
  int n = 20 * nlocal;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
 pack atom I's data for restart file including extra quantities
 xyz must be 1st 3 values, so that read_restart can test on them
 molecular types may be negative, but write as positive
 ------------------------------------------------------------------------- */

int AtomVecISPH::pack_restart(int i, double *buf) {
  int m = 1;

  // position
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];

  // velocity
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  // ale vector
  buf[m++] = xdot[i][0];
  buf[m++] = xdot[i][1];
  buf[m++] = xdot[i][2];

  // position history
  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    buf[m++] = xprev[i][k][0];
    buf[m++] = xprev[i][k][1];
    buf[m++] = xprev[i][k][2];
  }

  // velocity history
  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    buf[m++] = vprev[i][k][0];
    buf[m++] = vprev[i][k][1];
    buf[m++] = vprev[i][k][2];
  }

  // force
  buf[m++] = f[i][0];
  buf[m++] = f[i][1];
  buf[m++] = f[i][2];

  // pressure
  buf[m++] = pressure[i];

  // electric potential
  buf[m++] = psi[i];

  // concentration
  for (int k=0;k<ISPH_MAX_CONCENTRATION;++k)
    buf[m++] = concentration[k][i];

  // volume
  buf[m++] = vfrac[i];

  // mass
  buf[m++] = rmass[i];

  // viscosity
  buf[m++] = viscosity[i];

  // density
  buf[m++] = density[i];

  // part
  part[m++] = ubuf(part[i]).d;

  // tag
  buf[m++] = ubuf(tag[i]).d;

  // type
  buf[m++] = ubuf(type[i]).d;

  // mask
  buf[m++] = ubuf(mask[i]).d;

  // image
  buf[m++] = ubuf(image[i]).d;

  // molecule 
  buf[m++] = ubuf(molecule[i]).d;

  buf[m++] = ubuf(num_bond[i]).d;
  for (int k = 0; k < num_bond[i]; k++) {
    buf[m++] = ubuf(MAX(bond_type[i][k],-bond_type[i][k])).d;
    buf[m++] = ubuf(bond_atom[i][k]).d;
  }

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i, &buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
 unpack data for one atom from restart file including extra quantities
 ------------------------------------------------------------------------- */

int AtomVecISPH::unpack_restart(double *buf) {
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra, nmax, atom->nextra_store, "atom:extra");
  }

  int m = 1;

  // coord
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];

  // velocity
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  // ale vector
  xdot[nlocal][0] = buf[m++];
  xdot[nlocal][1] = buf[m++];
  xdot[nlocal][2] = buf[m++];

  // position history
  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    xprev[nlocal][k][0] = buf[m++];
    xprev[nlocal][k][1] = buf[m++];
    xprev[nlocal][k][2] = buf[m++];
  }

  // velocity history
  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    vprev[nlocal][k][0] = buf[m++];
    vprev[nlocal][k][1] = buf[m++];
    vprev[nlocal][k][2] = buf[m++];
  }

  f[nlocal][0] = buf[m++];
  f[nlocal][1] = buf[m++];
  f[nlocal][2] = buf[m++];

  pressure[nlocal] = buf[m++];
  psi[nlocal] = buf[m++];

  for (int k=0;k<ISPH_MAX_CONCENTRATION;++k)
    concentration[k][nlocal] = buf[m++];

  vfrac[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];

  viscosity[nlocal] = buf[m++];
  density[nlocal] = buf[m++];

  part[nlocal] = (int) ubuf(buf[m++]).i;

  tag[nlocal] = (int) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  molecule[nlocal] = (tagint) ubuf(buf[m++]).i;

  num_bond[nlocal] = (int) ubuf(buf[m++]).i;
  for (int k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    bond_atom[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int>(buf[0]) - m;
    for (int i = 0; i < size; i++)
      extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
 create one atom of itype at coord
 set other values to defaults
 ------------------------------------------------------------------------- */

void AtomVecISPH::create_atom(int itype, double *coord) {
  int nlocal = atom->nlocal;
  if (nlocal == nmax)
    grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;

  // coord
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  // velocity
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  // ale vector
  xdot[nlocal][0] = 0.0;
  xdot[nlocal][1] = 0.0;
  xdot[nlocal][2] = 0.0;

  // position history
  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    xprev[nlocal][k][0] = 0.0;
    xprev[nlocal][k][1] = 0.0;
    xprev[nlocal][k][2] = 0.0;
  }

  // velocity history
  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    vprev[nlocal][k][0] = 0.0;
    vprev[nlocal][k][1] = 0.0;
    vprev[nlocal][k][2] = 0.0;
  }

  f[nlocal][0] = 0.0;
  f[nlocal][1] = 0.0;
  f[nlocal][2] = 0.0;

  pressure[nlocal] = 0.0;

  eps[nlocal] = 0.0;
  sigma[nlocal] = 0.0;

  psi[nlocal] = 0.0;
  phi[nlocal] = 0.0;

  for (int k=0;k<ISPH_MAX_CONCENTRATION;++k) 
    concentration[k][nlocal] = 0.0;

  vfrac[nlocal] = 1.0; // volume will be recalculated in computeVolume
  rmass[nlocal] = 1.0; // mass will be recalculated in computeVolume to beconsistent with volume and density
  viscosity[nlocal] = 0.0;
  density[nlocal] = 1.0;

  // molecule
  molecule[nlocal] = 0;
  num_bond[nlocal] = 0;
  nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
 unpack one line from Atoms section of data file
 initialize other atom quantities
 ------------------------------------------------------------------------- */
#include "iostream"
void AtomVecISPH::data_atom(double *coord, imageint imagetmp, char **values) {
// MLP: Fill out later.
  int nlocal = atom->nlocal;
  if (nlocal == nmax)
    grow(0);

  tag[nlocal] = ATOTAGINT(values[0]); 
  if (tag[nlocal] <= 0)
    error->one(FLERR, "Invalid atom ID in Atoms section of data file");

  type[nlocal] = atoi(values[1]); 
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR, "Invalid atom type in Atoms section of data file");

  density[nlocal] = atof(values[2]); 
  if (density[nlocal] <= 0.0)
    error->one(FLERR, "Invalid density value");

  viscosity[nlocal] = atof(values[3]); 
  if (viscosity[nlocal] <= 0.0)
    error->one(FLERR, "Invalid viscosity value");

  pressure[nlocal] = atof(values[4]); 
  psi0[nlocal] = atof(values[5]); 

  eps[nlocal] = atof(values[6]); 
  if (eps[nlocal] <= 0.0)
    error->one(FLERR, "Invalid dielectric constant value");

  molecule[nlocal] = ATOTAGINT(values[7]); 

  // do not specify sigma for now

  // coord
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  // velocity
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  mask[nlocal] = 1;

  image[nlocal] = imagetmp;

  // ale vector
  xdot[nlocal][0] = 0.0;
  xdot[nlocal][1] = 0.0;
  xdot[nlocal][2] = 0.0;

  // position history
  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    xprev[nlocal][k][0] = 0.0;
    xprev[nlocal][k][1] = 0.0;
    xprev[nlocal][k][2] = 0.0;
  }

  // velocity history
  for (int k=0;k<ISPH_BDF_MAX_ORDER;++k) {
    vprev[nlocal][k][0] = 0.0;
    vprev[nlocal][k][1] = 0.0;
    vprev[nlocal][k][2] = 0.0;
  }

  f[nlocal][0] = 0.0;
  f[nlocal][1] = 0.0;
  f[nlocal][2] = 0.0;

  vfrac[nlocal] = 1.0; // volume will be recalculated in computeVolume
  rmass[nlocal] = 1.0; // mass will be recalculated in computeVolume to be consistent with volume and density

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
 return # of bytes of allocated memory
 ------------------------------------------------------------------------- */

bigint AtomVecISPH::memory_usage() {
// MLP: Fill out later.
  bigint bytes = 0;

  if (atom->memcheck("tag"))
    bytes += memory->usage(tag, nmax);
  if (atom->memcheck("type"))
    bytes += memory->usage(type, nmax);
  if (atom->memcheck("mask"))
    bytes += memory->usage(mask, nmax);
  if (atom->memcheck("image"))
    bytes += memory->usage(image, nmax);

  if (atom->memcheck("part"))
    bytes += memory->usage(part, nmax);

  if (atom->memcheck("x"))
    bytes += memory->usage(x, nmax, 3);
  if (atom->memcheck("v"))
    bytes += memory->usage(v, nmax, 3);

  if (atom->memcheck("xdot"))
    bytes += memory->usage(xdot, nmax, 3);

  if (atom->memcheck("vprev"))
    bytes += memory->usage(vprev, nmax, ISPH_BDF_MAX_ORDER, 3);
  if (atom->memcheck("xprev"))
    bytes += memory->usage(xprev, nmax, ISPH_BDF_MAX_ORDER, 3);

  if (atom->memcheck("f"))
    bytes += memory->usage(f, nmax * comm->nthreads, 3);

  if (atom->memcheck("pressure"))
    bytes += memory->usage(pressure, nmax);

  if (atom->memcheck("eps"))
    bytes += memory->usage(eps, nmax);
  if (atom->memcheck("sigma"))
    bytes += memory->usage(sigma, nmax);

  if (atom->memcheck("psi"))
    bytes += memory->usage(psi, nmax);
  if (atom->memcheck("phi"))
    bytes += memory->usage(phi, nmax);

  if (atom->memcheck("concentration"))
    bytes += memory->usage(concentration, ISPH_MAX_CONCENTRATION, nmax);    

  if (atom->memcheck("vfrac"))
    bytes += memory->usage(vfrac, nmax);
  if (atom->memcheck("rmass"))
    bytes += memory->usage(rmass, nmax);
  if (atom->memcheck("viscosity"))
    bytes += memory->usage(viscosity, nmax);
  if (atom->memcheck("density"))
    bytes += memory->usage(density, nmax);

  if (atom->memcheck("molecule")) bytes += memory->usage(molecule,nmax);
  if (atom->memcheck("nspecial")) bytes += memory->usage(nspecial,nmax,3);
  if (atom->memcheck("special"))
    bytes += memory->usage(special,nmax,atom->maxspecial);

  if (atom->memcheck("num_bond")) bytes += memory->usage(num_bond,nmax);
  if (atom->memcheck("bond_type"))
    bytes += memory->usage(bond_type,nmax,atom->bond_per_atom);
  if (atom->memcheck("bond_atom"))
    bytes += memory->usage(bond_atom,nmax,atom->bond_per_atom);

  return bytes;
}

void AtomVecISPH::pack_data(double **file) {
  error->all(FLERR, "AtomVecISPH::pack_data - not implemented");
}

void AtomVecISPH::write_data(FILE *file, int i, double **data) {
  error->all(FLERR, "AtomVecISPH::write_data - not implemented");
}

