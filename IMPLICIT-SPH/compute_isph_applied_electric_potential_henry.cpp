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
#include "domain.h"

#include <iostream>
#include "utils.h"

#include "pair_isph_corrected.h"
//#include "pair_isph_mls.h"

#include "compute_isph_applied_electric_potential_henry.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_AppliedElectricPotentialHenry::ComputeISPH_AppliedElectricPotentialHenry"
ComputeISPH_AppliedElectricPotentialHenry::ComputeISPH_AppliedElectricPotentialHenry(LAMMPS *lmp, int narg, char **arg) 
  : Compute(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Compute ipsh/efield-applied/analytic/henry command requires a matching atom_style, i.e., isph");

  peratom_flag = 1;
  size_peratom_cols = 3;

  if (narg < 8)
    error->all(FLERR,"Illegal number of arguments for compute isph command");

  // arg [3] efield applied to x direction
  _eapp = force->numeric(FLERR, arg[3]);

  // arg [4,5,6] centeroid
  _c[0] = force->numeric(FLERR, arg[4]);
  _c[1] = force->numeric(FLERR, arg[5]);
  _c[2] = force->numeric(FLERR, arg[6]);

  // arg [7] radius
  _a = force->numeric(FLERR, arg[7]);

  // arg [8] lambda
  _sratio = force->numeric(FLERR, arg[8]);

  if (comm->me == 0) 
    cout << " - ComputeISPH_AppliedElectricPotentialHenry - " << endl
         << "   applied efield = " << _eapp << endl
         << "   conductivity ratio (sigma_solid/sigma_fluid) = " << _sratio << endl
         << "   radius of sphere = " << _a << endl
         << "   center of sphere = " << _c[0] << ", " << _c[1] << ", " << _c[2] << endl;
  
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_AppliedElectricPotentialHenry::~ComputeISPH_AppliedElectricPotentialHenry"
ComputeISPH_AppliedElectricPotentialHenry::~ComputeISPH_AppliedElectricPotentialHenry() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_AppliedElectricPotentialHenry::init"
void ComputeISPH_AppliedElectricPotentialHenry::init() {
  FUNCT_ENTER(comm->me);

  int count = 0;
  const string str = "isph/applied-electric-potential";

  // one solution for electric potential is allowed
  for (int i=0;i<modify->ncompute;++i)
    count += (strncmp(modify->compute[i]->style, str.c_str(), str.size()) == 0);

  if (count > 1 && comm->me == 0)
    error->all(FLERR,"More than one compute isph/applied-electric-potential");

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_AppliedElectricPotentialHenry::compute_scalar"
double ComputeISPH_AppliedElectricPotentialHenry::compute_scalar() {
  FUNCT_ENTER(comm->me);
    
  const int dim = domain->dimension;
  const int nlocal = atom->nlocal;
  
  int *type = atom->type;
  double **x = atom->x;
  
  double *phi = atom->phi;
  double **phigrad = atom->phigrad;

  double cnt = 0.0, err_phi = 0.0, err_phigrad = 0.0, norm_phi = 0.0, norm_phigrad = 0.0;

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));

  FilterBinary filter;
  filter.setPairYes(PairISPH::Fluid - PairISPH::BufferDirichlet - PairISPH::BufferNeumann);

  for (int i=0;i<nlocal;++i) {
    double phi_at_i = 0.0, phigrad_at_i[3], dx[3] = {};
    for (int k=0;k<dim;++k) 
      dx[k] = x[i][k] - _c[k];

    if (filter.yes(pair->getParticleKind(type[i]))) {
      // compute analytic solutions
      computeHenrySolution(dim, dx, phi_at_i, phigrad_at_i);

      const double diff_phi = pow(phi_at_i - phi[i], 2); 
      const double diff_phigrad[3] = { pow(phigrad[i][0] - phigrad_at_i[0], 2),
                                       pow(phigrad[i][1] - phigrad_at_i[1], 2),
                                       pow(phigrad[i][2] - phigrad_at_i[2], 2) };
      err_phi += diff_phi;
      err_phigrad += (diff_phigrad[0] + diff_phigrad[1] + (dim > 2 ? diff_phigrad[2] : 0.0));
      
      norm_phi += pow(phi_at_i, 2);
      norm_phigrad += (pow(phigrad_at_i[0], 2) + pow(phigrad_at_i[1], 2) + (dim > 2 ? pow(phigrad_at_i[2], 2) : 0.0));

      ++cnt;
    }
  }
  
  double sbuf[5] = {}, rbuf[5] = {}; 
  sbuf[0] = cnt; 
  sbuf[1] = err_phi; 
  sbuf[2] = err_phigrad;
  sbuf[3] = norm_phi; 
  sbuf[4] = norm_phigrad;
  MPI_Reduce(sbuf, rbuf, 5, MPI_DOUBLE, MPI_SUM, 0, world);

  if (comm->me == 0) {
    cnt = rbuf[0];
    err_phi = sqrt(rbuf[1] / cnt);
    err_phigrad = sqrt(rbuf[2] / cnt);
    norm_phi = sqrt(rbuf[3] / cnt);
    norm_phigrad = sqrt(rbuf[4] / cnt);

    auto prec = cout.precision(8);
    
    cout << "ComputeISPH_AppliedElectricPotentialHenry:: # particles  = " << cnt  << endl
         << "                                            phi.error    = " << err_phi << endl
         << "                                            phi.norm     = " << norm_phi << endl
         << "                                            relative err = " << err_phi/norm_phi << endl
         << "                                            gradphi.error= " << err_phigrad << endl
         << "                                            gradphi.norm = " << norm_phigrad << endl
         << "                                            relative err = " << err_phigrad/norm_phigrad << endl;

    cout.unsetf(ios::scientific);
    cout.precision(prec);
  }

  FUNCT_EXIT(comm->me);  
  
  return 0.0;
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_AppliedElectricPotentialHenry::compute_peratom"
void ComputeISPH_AppliedElectricPotentialHenry::compute_peratom() {
  FUNCT_ENTER(comm->me);

  if (invoked_peratom != update->ntimestep) {
    invoked_peratom = update->ntimestep;
    
    const int dim = domain->dimension;
    const int nsize = (atom->nlocal + atom->nghost);

    double **x = atom->x;

    double *phi = atom->phi;
    double **phigrad = atom->phigrad;

    array_atom = phigrad;

    for (int i=0;i<nsize;++i) {
      double dx[3] = {};
      for (int k=0;k<dim;++k) 
        dx[k] = x[i][k] - _c[k];

      computeHenrySolution(dim, dx, phi[i], phigrad[i]);
    }
  }

  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_AppliedElectricPotentialHenry::computeHenrySolution"
void ComputeISPH_AppliedElectricPotentialHenry::computeHenrySolution(const int dim, const double *dx,
                                                                     double &phi_at_i, double *phigrad_at_i) {
  //FUNCT_ENTER(comm->me);

  double lambda = 0.0;
  if (dim > 2)
    lambda = (1.0 - _sratio)/(2.0 + _sratio);
  else
    lambda = (1.0 - _sratio)/(1.0 + _sratio);

  const double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

  // compute analytic solutions
  if (r < _a) {
    phigrad_at_i[0] = -_eapp*(1.0 + lambda);
    phigrad_at_i[1] = 0.0;
    phigrad_at_i[2] = 0.0;
    phi_at_i = -_eapp*(1.0 + lambda)*dx[0];
  } else {
    if (dim > 2) {
      const double a3 = pow(_a,3);
      const double r5 = pow(r,5);
      phigrad_at_i[0] = _eapp*(-1.0 + a3*lambda*(2*dx[0]*dx[0] - dx[1]*dx[1] - dx[2]*dx[2])/r5);
      phigrad_at_i[1] = 3*a3*_eapp*lambda*dx[0]*dx[1]/r5;
      phigrad_at_i[2] = 3*a3*_eapp*lambda*dx[0]*dx[2]/r5;
      phi_at_i = -_eapp*(1.0 + lambda*pow(_a/r,3))*dx[0];
    } else {
      const double a2 = pow(_a,2);
      const double r4 = pow(r,4);
      phigrad_at_i[0] = _eapp*(-1.0 + a2*lambda*(dx[0]*dx[0] - dx[1]*dx[1])/r4);
      phigrad_at_i[1] = 2*a2*_eapp*lambda*dx[0]*dx[1]/r4;
      phigrad_at_i[2] = 0.0;
      phi_at_i = -_eapp*(1.0 + lambda*pow(_a/r,2))*dx[0];
    }
  }

  //FUNCT_EXIT(comm->me);  
}

