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

#include "Teuchos_XMLParameterListCoreHelpers.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

#include "pair_isph_corrected.h"
//#include "pair_isph_mls.h"

#include "fix_isph_error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define DASHLINE "==================================================="

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Error::FixISPH_Error"
FixISPH_Error::FixISPH_Error(LAMMPS *lmp, int narg, char **arg) 
  : Fix(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1) 
    error->all(FLERR, "Fix ipsh/error command requires a matching atom_style, i.e., isph");
  if (narg < 3)
    error->all(FLERR,"Illegal number of arguments for fix isph command");

  if (g_params.isSublist("Analytic Solution")) {
    _param = g_params.sublist("Analytic Solution");
    _is_solution_known = 1;
  } else {
    _is_solution_known = 0;
  }
  time_integrate = 0;

  FUNCT_EXIT(comm->me);
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Error::~FixISPH_Error"
FixISPH_Error::~FixISPH_Error() {
  FUNCT_ENTER(comm->me);
  FUNCT_EXIT(comm->me);  
}
/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Error::setmask"
int FixISPH_Error::setmask() {
  FUNCT_ENTER(comm->me);

  int mask = 0;
  mask |= FINAL_INTEGRATE;

  FUNCT_EXIT(comm->me);
  return mask;
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Error::init"
void FixISPH_Error::init() {
  FUNCT_ENTER(comm->me);

  if (!_is_solution_known) {
    FUNCT_EXIT(comm->me);
    return;
  }
  
  _exact = Teuchos::rcp(new PG_RuntimeCompiler::Function(0, "Analytic Solution"));

  //if (comm->me == 0)
  //  cout << g_params << endl ;

  if (!_param.isSublist("Variable List") ||
      !_param.isSublist("Function List")) {
    if (comm->me == 0) {
      cout << "FixISPH_Error:: Variable and/or function lists are empty" << endl;      
      cout << "                error is not computed" << endl;
    }
    return;
  }
    
  // load variables
  _vars = _param.sublist("Variable List");   
  if (comm->me == 0) {
    cout << endl << "Variable List";
    cout << endl << "=============" << endl;
    _vars.print(cout);
    cout << endl;
  }
    
  long long int cnt = 0;
  for (auto iter=_vars.begin();iter!=_vars.end();++iter) {
    string name = _vars.name(iter);
    string type = _vars.get(name, "");
    _exact->addVar(type, name);
    _vars.set(name, to_string(cnt));

    // initialize all variables zero; necessary to work with RTC
    _exact->varValueFill(stoi(_vars.get(name,"")), 0);     
    ++cnt;
  }  

  // load funciton list
  _func = _param.sublist("Function List");
  if (comm->me == 0) {
    cout << endl << "Function List";
    cout << endl << "=============" << endl;
    _func.print(cout);
    cout << endl;
  }

  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Error::final_integrate"
void FixISPH_Error::final_integrate() {
  FUNCT_ENTER(comm->me);

  if (g_params.isSublist("Physics Configuration")) {
    auto& phys_param = g_params.sublist("Physics Configuration");

    if (phys_param.get("Poisson Boltzmann","Disabled") == "Enabled") {
      computePoissonBoltzmannError();
    } else {
      if (comm->me == 0) {
        cout << "FixISPH_Error:: Poisson Boltzmann is disabled" << endl;      
        cout << "                error is not computed" << endl;
        cout << DASHLINE << endl << endl;
      }
    }
    if (phys_param.get("Incompressible Navier Stokes","Disabled") == "Enabled") {
      computeIncompressibleNavierStokesError();
    } else {
      if (comm->me == 0) {
        cout << "FixISPH_Error:: Incomppressible Navier Stokes is disabled" << endl;      
        cout << "                error is not computed" << endl;
        cout << DASHLINE << endl << endl;
      }
    }

  } else {
    if (comm->me == 0) {
      cout << "FixISPH_Error:: Physics Configuration is empty in the top-level parameter list" << endl;      
      cout << "                error is not computed" << endl;
      cout << DASHLINE << endl << endl;
    }
  }
  FUNCT_EXIT(comm->me);
}

#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Error::computePoissonBoltzmannError"
void FixISPH_Error::computePoissonBoltzmannError() {
  FUNCT_ENTER(comm->me);

  int nlocal = atom->nlocal;
  int *type = atom->type;
  double **x = atom->x;
  double *vfrac = atom->vfrac;

  double *eps = atom->eps;
  double *psi = atom->psi;
  double **psigrad = atom->psigrad;
  
  bool func_psi_exist    = (_func.get("psi","none") != "none");
  bool func_grad_x_exist = (_func.get("psi.grad.x","none") != "none");
  bool func_grad_y_exist = (_func.get("psi.grad.y","none") != "none");
  bool func_grad_z_exist = (_func.get("psi.grad.z","none") != "none");

  bool func_grad_exist   = (func_grad_x_exist | func_grad_y_exist | func_grad_z_exist);
  bool func_exist        = (func_grad_exist | func_psi_exist);
  
  Epetra_Map mappy(-1, nlocal, 1, Epetra_MpiComm(world));
  Epetra_Vector 
    err_psi(mappy), sol_psi(mappy), 
    err_grad(mappy), sol_grad(mappy);

  string idx;
  // fill the optional variables
  if (g_params.isSublist("Poisson Boltzmann")) {
    auto& pb_param = g_params.sublist("Poisson Boltzmann");
    if ((idx = _vars.get("psiref","none")) != "none")
      _exact->varValueFill(stoi(idx), pb_param.get("psiref", 1.0));
    if ((idx = _vars.get("ezcb","none")) != "none")
      _exact->varValueFill(stoi(idx), pb_param.get("ezcb", 1.0));
  } else {
    if ((idx = _vars.get("psiref","none")) != "none")
      _exact->varValueFill(stoi(idx), 1.0);
    if ((idx = _vars.get("ezcb","none")) != "none")
      _exact->varValueFill(stoi(idx), 1.0);
  }

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));  
  if (pair == NULL) 
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");

  int ntotal = 0;
  double volume = 0.0;
  for (int i=0;i<nlocal;++i) {
    
    if (pair->getParticleKind(type[i]) == PairISPH_Corrected::Solid)
      continue;

    // fill the essential variables
    _exact->varValueFill(stoi(_vars.get("pt.x","")), x[i][0]);
    _exact->varValueFill(stoi(_vars.get("pt.y","")), x[i][1]);
    _exact->varValueFill(stoi(_vars.get("pt.z","")), x[i][2]);
    _exact->varValueFill(stoi(_vars.get("eps","")),  eps[i]);
    
    if (func_psi_exist) {
      _exact->addBody(_func.get("psi",""));
      _exact->execute();
      
      double val  = _exact->getValueOfVar("psi");
      double diff = abs(psi[i] - val);

      err_psi[i] += diff*diff;
      sol_psi[i] += val*val;
    }
    if (func_grad_x_exist) {
      _exact->addBody(_func.get("psi.grad.x",""));
      _exact->execute();
      
      double val  = _exact->getValueOfVar("psi.grad.x");
      double diff = abs(psigrad[i][0] - val);
      
      err_grad[i] += diff*diff;
      sol_grad[i] += val*val;
    }
    if (func_grad_y_exist) {
      _exact->addBody(_func.get("psi.grad.y",""));
      _exact->execute();
      
      double val  = _exact->getValueOfVar("psi.grad.y");
      double diff = abs(psigrad[i][1] - val);
      
      err_grad[i] += diff*diff;
      sol_grad[i] += val*val;
    }
    if (func_grad_z_exist) {
      _exact->addBody(_func.get("psi.grad.x",""));
      _exact->execute();
      
      double val  = _exact->getValueOfVar("psi.grad.z");
      double diff = abs(psigrad[i][2] - val);
      
      err_grad[i] += diff*diff;
      sol_grad[i] += val*val;
    }

    // consider volume fraction occupied by the particle
    switch (pair->getDiscretizationInfo()) {
    case PairISPH::Corrected:
      volume += vfrac[i];
      break;
    case PairISPH::MLS:
      volume += 1.0;
      break;
    }
    ++ntotal;
  }

  double volume_all = 0.0;
  int ntotal_all = 0;
  MPI_Allreduce(&volume, &volume_all, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&ntotal, &ntotal_all, 1, MPI_INT, MPI_SUM, world);

  double err_sum_psi_sq = 0.0; err_psi.Norm1(&err_sum_psi_sq);
  double sol_sum_psi_sq = 0.0; sol_psi.Norm1(&sol_sum_psi_sq);

  double err_sum_grad_sq = 0.0; err_grad.Norm1(&err_sum_grad_sq);
  double sol_sum_grad_sq = 0.0; sol_grad.Norm1(&sol_sum_grad_sq);

  // let's compare apple to apple
  double denominator(ntotal_all); // denominator(volume_all);
  
  double err_norm_psi = sqrt(err_sum_psi_sq/denominator);
  double sol_norm_psi = sqrt(sol_sum_psi_sq/denominator);

  double err_norm_grad = sqrt(err_sum_grad_sq/denominator);
  double sol_norm_grad = sqrt(sol_sum_grad_sq/denominator);
  
  if (comm->me == 0) {
    streamsize prec = cout.precision();
    cout.precision(15);
    cout << scientific;
    if (func_exist) {
      cout << "FixISPH_Error:: analytic solution is known for Poisson Boltzmann" << endl;
      cout << "                error is computed against the analytic solution" << endl;
      cout << "                total # of particles = " << ntotal_all << endl;
      cout << "                total volume         = " << volume_all << endl;
      cout << endl;
      cout << "                note: if MLS is used total volume is same as the total number of particles" << endl;
      if (func_psi_exist) {
        cout << endl;
        cout << "                sol.psi.norm2        = " << sol_norm_psi << endl;
        cout << "                err.psi.norm2        = " << err_norm_psi << endl;
        cout << "                relative error       = " << err_norm_psi/sol_norm_psi << endl;
      }
      if (func_grad_exist) {
        cout << endl;
        cout << "                sol.psi.grad.norm2   = " << sol_norm_grad << endl;
        cout << "                err.psi.grad.norm2   = " << err_norm_grad << endl;
        cout << "                relative error       = " << err_norm_grad/sol_norm_grad << endl;
      }
      cout << DASHLINE << endl << endl;
    } else {
      cout << "FixISPH_Error:: analytic solution for Poisson Boltzmann is unknown" << endl;      
      cout << "                error is not computed" << endl;
      cout << DASHLINE << endl << endl;
    }
    cout.unsetf(ios::scientific);
    cout.precision(prec);
  }
  FUNCT_EXIT(comm->me);
}

#undef  __FUNCT__
#define __FUNCT__ "FixISPH_Error::computeIncompressibleNavierStokesError"
void FixISPH_Error::computeIncompressibleNavierStokesError() {
  FUNCT_ENTER(comm->me);

  PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
  if (pair == NULL) 
    error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH");

  int nlocal = atom->nlocal;
  int *type = atom->type;
  double **x = atom->x;
  double *vfrac = atom->vfrac;

  double *p = atom->pressure;
  double *nu = atom->viscosity;
  double *rho = atom->density;
  double **v;
  if (pair->isAleEnabled())
    v = atom->v;
  else
    v = pair->getVstar();

  // we compute 
  // error = | u_exact{x^n, t^{n+1}} - u_approx{x^n, t^{n+1}}

  bool func_u_x_exist    = (_func.get("u.x","none") != "none");
  bool func_u_y_exist    = (_func.get("u.y","none") != "none");
  bool func_u_z_exist    = (_func.get("u.z","none") != "none");
  bool func_p_exist      = (_func.get("p","none") != "none");
  
  bool func_u_exist      = (func_u_x_exist | func_u_y_exist | func_u_z_exist);
  bool func_exist        = (func_u_exist | func_p_exist);

  Epetra_Map mappy(-1, nlocal, 1, Epetra_MpiComm(world));
  Epetra_Vector 
    err_u(mappy), sol_u(mappy), 
    err_p(mappy), sol_p(mappy);
  
  if (g_params.isSublist("Incompressible Navier Stokes")) {
    auto& ns_param = g_params.sublist("Incompressible Navier Stokes");
    if (_vars.get("umax","none") != "none")
      _exact->varValueFill(stoi(_vars.get("umax","")), ns_param.get("umax", 0.1));
    if (_vars.get("g.x","none") != "none")
      _exact->varValueFill(stoi(_vars.get("g.x","")), ns_param.get("g.x", 0.0));
    if (_vars.get("g.y","none") != "none")
      _exact->varValueFill(stoi(_vars.get("g.y","")), ns_param.get("g.y", 0.0));
    if (_vars.get("g.z","none") != "none")
      _exact->varValueFill(stoi(_vars.get("g.z","")), ns_param.get("g.z", 0.0));
  } else {
    if (_vars.get("umax","none") != "none")
      _exact->varValueFill(stoi(_vars.get("umax","")), 0.1);
    if (_vars.get("g.x","none") != "none")
      _exact->varValueFill(stoi(_vars.get("g.x","")), 0.0);
    if (_vars.get("g.y","none") != "none")
      _exact->varValueFill(stoi(_vars.get("g.y","")), 0.0);
    if (_vars.get("g.z","none") != "none")
      _exact->varValueFill(stoi(_vars.get("g.z","")), 0.0);
  }

  double mean_p(0.);
  mean_p = pair->computeZeroMeanPressure(p, false);

  int ntotal = 0;
  double volume = 0.0;
  for (int i=0;i<nlocal;++i) {

    if (pair->getParticleKind(type[i]) == PairISPH_Corrected::Solid)
      continue;

    // fill the essential variables
    _exact->varValueFill(stoi(_vars.get("pt.x","")), x[i][0]);
    _exact->varValueFill(stoi(_vars.get("pt.y","")), x[i][1]);
    _exact->varValueFill(stoi(_vars.get("pt.z","")), x[i][2]);
    
    // fill the optional variables
    if (_vars.get("nu","none") != "none")
      _exact->varValueFill(stoi(_vars.get("nu","")),   nu[i]);

    if (_vars.get("rho","none") != "none")
      _exact->varValueFill(stoi(_vars.get("rho","")),  rho[i]);

    if (_vars.get("t","none") != "none")
      _exact->varValueFill(stoi(_vars.get("t","")),    update->dt*update->ntimestep);

    if (func_u_x_exist) {
      _exact->addBody(_func.get("u.x",""));
      _exact->execute();
          
      double val  = _exact->getValueOfVar("u.x");
      double diff = abs(v[i][0] - val);
      
      err_u[i] += diff*diff;
      sol_u[i] += val*val;
    }
    if (func_u_y_exist) {
      _exact->addBody(_func.get("u.y",""));
      _exact->execute();
          
      double val  = _exact->getValueOfVar("u.y");
      double diff = abs(v[i][1] - val);
          
      err_u[i] += diff*diff;
      sol_u[i] += val*val;
    }
    if (func_u_z_exist) {
      _exact->addBody(_func.get("u.z",""));
      _exact->execute();
      
      double val  = _exact->getValueOfVar("u.z");
      double diff = abs(v[i][2] - val);
      
      err_u[i] += diff*diff;
      sol_u[i] += val*val;
    }
    if (func_p_exist) {
      _exact->addBody(_func.get("p",""));
      _exact->execute();
      
      double val  = _exact->getValueOfVar("p");
      double diff = abs(p[i] - mean_p - val);
      
      err_p[i] += diff*diff;
      sol_p[i] += val*val;
    }

    // consider volume fraction occupied by the particle
    switch (pair->getDiscretizationInfo()) {
    case PairISPH::Corrected:
      volume += vfrac[i];
      break;
    case PairISPH::MLS:
      volume += 1.0;
      break;
    }
    ++ntotal;
  }

  double volume_all = 0.0;
  int ntotal_all = 0;

  MPI_Allreduce(&volume, &volume_all, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&ntotal, &ntotal_all, 1, MPI_INT, MPI_SUM, world);

  double err_sum_u_sq = 0.0; err_u.Norm1(&err_sum_u_sq);
  double sol_sum_u_sq = 0.0; sol_u.Norm1(&sol_sum_u_sq);

  double err_sum_p_sq = 0.0; err_p.Norm1(&err_sum_p_sq);
  double sol_sum_p_sq = 0.0; sol_p.Norm1(&sol_sum_p_sq);

  // let's compare apple to apple
  double denominator(ntotal_all); // denominator(volume_all);
      
  double err_norm_u = sqrt(err_sum_u_sq/denominator);
  double sol_norm_u = sqrt(sol_sum_u_sq/denominator);

  double err_norm_p = sqrt(err_sum_p_sq/denominator);
  double sol_norm_p = sqrt(sol_sum_p_sq/denominator);
      

  if (comm->me == 0) {
    streamsize prec = cout.precision();
    cout.precision(15);
    cout << scientific;
    if (func_exist) {
      cout << "FixISPH_Error:: analytic solution is known for Incompressible Navier Stokes" << endl;
      cout << "                error is computed against the analytic solution" << endl;
      cout << "                total # of particles = " << ntotal_all << endl;
      cout << "                total volume         = " << volume_all << endl;
      cout << endl;
      cout << "                note: if MLS is used total volume is same as the total number of particles" << endl;
      if (func_u_exist) {
        cout << endl;
        cout << "                sol.u.norm2          = " << sol_norm_u << endl;
        cout << "                err.u.norm2          = " << err_norm_u << endl;
        cout << "                relative error       = " << err_norm_u/sol_norm_u << endl;
      }
      if (func_p_exist) {
        cout << endl;
        cout << "                sol.p.norm2          = " << sol_norm_p << endl;
        cout << "                err.p.norm2          = " << err_norm_p << endl;
        cout << "                relative error       = " << err_norm_p/sol_norm_p << endl;
      }
      cout << DASHLINE << endl << endl;
    } else {
      cout << "FixISPH_Error:: analytic solution for Incompressible Navier Stokes is unknown" << endl;      
      cout << "                error is not computed" << endl;
      cout << DASHLINE << endl << endl;
    }
    cout.unsetf(ios::scientific);
    cout.precision(prec);
  }
  FUNCT_EXIT(comm->me);
} 
