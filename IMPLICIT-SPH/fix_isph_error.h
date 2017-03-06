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

#ifdef FIX_CLASS

FixStyle(isph/error,FixISPH_Error)

#else

#ifndef LMP_FIX_ISPH_ERROR_H
#define LMP_FIX_ISPH_ERROR_H

#include "fix.h"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "RTC_FunctionRTC.hh"

namespace LAMMPS_NS {

  class FixISPH_Error: public Fix {
    
  public:
    FixISPH_Error(class LAMMPS *, int, char **);
    virtual~FixISPH_Error();

    virtual int setmask();
    virtual void init();
    virtual void final_integrate();

    void computeIncompressibleNavierStokesError();
    void computePoissonBoltzmannError();

  protected:

    Teuchos::ParameterList _param, _vars, _func;
    int _is_solution_known;
    Teuchos::RCP<PG_RuntimeCompiler::Function> _exact;

  };

}

#endif
#endif
