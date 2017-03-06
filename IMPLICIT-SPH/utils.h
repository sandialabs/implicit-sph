#ifndef __UTILS_H__
#define __UTILS_H__

#include "macrodef.h"
#include "utils_reference.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"  
#include "Teuchos_TimeMonitor.hpp"

namespace LAMMPS_NS {

  extern int g_funct_counter;
  extern class Teuchos::ParameterList g_params;

  typedef class UtilsReference Utils;
  extern Utils util;

  extern class Teuchos::RCP<Teuchos::Time> g_timer_init, 
    g_timer_compute_volumes, 
    g_timer_compute_gradient_correction, 
    g_timer_compute_laplacian_correction,
    g_timer_compute_normals,
    g_timer_compute_curlcurl,
    g_timer_compute_f_poisson_boltzmann,
    g_timer_compute_jacobian_poisson_boltzmann,
    g_timer_compute_poisson_boltzmann,
    g_timer_compute_applied_electric_field,
    g_timer_compute_applied_electric_potential,
    g_timer_solve_applied_electric_potential,
    g_timer_compute_surface_tension,
    g_timer_compute_solute_transport,
    g_timer_compute_solute_transport_species,
    g_timer_solve_solute_transport_species,
    g_timer_compute_helmholtz,
    g_timer_compute_block_helmholtz,
    g_timer_solve_helmholtz,
    g_timer_compute_poisson,
    g_timer_solve_poisson,
    g_timer_correct_velocity,
    g_timer_correct_pressure,
    g_timer_compute_imcompressible_navier_stokes,
    g_timer_compute_mass_matrix,
    g_timer_compute_mass_matrix_compact_poisson,
    g_timer_ale_predict_velocity,
    g_timer_ale_correct_velocity;

  // need class for supported type

  inline bool is_same     (int itype, int jtype)  { return (itype == jtype); }
  inline bool is_different(int itype, int jtype)  { return (itype != jtype); }
  inline bool is_any      (int itype, int jtype)  { return true; }

}

#endif
