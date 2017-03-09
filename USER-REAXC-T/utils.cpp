#include "utils.h"
#include "utils_reference.h"

#include "Teuchos_ParameterList.hpp"

// global declaration
namespace LAMMPS_NS {

  int g_funct_counter = 0;
  Utils util;

  using namespace Teuchos;

  ParameterList g_params;

  RCP<Teuchos::Time> g_timer_init = TimeMonitor::getNewTimer("ISPH: init"),
    g_timer_compute_volumes = TimeMonitor::getNewTimer("ISPH: computeVolumes"),
    g_timer_compute_gradient_correction = TimeMonitor::getNewTimer("ISPH: computeGradientCorrection"),
    g_timer_compute_laplacian_correction = TimeMonitor::getNewTimer("ISPH: computeLaplacianCorrection"),
    g_timer_compute_normals = TimeMonitor::getNewTimer("ISPH: computeNormals"),
    g_timer_compute_curlcurl = TimeMonitor::getNewTimer("ISPH: computeCurlCurl"),
    g_timer_compute_f_poisson_boltzmann = TimeMonitor::getNewTimer("ISPH: computeF"),
    g_timer_compute_jacobian_poisson_boltzmann = TimeMonitor::getNewTimer("ISPH: computeJacobian"),
    g_timer_compute_poisson_boltzmann = TimeMonitor::getNewTimer("ISPH: computePoissonBoltzmann"),
    g_timer_compute_applied_electric_field = TimeMonitor::getNewTimer("ISPH: computeAppliedElectricField"),
    g_timer_compute_applied_electric_potential = TimeMonitor::getNewTimer("ISPH: computeAppliedElectricPotential"),
    g_timer_solve_applied_electric_potential = TimeMonitor::getNewTimer("ISPH: solveAppliedElectricPotential"),
    g_timer_compute_surface_tension = TimeMonitor::getNewTimer("ISPH: computeSurfaceTension"),
    g_timer_compute_solute_transport = TimeMonitor::getNewTimer("ISPH: computeSoluteTransport"),
    g_timer_compute_solute_transport_species = TimeMonitor::getNewTimer("ISPH: computeSoluteTransportSpecies"),
    g_timer_solve_solute_transport_species = TimeMonitor::getNewTimer("ISPH: solveSoluteTransportSpecies"),
    g_timer_correct_velocity = TimeMonitor::getNewTimer("ISPH: correctVelocity"),
    g_timer_correct_pressure = TimeMonitor::getNewTimer("ISPH: correctPressure"),
    g_timer_compute_helmholtz = TimeMonitor::getNewTimer("ISPH: computeHelmholtz"),
    g_timer_compute_block_helmholtz = TimeMonitor::getNewTimer("ISPH: computeBlockHelmholtz"),
    g_timer_solve_helmholtz = TimeMonitor::getNewTimer("ISPH: solveHelmholtz"),
    g_timer_compute_poisson = TimeMonitor::getNewTimer("ISPH: computePoisson"),
    g_timer_solve_poisson = TimeMonitor::getNewTimer("ISPH: solvePoisson"),
    g_timer_compute_imcompressible_navier_stokes = TimeMonitor::getNewTimer("ISPH: computeImcompressibleNavierStokes"),
    g_timer_compute_mass_matrix = TimeMonitor::getNewTimer("ISPH: computeMassMatrix"),
    g_timer_compute_mass_matrix_compact_poisson = TimeMonitor::getNewTimer("ISPH: computeMassMatrixCompactPoisson"),
    g_timer_ale_predict_velocity = TimeMonitor::getNewTimer("ISPH: predictAleVelocity"),
    g_timer_ale_correct_velocity = TimeMonitor::getNewTimer("ISPH: correctAleVelocity");
}
