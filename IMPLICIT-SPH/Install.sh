# Install/unInstall package files in LAMMPS

myfiles="
    atom.h \
    dump_custom.h \
    set.h \
    \
    bond_isph.h \
    bond_isph_fene.h \
    bond_isph_fene_expand.h \
    bond_isph_harmonic.h \
    \
    atom_vec_isph.h \
    \
    color.h \
    kernel.h \
    kernel_cubic.h \
    kernel_quintic.h \
    kernel_wendland.h \
    pairwise_force.h \
    scaled_taylor_monomial.h \
    \
    pair_isph.h \
    pair_isph_corrected.h \
    \
    compute_isph_applied_electric_potential_henry.h \
    compute_isph_concentration.h \
    compute_isph_cylinder_porous.h \
    compute_isph_status.h \
    compute_isph_sphere_porous.h \
    compute_isph_velocity_curl.h \
    compute_isph_velocity_divergence.h \
    \
    fix_isph.h \
    fix_isph_error.h \
    fix_isph_shift.h \
    fix_isph_tgv.h \
    fix_isph_modify_type.h \
    fix_isph_modify_concentration.h \
    fix_isph_modify_velocity.h \
    fix_isph_modify_phi.h \
    fix_isph_ignore_phase_gradient.h \
    fix_isph_quit.h \
    \
    macrodef.h \
    utils.h \
    utils_reference.h \
    \
    time_bdf.h \
    \
    precond.h \
    precond_ifpack.h \
    precond_ml.h \
    \
    solver_nox.h \
    solver_nox_impl.h \
    solver_nox_aztecOO.h \
    solver_nox_stratimikos.h \
    solver_lin.h \
    solver_lin_belos.h \
    \
    pair_for.h \
    filter.h \
    functor.h \
    functor_normalize_vector.h \
    functor_exact_solution.h \
    \
    functor_graph.h \
    functor_graph_boundary.h \
    functor_volume.h \
    functor_normal.h \
    functor_gradient_correction.h \
    functor_laplacian_correction.h \
    \
    functor_smooth_field.h \
    functor_gradient.h \
    functor_gradient_operator.h \
    functor_gradient_dot_operator_matrix.h \
    functor_divergence.h \
    functor_laplacian.h \
    functor_laplacian_matrix.h \
    \
    functor_curl.h \
    functor_curlcurl.h \
    \
    functor_phase_gradient.h \
    functor_phase_divergence.h \
    functor_phase_divergence_adami.h \
    functor_correct_phase_normal.h \
    \
    functor_boundary_morris_holmes.h \
    functor_boundary_morris_normal.h \
    functor_boundary_dirichlet.h \
    functor_boundary_navier_slip.h \
    functor_boundary_neumann.h \
    \
    functor_poisson_boltzmann_f.h \
    functor_poisson_boltzmann_extra_f.h \
    functor_poisson_boltzmann_jacobian.h \
    functor_electrostatic_force.h \
    \
    functor_continuum_surface_force.h \
    functor_pairwise_force.h \
    \
    functor_solute_transport.h \
    functor_applied_electric_potential.h \
    \
    functor_random_stress.h \
    \
    functor_incomp_navier_stokes_helmholtz.h \
    functor_incomp_navier_stokes_block_helmholtz.h \
    functor_incomp_navier_stokes_poisson.h \
    functor_correct_pressure.h \
    functor_correct_velocity.h \
    functor_advance_time_begin.h \
    functor_advance_time_end.h \
    functor_compute_shift.h \
    functor_apply_shift.h \
    \
    functor_traction_vector.h \
    \
    mirror.h \
    mirror_morris_holmes.h \
    mirror_morris_normal.h \
    \
    atom.cpp \
    dump_custom.cpp \
    set.cpp \
    \
    bond_isph_fene.cpp \
    bond_isph_fene_expand.cpp \
    bond_isph_harmonic.cpp \
    \
    atom_vec_isph.cpp \
    \
    pair_isph.cpp \
    pair_isph_corrected.cpp \
    \
    compute_isph_applied_electric_potential_henry.cpp \
    compute_isph_concentration.cpp \
    compute_isph_cylinder_porous.cpp \
    compute_isph_status.cpp \
    compute_isph_sphere_porous.cpp \
    compute_isph_velocity_curl.cpp \
    compute_isph_velocity_divergence.cpp \
    \
    fix_isph.cpp \
    fix_isph_error.cpp \
    fix_isph_shift.cpp \
    fix_isph_tgv.cpp \
    fix_isph_modify_type.cpp \
    fix_isph_modify_concentration.cpp \
    fix_isph_modify_velocity.cpp \
    fix_isph_modify_phi.cpp \
    fix_isph_ignore_phase_gradient.cpp \
    fix_isph_quit.cpp \
    \
    utils.cpp \
    utils_reference.cpp \
    \
    solver_lin.cpp 
"

if (test $1 = 1) then
    for file in $myfiles; do cp -p $file .. ; done
elif (test $1 = 0) then
    for file in $myfiles; do rm -f ../$file ; done
fi

# #Hack for using old uncorrected scheme
# cp functor_uncorrected_gradient.h ../functor_gradient.h
# cp functor_uncorrected_divergence.h ../functor_divergence.h
# cp functor_uncorrected_laplacian.h ../functor_laplacian.h
# cp functor_uncorrected_laplacian_matrix.h ../functor_laplacian_matrix.h


    # pair_isph_mls.cpp \
    # kernel_mls.h \
    # pair_isph_mls.h \

    # functor_mls_mass_matrix.h \
    # functor_mls_mass_matrix_compact_poisson.h \
    # functor_mls_normal.h \
    # functor_mls_helper.h \
    # functor_mls_helper_compact_poisson.h \
    # \
    # functor_mls_gradient.h \
    # functor_mls_gradient_operator.h \
    # functor_mls_gradient_compact_poisson.h \
    # functor_mls_divergence.h \
    # functor_mls_laplacian.h \
    # functor_mls_laplacian_compact_poisson.h \
    # functor_mls_laplacian_matrix.h \
    # functor_mls_laplacian_matrix_compact_poisson.h \
    # functor_mls_curl.h \
    # functor_mls_curlcurl.h \
    # \
    # functor_mls_boundary_neumann.h \
    # \
    # test_mls.h \
    # test_mls.cpp \
    # test_mls_gradient_compact_poisson.cpp \
    # test_mls_laplacian_compact_poisson.cpp \
    # test_mls_laplacian_matrix_compact_poisson.cpp
    # compute_isph_status_flow_past_cylinder.h \
    # compute_isph_status_flow_past_cylinder.cpp \
    # \
    # functor_ale_advection.h \
    # functor_ale_advection_matrix.h \
    # functor_ale_predict_velocity.h \
    # functor_ale_predict_velocity_curlcurl.h \
    # functor_ale_correct_velocity.h \
    # functor_ale_incomp_navier_stokes_poisson.h \
    # functor_ale_incomp_navier_stokes_poisson_boundary.h \
    # functor_ale_incomp_navier_stokes_compact_poisson_boundary.h \
    # functor_ale_incomp_navier_stokes_helmholtz.h \
    # functor_ale_incomp_navier_stokes_helmholtz_curlcurl.h \
    # functor_ale_track_particles.h \
    # functor_ale_apply_shift.h \
