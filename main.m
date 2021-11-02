model = createpde();
importGeometry(model, 'template-model\files\CORNER.STL');
generateMesh(model);
pdeplot3D(model)

% M = compute_body_inertia_coefficients_matrix()
% C = compute_hydrostatic_restoring_terms_matrix()
% [A, B] = calculate_added_mass_and_damping_coefficient_matrices()
% F = calculate_exciting_forces_vector()
