fe_file = load("geometry/distressed_robotD.mat");
nq = 9; p = [2,2]; k = 1;
f = @(x) 1;

problem_B = @(x, u, u_der_1, u_der_2, w, w_der_1, w_der_2) ...
            w_der_1 * k * u_der_1 + w_der_2 * k * u_der_1;
problem_L = @(x, w, w_der_1, w_der_2) w * f(x);

ref_data = create_ref_data(nq, p, 'integrate');
geom_map = create_geometric_map(fe_file.fe_geometry, ref_data);
fe_space = fe_file.fe_space;

[A, b] = assemble_fe_problem(fe_space, ref_data, ...
            geom_map, problem_B, problem_L);
