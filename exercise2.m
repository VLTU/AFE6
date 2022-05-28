%fe_file = load("geometry/distressed_robotD.mat");
fe_file = load("geometry/distressed_robotDN.mat");
nq = 8; p = [2,2]; k = 1;
%f = @(x) 1;
f = @(x) sin(x(1)) * cos(x(2));

problem_B = @(x, u, u_der_1, u_der_2, w, w_der_1, w_der_2) ...
            w_der_1 * k * u_der_1 + w_der_2 * k * u_der_2;
problem_L = @(x, w, w_der_1, w_der_2) w * f(x);

ref_data = create_ref_data(nq, p, 'integrate');
geom_map = create_geometric_map(fe_file.fe_geometry, ref_data);
fe_space = fe_file.fe_space;

[A, b] = assemble_fe_problem(fe_space, ref_data, ...
                 geom_map, problem_B, problem_L);

temp = linsolve(A, b);
n = fe_space.n;
u = zeros(n, 1);
count = 1;

for i = 1:n
    if ~ismember(i, fe_space.boundary_bases)
        u(i) = temp(count);
        count = count + 1;
    end
end

ref_data = create_ref_data(nq, p, 'plot');
plot_fe_problem(fe_space, ref_data, geom_map, u);
