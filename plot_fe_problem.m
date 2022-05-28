function plot_fe_problem(fe_space, ref_data, geom_map, u)

m = size(geom_map.map, 3);
n = fe_space.n;
nq = size(ref_data.evaluation_points, 1);
x1 = zeros(1, m * nq);
x2 = zeros(1, m * nq);
y = zeros(1, m * nq);
y_der_1 = zeros(1, m * nq);
y_der_2 = zeros(1, m * nq);

for l = 1:m

    I = fe_space.support_and_extraction(l).supported_bases;
    for i = 1:size(I, 1)
        
        Ni = fe_space.support_and_extraction(l).extraction_coefficients(i, :) ...
             * ref_data.reference_basis;
        
        temp1 = fe_space.support_and_extraction(l).extraction_coefficients(i, :) ...
               * ref_data.reference_basis_derivatives(:, :, 1);
        temp2 = fe_space.support_and_extraction(l).extraction_coefficients(i, :) ...
               * ref_data.reference_basis_derivatives(:, :, 2);

        Ni_der_1 = temp1 .* geom_map.imap_derivatives(:, 1, l)' ...
                   + temp2 .* geom_map.imap_derivatives(:, 2, l)'; 
        Ni_der_2 = temp1 .* geom_map.imap_derivatives(:, 3, l)' ...
                   + temp2 .* geom_map.imap_derivatives(:, 4, l)';
        
        x1((l-1) * nq + 1:l * nq) = geom_map.map(:, 1, l);

        x2((l-1) * nq + 1:l * nq) = geom_map.map(:, 2, l);

        y((l-1) * nq + 1:l * nq) = y((l-1) * nq + 1:l * nq) + u(I(i)) .* Ni;

        y_der_1((l-1) * nq + 1:l * nq) = y_der_1((l-1) * nq + 1:l * nq) + u(I(i)) .* Ni_der_1;

        y_der_2((l-1) * nq + 1:l * nq) = y_der_2((l-1) * nq + 1:l * nq) + u(I(i)) .* Ni_der_2;
    end
end

fig1 = figure;
figure(fig1);
scatter(x1, x2, [], y, 'filled');
title('Solution plot'); 

fig2 = figure;
figure(fig2);
scatter(x1, x2, [], y_der_1, 'filled');
title('Derivative first component plot');

fig2 = figure;
figure(fig2);
scatter(x1, x2, [], y_der_2, 'filled');
title('Derivative second component plot');

end
