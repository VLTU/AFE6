function [A, b] = assemble_fe_problem(fe_space, ref_data, geom_map, problem_B, ...
                        problem_L)


m = size(geom_map.map, 3);
n = fe_space.n;
nq = size(ref_data.evaluation_points, 1);

A = zeros(n);
b = zeros(n, 1);

for l = 1:m

    I = fe_space.support_and_extraction(l).supported_bases;
    for i = 1:size(I, 2)
        
        Ni = fe_space.support_and_extraction(l).extraction_coefficients(i, :) ...
             * ref_data.reference_basis; 
                
        temp1 = fe_space.support_and_extraction(l).extraction_coefficients(i, :) ...
               * ref_data.reference_basis_derivatives(:, :, 1);
        temp2 = fe_space.support_and_extraction(l).extraction_coefficients(i, :) ...
               * ref_data.reference_basis_derivatives(:, :, 2);

        Ni_der_1 = temp1 * geom_map.imap_derivatives(:, 1, l) ...
                   + temp2 * geom_map.imap_derivatives(:, 2, l); 
        Ni_der_2 = temp1 * geom_map.imap_derivatives(:, 3, l) ...
                   + temp2 * geom_map.imap_derivatives(:, 4, l); 
    
        for r = 1:nq
            
            b(I(i)) = b(I(i)) + problem_L(geom_map.map(r,:,l), Ni(r)) ...
                                * ref_data.quadrature_weights(r) ...
                                * det(reshape(geom_map.map_derivatives(r,:,l), 2,2));

        end

        for j = 1:size(I, 2)

            Nj = fe_space.support_and_extraction(l).extraction_coefficients(j, :) ...
                 * ref_data.reference_basis; %Should be [1 X n_q]
                    
            temp1 = fe_space.support_and_extraction(l).extraction_coefficients(j, :) ...
                   * ref_data.reference_basis_derivatives(:, :, 1);
            temp2 = fe_space.support_and_extraction(l).extraction_coefficients(j, :) ...
                   * ref_data.reference_basis_derivatives(:, :, 2);

            Nj_der_1 = temp1 * geom_map.imap_derivatives(:, 1, l) ...
                       + temp2 * geom_map.imap_derivatives(:, 2, l); 
            Nj_der_2 = temp1 * geom_map.imap_derivatives(:, 3, l) ...
                       + temp2 * geom_map.imap_derivatives(:, 4, l); 

            for r = 1:nq    

                A(I(i),I(j)) = A(I(i), I(j)) ...
                               + problem_B(geom_map.map(r,:,l), ...
                                           Ni, ...
                                           Ni_der_1, ...
                                           Ni_der_2, ... 
                                           Nj, ...
                                           Nj_der_1, ...
                                           Nj_der_2) ...
                                * ref_data.quadrature_weights(r) ...
                                * det(reshape(geom_map.map_derivatives(r,:,l), 2,2));

            end
        end 
    end
end

A(fe_space.boundary_bases,:) = [];
A(:,fe_space.boundary_bases) = [];
b(fe_space.boundary_bases) = [];

end
