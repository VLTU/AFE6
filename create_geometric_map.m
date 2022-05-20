function geometric_map = create_geometric_map(fe_geometry, ref_data)
%Write about i/o

%map: (nq x 2 x m): (i,:,k)-th row: phi_{k}(xi^1_i,xi^2_i)
Xij = fe_geometry.map_coefficients;
%Generating map by matrix product. 

map = zeros(size(ref_data.evaluation_points,1),2,fe_geometry.m);
map_derivatives = zeros(size(ref_data.evaluation_points,1),4,fe_geometry.m);
%Generating map 
for k = 1:fe_geometry.m
    temp = (Xij(:,:,k).' * ref_data.reference_basis);
    map(:,:,k) = temp.';
end 

%Generating map derivatives 
%map derivative row j contains [dphi_1/dxi_1,dphi_2/dxi_1,dphi_1/dxi_2,dphi_2/dxi2]
%evaluated at row corresponding to evaluation point j. 
for k = 1:fe_geometry.m
    %Generating derivative by first coordinate 
    temp_1 = (Xij(:,:,k).' * ref_data.reference_basis_derivatives(:,:,1));
    %Second coordinate 
    temp_2 = (Xij(:,:,k).' * ref_data.reference_basis_derivatives(:,:,2));
    
    map_derivatives(:,1:2,k) = temp_1.';
    map_derivatives(:,3:4,k) = temp_2.';
end 

  
%Generating the inverse map derivative:
%imap derivative row j contains [dphi_inv_1/dx_1,dphi_inv_2/dx_1,dphi_inv_1/dx_2,dphi_inv_2/dxi2]
%Generating determinant for inverse jacobian
determinant = zeros(size(ref_data.evaluation_points,1),1,fe_geometry.m);
imap_derivatives = zeros(size(ref_data.evaluation_points,1),4,fe_geometry.m);
Jac_2D = zeros(2,2);

for k = 1:fe_geometry.m
    %Check for first evaluation point for verification purposes
    %Jac_2D(1,1) = map_derivatives(1,1,k);
    %Jac_2D(2,1) = map_derivatives(1,2,k);
    %Jac_2D(1,2) = map_derivatives(1,3,k);
    %Jac_2D(2,2) = map_derivatives(1,4,k);
    %Jac_2D = inv(Jac_2D);
    %Jac_2D(:).'
    
    
    %Filling Jacobian determinant 
    determinant(:,1,k) = map_derivatives(:,1,k) .* map_derivatives(:,4,k) - map_derivatives(:,2,k) .* map_derivatives(:,3,k);

    
    %Generating derivatives (Inverse jacobian)
    dinv_11 = 1./determinant(:,1,k) .* map_derivatives(:,4,k);
    dinv_21 = -1./determinant(:,1,k) .* map_derivatives(:,2,k);
    dinv_12 = -1./determinant(:,1,k) .* map_derivatives(:,3,k);
    dinv_22 =  1./determinant(:,1,k) .* map_derivatives(:,1,k);
    
    %Filling derivatives
    imap_derivatives(:,1,k) = dinv_11;
    imap_derivatives(:,2,k) = dinv_21;
    imap_derivatives(:,3,k) = dinv_12;
    imap_derivatives(:,4,k) = dinv_22;
    imap_derivatives(1,:,k);

end 



geometric_map = struct('map',map,'map_derivatives',map_derivatives, ...
    'imap_derivatives',imap_derivatives);

    
end 