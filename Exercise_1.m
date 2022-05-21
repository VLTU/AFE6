
%Setting parameters 
%quadrature points in single dimension Note: ref has nq = nq_1^2
nq_1 = 30;
n_q_1 = nq_1;
%polynomial bi-degree 
p = [2,2];


geometry_file = load("star5.mat");
geometry_file.p
ref_data = create_ref_data(nq_1,p,'plot');''
%geometry_file.fe_geometry.map_coefficients

A = create_geometric_map(geometry_file.fe_geometry,ref_data);
%Plotting basis functions 
basis_index = 3;
coordinates = A.map;
ref_data.reference_basis;

geometry_file.fe_space.support_and_extraction.supported_bases
S = full(geometry_file.fe_space.support_and_extraction(1).extraction_coefficients)



for k=1:geometry_file.fe_geometry.m
 relevant_space = geometry_file.fe_space.support_and_extraction(k);
 meta_index = find(relevant_space.supported_bases == basis_index);
 full(relevant_space.extraction_coefficients(meta_index,:))
 ref_data.reference_basis
 data = relevant_space.extraction_coefficients(meta_index,:) * ref_data.reference_basis
 x = A.map(:,1,k).'
 y = A.map(:,2,k).'
 size1d = sqrt(size(x,2))

 x = reshape(x,size1d,size1d)
 y = reshape(y,size1d,size1d)
 data  = reshape(data,size1d,size1d)
 figure(k)
 pcolor(x,y,data);

 colorbar;
 shading interp
  
 %ref_data.evaluation_points
 %A.map(:,:,k)
end


%A.map
%A.map_derivatives
