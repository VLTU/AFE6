
%Setting parameters 
%quadrature points in single dimension Note: ref has nq = nq_1^2
nq_1 = 2;
%polynomial bi-degree 
p = [2,2];


geometry_file = load("star3.mat");
ref_data = create_ref_data(nq_1,p,'plot')
%geometry_file.fe_geometry.map_coefficients

A = create_geometric_map(geometry_file.fe_geometry,ref_data);