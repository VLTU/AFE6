function ref_data = create_ref_data(nq1, p, kind)
    % MAKE SURE THAT YOUR MATLAB PATH CONTAINS THE FOLDERS:
    % ./quadrature
    % ./b_splines

    % INPUT
    % nq1: number of quadrature points in each direction
    % p: polynomial degree (p1,p2)
    % kind: 'plot' or 'integrate'

    % OUTPUT
    % ref_data:
    % - p
    % - evaluation_points
    % - quadrature_weights (empty if kind equals 'plot')
    % - reference_basis
    % - reference_basis_derivatives

    %% build all fields
    if ~exist('kind','var')
        kind = 'integrate';
    end

    if strcmp(kind,'plot')
        x = linspace(0, 1, nq1);
        x = x(:);
        quadrature_weights = [];

    elseif strcmp(kind,'integrate')
        [x,w] = lgwt(nq1, 0, 1);
        x = flipud(x); w = flipud(w);
        quadrature_weights = kron(w,w);

    end
    y = x;
    tmpx = B_ders_basis_funs_global(x, ...
                                   p(1), ...
                                   [zeros(1,p(1)+1) ones(1,p(1)+1)], ...
                                   1);
    tmpy = B_ders_basis_funs_global(y, ...
                                   p(2), ...
                                   [zeros(1,p(2)+1) ones(1,p(2)+1)], ...
                                   1);

    % evaluate tensor product basis and derivatives
    reference_basis = kron(tmpy(:,:,1),tmpx(:,:,1));
    reference_basis_derivatives = cat(3, ...
                                      kron(tmpy(:,:,1),tmpx(:,:,2)), ... % dx
                                      kron(tmpy(:,:,2),tmpx(:,:,1))); % dy

    % create evaluation point grid and save
    [y, x] = meshgrid(y, x);
    evaluation_points = [x(:), y(:)];

    %% save and return
    ref_data = struct('p', p, ...
                      'evaluation_points', evaluation_points, ...
                      'quadrature_weights', quadrature_weights, ...
                      'reference_basis', reference_basis, ...
                      'reference_basis_derivatives', reference_basis_derivatives);

end