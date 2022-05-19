function N = B_ders_basis_funs_global(u, p, U, varargin)

    % Compute the non-vanishing basisfunction and it's derivative
    % Input: u, p, U
    % Output: ders

    % number of basis functions - 1
    n = length(U) - p - 2;

    % number of evaluation points
    np = length(u);

    % allocate memory for 1 derivative
    if isempty(varargin)
        nDer = 1;
    else
        nDer = varargin{1};
    end
    N = zeros(n+1,np,nDer+1);

    % loop over evaluation points
    for i = 1:length(u)
        % point
        u_i = u(i);
        % evalaute derivatives
        [span,B] = B_ders_basis_funs(u_i,p,U,nDer);
        % which basis functions do these correspond to?
        ibasis = ((span+1)-p):(span+1);
        N(ibasis,i,:) = B;
    end
end