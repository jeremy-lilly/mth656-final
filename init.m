function [P_gq, u_e, u_gq, qdf0] = init(J, xjs, max_deg_qh, xi, w, L, u)
% DESCRIPTION: Compute the initial degrees of freedom for the DG method.
% INPUTS:
%   J: Number of grid cells.
%   xjs: Centers of grid cells.
%   max_deg_qh: Maximum degree of the piece-wise polynomial approximation
%       q_h to compute.
%   xi: Gaussian Quadrature points in [-1, 1].
%   w: Gaussian Quadrature wieghts corresponding to xi.
%   u: Function defining fluid velocity on [0,1].
% RETURNS:
%   u_e: Fluid velocity at cell edges.
%   u_gq: Fluid velocity at GQ points.
%   P_gq: 2D matrix of Legendre polynomials, rows correspond to degrees of
%       the polynomials.
%   qdf0: 2D matrix containing the initial degrees of freedom; 
%       rows correspond to modes, columns to cells.


    % create matrix of x values
    x_gq = xMatrixOfGQPoints(xjs, xi);
    
    % get fluid velocity at GQ points
    u_gq = u(x_gq);
    
    % get fluid velocity at cell edges
    dx = xjs(2) - xjs(1);
    edges = xjs - dx/2;
    edges = [edges (edges(end) + dx)];
    u_e = u(edges);
    

    % build initial function at quadrature points from each grid cell
    Q_gq = cos((pi/2)*((x_gq - 0.5)/L)).^2 .* (abs(x_gq - 0.5) < L);

    % get legendre polynomials at quadrature points
    P_gq = legendrePolys(max_deg_qh, xi);

    % calculate inital values for degrees of freedom using gaussian 
    % quadrature
    qdf0 = zeros(max_deg_qh+1, J);  % (modes, cells)
    for m = 1:max_deg_qh+1
        for j = 1:J
            qdf0(m, j) = ((2*m - 1)/2) * sum(w.*Q_gq(j, :).*P_gq(m, :));
        end
    end


end
