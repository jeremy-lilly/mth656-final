function [finex, fineQ, q_h] = dg(M, int, J, T, N, u)
% AUTHOR: Jeremy Lilly
% DESCRIPTION: Computes approximations of the solution q(x, t) to
%   q_t + (uq)_x = 0 on [a,b] where
%
%       q(x, 0) = cos( pi/2 * (x - 0.5)/L )     when |x - 0.5| < L
%               = 0                             otherwise
%
% INPUTS:
%   M: Degree of polynomial approximations.
%   int: Spatial interval of the form [a, b].
%   J: Number of grid cells to divide [0, 1] into.
%   T: Max time to step out to.
%   N: dt = 1/N.
%   u: Function defining fluid velocity on [0,1].
% RETURNS:
%   finex: 2D matrix of points that q_h is calculated at, each row
%       corresponds to a grid cell.
%   fineQ: 2D matrix where each row is the inital condition on the i-th
%       grid cell.
%   q_h: 2D matrix where each row is the solution on the i-th grid cell.


% define constants
dx = (int(2) - int(1))/J;  % length of grid cells
L = 0.2;  % constant used in definition of initial function
fineN = 100;  % number of points in each grid cell to plot at

% define quadrature points and weights
xi = [-0.906179845938664, ...
      -0.538469310105683, ...
       0, ...
       0.538469310105683, ...
       0.906179845938664];
w = [0.236926885056189, ...
     0.478628670499366, ...
     0.568888888888889, ...
     0.478628670499366, ...
     0.236926885056189];


% centers of grid cells
xjs = int(1) + dx/2 + dx*((1:J) - 1);

% initialize the degrees of freedom
[P_gq, u_e, u_gq, qdf] = init(J, xjs, M, xi, w, L, u);

% calculate values of derivatives of Legendre polys at GQ points
Pprime_gq = legendrePolysPrime(P_gq, xi);

% values of x to plot at
finex = fineXMatrix(xjs, fineN);

% Legendre polys evaluated at points in [-1, 1] corresponding to those
% in finex
fineP = legendrePolys(M, -1:2/(fineN-1):1);

% initial function plotted at points from finex
fineQ = zeros(J, fineN);
for j = 1:J
    fineQ(j, :) = cos((pi/2)*((finex(j, :) - 0.5)/L)).^2 ...
        .* (abs(finex(j, :) - 0.5) < L);
end


dt = 1/N;
% time-stepping
for i = 1:T*N  
    qdf = RKstep(qdf, P_gq, Pprime_gq, w, dx, u_e, u_gq, dt);    
end

% construct solution from degrees of freedom
q_h = qh(qdf, fineP);
%q_h = @(P) qh(qdf, P);
