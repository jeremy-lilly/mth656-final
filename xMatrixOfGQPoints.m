function x = xMatrixOfGQPoints(xjs, xi)
% DESCRIPTION:
%   Generates a matrix of points in grid cells corresponding to 
%   quadrature points from [-1, 1]. Assumes uniform cell length.
% INPUTS:
%   xjs: 1D vector containing centers of grid cells.
%   xi: 1D vector containing GQ points in [-1, 1].
% RETURNS:
%   x: 2D matrix whose rows correspond to different grid cells.


    % get meta data
    J = length(xjs);
    N = length(xi);
    dx = xjs(2) - xjs(1);
    
    % build x
    x = zeros(J, N);
    for j = (1:J)
        x(j, :) = xjs(j) + (dx/2)*xi;
    end


end
