function x = fineXMatrix(xjs, N)
% DESCRIPTION: Creates a maxtix of x values, where each row in the
%   maxtirx corresponds to a grid cell (assumes grid cells have 
%   uniform length).
% INPUTS:
%   J: Number of grid cells.
%   N: Number of points per grid cell.
%   xjs: Centers of grid cells.
% RETURNS:
%   x: 2D matrix of x values, where each row corresponds to a different 
%   grid cell.


    % get number of grid cells and grid length
    J = length(xjs);
    dx = xjs(2) - xjs(1);
    
    x = zeros(J, N);
    for j = 1:J
        x(j, :) = (xjs(j) - dx/2):dx/(N - 1):(xjs(j) + dx/2);
    end


end
