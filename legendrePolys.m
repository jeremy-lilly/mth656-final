function P = legendrePolys(max_deg, x)
% DESCRIPTION: Generate Legendre polynomials eveluated at points from x.
% INPUT:
%   max_deg: Largest degree Legendre Polynomial to compute.
%   x: 1D vector of x values to eveluate polynomials at.
% RETURNS:
%   P: 2D matrix of Legendre polynomialsevaluated at x, rows 
%       correspond to degrees of the polynomials.
    
    
    P = zeros(max_deg+1, length(x));
    P(1, :) = 1; P(2, :) = x;
    for m = 2:max_deg
        P(m+1, :) = ((2*(m-1) + 1)*x.*P(m, :) - (m-1)*P(m-1, :)) / m;
    end
    
    
end
