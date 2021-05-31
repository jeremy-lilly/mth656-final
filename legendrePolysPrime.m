function P_prime = legendrePolysPrime(P, x)
% DESCRIPTION: Generate derivatives of Legendre polynomials eveluated 
%   at points from x.
% INPUT:
%   P: 2D martix of Legendre polynomials (where rows correspond to degree)
%       to calculate the derivatives of.
%   x: 1D vector of values that P corresponds to.
% RETURNS:
%   P_prime: 2D matrix of derivaties of Legendre polynomials, rows 
%       correspond to degrees of the polynomials.


    % get meta data
    [M, N] = size(P);
    
    P_prime = zeros(M, N);
    P_prime(1, :) = 0;
    for m = 2:M
        P_prime(m, :) = ((m-1)./(x.^2 - 1)) .* (x.*P(m, :) - P(m-1, :));
    end


end
