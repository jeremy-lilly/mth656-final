function q_h = qh(qdf, P)
% DESCRIPTION: Calculates the piecewise polynomial approximation q_h.
% INPUTS:
%   qdf: The degrees of freedom.
%   P: 2D matrix of Legendre polynomials evaluated at points corresponding
%       to those in x. Each row corresponds to a different degree.
% RETURNS:
%   qh: 2D matrix of qh eveluated at points from x.


    % get meta constants
    [Mtilde, J] = size(qdf);
    M = Mtilde - 1;
    
    [~, len] = size(P);
    
    q_h = zeros([J, len]);
    for j = 1:J
        for m = 1:M+1
            q_h(j, :) = q_h(j, :) + qdf(m, j)*P(m, :);
        end
    end
    
    
end
