function [q_h] = qh_func(qdf, x)
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
    
    x = reshape(x, M, []);
    
    [~, len] = size(x);
    
    Px = -1:2/(len-1):1;
    P = zeros(M+1, length(Px));
    P(1, :) = 1; P(2, :) = Px;
    for m = 2:max_deg
        P(m+1, :) = ((2*(m-1) + 1)*Px.*P(m, :) - (m-1)*P(m-1, :)) / m;
    end
    
    q_h = zeros([J, len]);
    for j = 1:J
        for m = 1:M+1
            q_h(j, :) = q_h(j, :) + qdf(m, j)*P(m, :);
        end
    end


end

