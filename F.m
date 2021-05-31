function forcing = F(qdf, P_gq, Pprime_gq, w, dx, u_e, u_gq)
% DESCRIPTION: Calculates the forcing term for our system of ODEs.
% INPUTS:
%   qdf: 2D matrix of the degrees of freedom.
%   P_gq: 2D matrix of Legendre Polynomials eveluated at Gaussian 
%       Quadrature points.
%   Pprime_gq: 2D matrix of derivatives of Legendre Polynomials eveluated 
%       at Gaussian Quadrature points.
%   w: Gaussian Quadrature weights.
%   dx: Cell length.
%   u_e: Fluid velocity at cell edges, should have length J+1.
%   u_gq: 2D matrix of fluid velocity evaluated at Gaussian Quadrature
%       points.
% RETURNS:
%   forcing: 2D matrix of the forcing terms where the rows correspond
%       to modes and columns to cells.


    % get meta constants
    [Mtilde, J] = size(qdf);
    M = Mtilde-1;  % max mode number (number of modes - 1)
    [~, N] = size(P_gq);  % number of GQ points
    

    % calculate fluxes at interior cell edges
    fhat = zeros(1, J+1);
    for j = 2:J
        if u_e(j) > 0
            fhat(j) = u_e(j) * sum(qdf(:, j-1));
        else
            fhat(j) = u_e(j) * sum((-1).^(0:M) * qdf(:, j));
        end
    end
    
    % no boundary restrictions
%     for j = [1 J+1]
%         if u_e(j) > 0
%             fhat(j) = u_e(j) * sum(qdf(:, j-1));
%         else
%             fhat(j) = u_e(j) * sum((-1).^(0:M) * qdf(:, j));
%         end
%     end
    
    % periodic boundary conditions (boundary cell edges)
    if u_e(end) > 0
        fhat(end) = u_e(end) * sum(qdf(:, end));
        fhat(1) = fhat(end);
    else
        fhat(1) = u_e(1) * sum((-1).^(0:M) * qdf(:, 1));
        fhat(end) = fhat(1);
    end
    
    
    % calculate forcing terms
    forcing = zeros(M+1, J);
    for m = 1:M+1
        for j = 1:J
            flux_term = ((-1)^(m-1))*fhat(j) - fhat(j+1);
            
            
            integral_term = 0;
            for i = 1:N
               integral_term = integral_term + (w(i) * u_gq(j, i) * ...
                   sum(qdf(:, j).*P_gq(:, i)) * Pprime_gq(m, i)); 
            end

            
            forcing(m, j) = ((2*m - 1)/dx) * (flux_term + integral_term);
        end
    end


end
