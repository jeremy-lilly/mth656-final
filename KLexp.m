function Bx = KLexp(xi, x)
% DESCRIPTION: Computes the KL expansion of a Brownian Bridge Bx and
%   returns abs(Bx) + 1
% INPUTS:
%   xi: A 1D of vector of realizations of the N(0, 1) random variables
%       used in the KL expansion.
%   x: Values to evaluate the KL expansion at.
% RETURNS:
%   next_qdf: The degrees of freedom at the next time step.


    Bx =  0;
    for n = 1:size(xi)  % compute truncated expansion
        Bx = Bx + xi(n)*(1/(pi*n)).*sin((pi*n).*x);  % brownian bridge
    end

    Bx = abs(Bx) + 1;


end
