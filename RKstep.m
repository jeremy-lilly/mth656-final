function next_qdf = RKstep(qdf, P, Pprime, w, dx, u_e, u_gq, dt)
% DESCRIPTION: Advances one step in time using the classical fourth-order
%   four-stage Runge-Kutta method.
% INPUTS:
%   qdf: 2D matrix of the degrees of freedom.
%   P_gq: 2D matrix of Legendre Polynomials eveluated at Gaussian 
%       Quadrature points.
%   Pprime_gq: 2D matrix of derivatives of Legendre Polynomials eveluated 
%       at Gaussian Quadrature points.
%   w: Gaussian Quadrature weights.
%   dx: Cell length.
%   u: Constant fluid velocity.
%   dt: size of time-step.
% RETURNS:
%   next_qdf: The degrees of freedom at the next time step.

    
    Q0 = qdf;
    
    F0 = F(Q0, P, Pprime, w, dx, u_e, u_gq);
    Q1 = qdf + (dt/2)*F0;
    
    F1 = F(Q1, P, Pprime, w, dx, u_e, u_gq);
    Q2 = qdf + (dt/2)*F1;
    
    F2 = F(Q2, P, Pprime, w, dx, u_e, u_gq);
    Q3 = qdf + dt*F2;
    
    F3 = F(Q3, P, Pprime, w, dx, u_e, u_gq);
    next_qdf = qdf + (dt/6)*(F0 + 2*F1 + 2*F2 + F3);


end
