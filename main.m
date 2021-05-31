% AUTHOR: Jeremy Lilly
% DESCRIPTION: Use discontinous Galerkin, KL expansion, and Clenshaw-Curtis 
% stochatic collocation to compute the expected value of the solution to
% the advection equation
%   q_t + (U_x q)_x = 0
% for a random process U_x.

close all; clear;

rng('shuffle')  % seed rng
make_plots = 1;
save_data = 0;


%% Define parameters

% set number KL terms to keep
NKL = 5;
percent_var = sum( 1/pi^2 * 1./(1:NKL).^2 ) / (1/6);

% set 1D level for Clenshaw-Curtis
cc_level = 2;

% define parameters for DG
M = 5;  % degree of polynomial approximations
I = [0 1];  % spatial interval to compute soln on
J = 50;  % number of grid cells
T = 1;  % time T to compute soln to 
N = 900;  % dt = 1/N


%% Perform quadrature

% get quadrature nodes
[quad_grid, quad_weights] = cc_grid_dataset(NKL, 2^(cc_level-1)+1);
%[quad_grid, quad_weights] = sparse_grid_cc_dataset(NKL, cc_level-1);

q_exp = 0;

for i = 1:size(quad_grid, 2)
    
    % set advection velocity to KLexp at current quadrature node
    Ux = @(x) KLexp(quad_grid(:, i), x);
    
    % use DG to get soln at current quadrature node
    [x, q_0, q_h] = dg(M, I, J, T, N, Ux);
    
    % PLOTS for each quadrature node, plot the realization of the soln
    if make_plots == 1
        figure(1); hold on;
        if i == 1  % plot inital condidtion
            for j = 1:J
                init_plot = plot(x(j, :), q_0(j, :), ...
                                 'k--', 'LineWidth', 1.5);
            end
        end
        
        c1 = rand(1); c2 = rand(1); c3 = rand(1);
        for j = 1:J
            plot(x(j, :), q_h(j, :), 'Color', [c1 c2 c3])
        end
        hold off;
    end
    % PLOTS end
    
    % calculate expected value of q_h using quadrature
    q_exp = q_exp + quad_weights(i)*q_h;
end

if make_plots == 1
    figure(1); hold on;
    for j = 1:J
        exp_plot = plot(x(j, :), q_exp(j, :),'r-.', 'LineWidth', 2.0);
    end
    title('Realizations at quadrature nodes', ...
          ['J = ', num2str(J), ', N = ', num2str(NKL)]);
    legend([init_plot exp_plot], {'Inital Condition', 'Expected Value'});
    hold off;
end

if save_data == 1
    save(['m', num2str(M), 'j', num2str(J), 'nkl', ...
          num2str(NKL), 'l', num2str(cc_level) '.mat'], ...
         'q_exp', 'M', 'J', 'NKL', 'x');
end
 
