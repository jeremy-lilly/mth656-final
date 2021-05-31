% AUTHOR: Jeremy Lilly
% DESCRIPTION: Create a gif of the solution to the advection equation as
% solved by a DG method

close all; clear;

% define parameters for DG
M = 5;  % degree of polynomial approximations
I = [0 1];  % spatial interval to compute soln on
J = 20;  % number of grid cells
N = 350;  % dt = 1/N

u = @(x) ones(size(x));

fig = figure;
axis tight manual
filename = 'advection.gif';


for T = 0:0.05:1
    [x, q_0, q_h] = dg(M, I, J, T, N, u);
    
    clf(fig);
    
    hold on;
    for j = 1:J
        init_plot = plot(x(j, :), q_0(j, :), ...
            'k--', 'LineWidth', 1.5);
    end
    
    
    for j = 1:J
        plot(x(j, :), q_h(j, :), 'Color', [0.3 0.65 1])
    end
    title(['DG solution to {q_t + q_x = 0}'], ['t = ', num2str(T)]);
    xlabel('x')
    ylabel('{q_h(x, t)}')
    axis([0 1 -0.2 1.2])
    hold off;
    
    frame = getframe(fig);
    img = frame2im(frame);
    [imgind, color_map] = rgb2ind(img, 256);
    
    if T == 0
        imwrite(imgind, color_map, filename, ...
                'gif', 'LoopCount', Inf, 'DelayTime', 0.25)
    else
        imwrite(imgind, color_map, filename, ... 
                'WriteMode', 'append', 'DelayTime', 0.25)
    end
    
end