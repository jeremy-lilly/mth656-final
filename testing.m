close all
clear all

M = 5;
N = 1000;
J = 20;
T = 1;
int = [0 1];

NKL = 100;

u = @(x) KLexp(100, x);

mat = [];

%for j = 1:1000
    [x, Q, q_h] = dg(M, int, J, T, N, u);
    mat = [mat ; reshape(q_h', 1, [])];
%end

%create plot
for j = 1:J
    plot(x(j, :), Q(j, :), '--', ...
         x(j, :), q_h(j, :))
    %axis([0 1 -0.2 1.2])
    title("M = " + num2str(M) + ...
        ", T = " + num2str(T) + ", \Delta t = 1/" + num2str(N))
    hold on;
end

var = std(mat);
avg = mean(mat);
x = reshape(x', 1, []);

%plot(x, avg, x, var)

%reshape(x', 1, [])
