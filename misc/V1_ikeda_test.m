clear all

% u = ikeda parameter
% option = what to plot
%  'trajectory' - plot trajectory of random starting points
%  'limit' - plot the last few iterations of random starting points
function ikeda(u, option)
points = 3000; % how many starting points
N = 10000; % how many iterations
Nlimit = 20; % plot these many last points for 'limit' option
 
x = rand(1, points);  % the random starting points
y = rand(1, points); 
 
%P_range=0:1:30;
T_range=0:1:5;

s=size(T_range,2);

for i =1:s
    T=T_range(i);
for j = 1:points

    X = compute_ikeda_trajectory(u, x(j), y(j), N,T);
     
    switch option
        case 'trajectory' % plot the trajectories of a bunch of points
            plot_ikeda_trajectory(X); hold on;

        case 'limit'
            plot_limit(X, Nlimit); hold on;

        otherwise
            disp('Not implemented');
    end
end
end
end

% Plot the last n points of the curve - to see end point or limit cycle
function plot_limit(X, n)
plot(X(end - n:end, 1), X(end - n:end, 2), 'ko');
end

% Plot the whole trajectory
function plot_ikeda_trajectory(X)
plot(X(:, 1), X(:, 2));
% hold on; plot(X(1,1), X(1,2), 'bo', 'markerfacecolor', 'g'); hold off
title('Map with Dynamic T [0,5]')
end

% u is the ikeda parameter
% x,y is the starting point
% N is the number of iterations

function [X] = compute_ikeda_trajectory(u, x, y, N,T)
n_kerr=0.55;
tao=18.5;
tao_theta=185;
sigma_FCD=7.2;
xi_T=0.074;
eta_lin=0.4;
eta_c=1;
alpha_TPA=0.11;
gamma_FCA=0.2;

X = zeros(N, 2);
X(1, :) = [x y];


P=0.25;
delta=-3;
%T=5;
n=3;
for j = 2:N
    
    %x1 = 1 + u.* (x.* cos(0.4 - 6. / (1 + x.^2 + y.^2)) - y.* sin(0.4 - 6. / (1 + x.^2 + y.^2)));
    %y1 = u.* (x.* sin(0.4 - 6. / (1 + x.^2 + y.^2)) + y.* cos(0.4 - 6. / (1 + x.^2 + y.^2)));
    x1=x + u.* (sqrt(P) - x.*(1+alpha_TPA.*(x^2+y^2) + gamma_FCA.*n) - y.*(delta-n_kerr.*(x^2+y^2) + (n+sigma_FCD.*n.^0.8) - T)  );
    y1=y + u.* (x.*(delta - n_kerr.*(x^2+y^2)+(n+sigma_FCD.*n.^0.8) - T) -y.*(1+alpha_TPA.*(x^2+y^2) + gamma_FCA.*n) );
    x = x1;
    y = y1;

    X(j, :) = [x y];
end
end

ikeda(0.9, 'trajectory')