% Interpolation of Runge function
%   Runge function: f(x) = 1 / (1 + x^2)
%   Interpolation method:
%       Polynomial Interpolation,
%       Cubic Spline Interpolation,
%       Polynomial Interpolation using Chebyshev zeros

clc, clear, close all
%% Initialization
N = 21;
bd = [-1, 1];
runge = @(x) 1 ./ (1 + 25 * x.^2);
% first derivative of Runge function
rungeD1 = @(x) -50 * x ./ ((1 + 25 * x.^2)).^2;

% equal-space nodes in range [bd(1), bd(2)]
x1 = linspace(bd(1), bd(2), N)';
y1 = runge(x1);

% Chebyshev zeros in range [bd(1), bd(2)]
% in range [0, 1]
x2 = cos((2*(1:N) - ones(1,N)) / (2*N) * pi);
x2 = x2';
% Transform to range [bd(1), bd(2)]
x2 = ((bd(2)-bd(1)) * x2 + bd(1) + bd(2)) / 2;
y2 = runge(x2);

%% Interpolation

% ---- Polynomial (Newton) ------------------------------------------------
% equal-space nodes
P1 = NewtonInterp(x1, y1);

% Chebyshev zeros
P2 = NewtonInterp(x2, y2);

% ---- Cubic Spline -------------------------------------------------------
S = CubicSpline(x1, y1, rungeD1(x1(1)), rungeD1(x1(end)));

%% Display
x3 = (bd(1):0.005:bd(2))';

figure, set(gcf, 'outerposition', get(0,'screensize'));
plot(x3,P1.P(x3), 'LineWidth', 2), hold on
plot(x3,P2.P(x3), 'LineWidth', 2), hold on
plot(x3,ComputeSpline(S,x3), 'LineWidth', 2), hold on
plot(x3, runge(x3), 'LineWidth', 2), hold on
scatter(x1,y1, 'Filled'), hold on
scatter(x2,y2,'Filled')
legend('Newton, equal-space nodes',...
    'Newton, Chebyshev zeros',...
    'Cubic Spline, equal-space nodes',...
    'Runge function',...
    'equal-space nodes',...
    'Chebyshev zeros')
set(gca,'FontSize',20);
title(['Interpolation of Runge Function, N = ',num2str(N-1)], 'Fontsize', 20);
saveas(gcf,['results\N',num2str(N),'-Interpolation of Runge Function.png'])

if N > 12
    figure, set(gcf, 'outerposition', get(0,'screensize'));
    plot(x3,P1.P(x3), 'LineWidth', 2), hold on
    plot(x3, runge(x3), 'LineWidth', 2), hold on
    scatter(x1,y1,'Filled')
    legend('Newton, equal-space nodes',...
        'Runge function',...
        'equal-space nodes')
    set(gca,'FontSize',16);
    title(['Newton (equal-space nodes), N = ',num2str(N-1)], 'Fontsize', 20);
    saveas(gcf,['results\N',num2str(N),'-Newton (equal-space nodes).png'])

    figure, set(gcf, 'outerposition', get(0,'screensize'));
    plot(x3,P2.P(x3), 'LineWidth', 2), hold on
    plot(x3,ComputeSpline(S,x3), 'LineWidth', 2), hold on
    plot(x3, runge(x3), 'LineWidth', 2), hold on
    scatter(x1,y1,'Filled'), hold on
    scatter(x2,y2,'Filled')
    legend('Newton, Chebyshev zeros',...
        'Cubic Spline, equal-space nodes',...
        'Runge function',...
        'equal-space nodes',...
        'Chebyshev zeros')
    set(gca,'FontSize',16);
    title(['Newton (Chebyshev zeros) and Spline (equal-space nodes) N = ',num2str(N-1)], 'Fontsize', 20);
    saveas(gcf,['results\N',num2str(N),'-Newton (Chebyshev zeros) and Spline (equal-space nodes).png'])
end

%% error
t1 = (x1(N) + x1(N-1)) / 2;
t2 = (x2(N) + x2(N-1)) / 2;

f1 = runge(t1);
pe1 = P1.P(t1);
pc1 = P2.P(t2);
se1 = ComputeSpline(S,t1);

er_pe1 = f1 - pe1;
er_pc1 = f1 - pc1;
er_se1 = f1 - se1;

disp([num2str(f1), ' | ', num2str(pe1), ' | ' , num2str(er_pe1), ' |'])
disp([num2str(f1), ' | ', num2str(pc1), ' | ' , num2str(er_pc1), ' |'])
disp([num2str(f1), ' | ', num2str(se1), ' | ' , num2str(er_se1), ' |'])


