% Solve Poisson Equation
% -(\frac{\partial^2{u}}{\partial{x}^2} + \frac{\partial^2{u}}{\partial{y}^2}) = f(x,y)
% u(0,y) = 1, u(1,y) = e^y
% u(x,0) = 1, u(x,1) = e^x
%
% A * u + d = b
%   A: five-point difference scheme of Poisson equation
%   u: expected solution
%   d: boundary value, u(0,y), u(1,y), u(x,0), u(x,1)
%   b: f(x,y)

clc, clear, close all

%% Generate coefficient matrix A and boundary d
N = 15;
bd = [0,1];
h = (bd(2)-bd(1))/(N+1);
u_0y = @(x,y)1;
u_1y = @(x,y)exp(y);
u_x0 = @(x,y)1;
u_x1 = @(x,y)exp(x);
[A,d] = my_poisson_matrix(h,bd,u_0y,u_1y,u_x0,u_x1);

figure, set(gcf, 'outerposition', [11,11,1034,1034]);
imagesc(A), colorbar;
title(['Five-point difference scheme of Poisson equation (N = ', num2str(N),')'])
saveas(gcf, ['results\N', num2str(N),'-1-PoissonMatrix.png']);

%% Generate b
f = @(x,y) (x.^2 + y.^2) .* exp(x.*y);

ind = 1:N;
indx = repmat(ind',1,N)*h + ones(N)*bd(1);
indy = repmat(ind,N,1)*h + ones(N)*bd(1);

b = f(indx, indy);
b = reshape(b,N^2,1) * h^2;

%% Real solution
% u(x,y) = e^{x*y};
u_xy = @(x,y)exp(x.*y);
u0 = u_xy(indx,indy);

%% Solve A*x = b-d with LU Doolittle decomposition method
[x,y] = LE_doolittle(A,b-d);

absE = abs(A*x-(b-d));
max(absE(:))

u1 = reshape(x,N,N); max(max(abs(u0-u1)./u0))
figure, set(gcf, 'outerposition', [11,11,1180,1034]);
subplot(221), imagesc(u0), colorbar, title('real solution')
subplot(222), imagesc(u1), colorbar, title('solution by LU')
subplot(223), imagesc(abs(u0-u1)./u0), colorbar,
title('relative error between real solution and solution by LU')
subplot(224), imagesc(reshape(absE,N,N)./reshape(b-d,N,N)), colorbar, 
title('relative error of LU')
saveas(gcf, ['results\N', num2str(N),'-1-LU.png']);

%% Solve A*x = b-d with Jacobi iteration method
[x,count] = my_Jacobi(A,b-d,1e-15,1e6);
disp(['Jacobi: ', num2str(count)])

absE = abs(A*x-(b-d));
max(absE(:))

u2 = reshape(x,N,N); max(max(abs(u0-u2)./u0))
figure, set(gcf, 'outerposition', [11,11,1180,1034]);
subplot(221), imagesc(u0), colorbar, title('real solution')
subplot(222), imagesc(u2), colorbar, title('solution by Jacobi')
subplot(223), imagesc(abs(u0-u2)./u0), colorbar,
title('relative error between real solution and solution by Jacobi')
subplot(224), imagesc(reshape(absE,N,N)./reshape(b-d,N,N)), colorbar,
title('relative error of Jacobi')
saveas(gcf, ['results\N', num2str(N),'-1-Jacobi.png']);

%% Solve A*x = b-d with Gauss-Seidel iteration method
[x,count] = my_SOR(A,b-d,1,1e-15,1e6);
disp(['GS: ', num2str(count)])

absE = abs(A*x-(b-d));
max(absE(:))

u3 = reshape(x,N,N); max(max(abs(u0-u3)./u0))
figure, set(gcf, 'outerposition', [11,11,1180,1034]);
subplot(221), imagesc(u0), colorbar, title('real solution')
subplot(222), imagesc(u3), colorbar, title('solution by GS')
subplot(223), imagesc(abs(u0-u3)./u0), colorbar,
title('relative error between real solution and solution by GS')
subplot(224), imagesc(reshape(absE,N,N)./reshape(b-d,N,N)), colorbar,
title('relative error of GS')
saveas(gcf, ['results\N', num2str(N),'-1-GS.png']);

%% SOR, with different values of omega
omega = 1.7;
[x,count] = my_SOR(A,b-d,omega,1e-15,1e6);
disp(['SOR: ', num2str(count)])

absE = abs(A*x-(b-d));
max(absE(:))

u4 = reshape(x,N,N); max(max(abs(u0-u4)./u0))
figure, set(gcf, 'outerposition', [11,11,1180,1034]);
subplot(221), imagesc(u0), colorbar, title('real solution')
subplot(222), imagesc(u4), colorbar, title(['solution by SOR (w = ',num2str(omega),')'])
subplot(223), imagesc(abs(u0-u4)./u0), colorbar,
title('relative error between real solution and solution by SOR')
subplot(224), imagesc(reshape(absE,N,N)./reshape(b-d,N,N)), colorbar,
title('relative error of SOR')
saveas(gcf, ['results\N', num2str(N),'-1-SOR-w',num2str(omega),'.png']);

%% Transform mat to cell
AC = mat2cell(A, N*ones(N,1), N*ones(1,N));
bdc = mat2cell(b-d, N*ones(N,1));

%% Block Jacobi
[x,count] = my_BJacobi(AC, bdc, 1e-15, 1e6);
disp(['BJacobi: ', num2str(count)])

absE = abs(A*cell2mat(x) - (b-d));
max(absE(:))

u5 = reshape(cell2mat(x),N,N); max(max(abs(u0-u5)./u0))
figure, set(gcf, 'outerposition', [11,11,1180,1034]);
subplot(221), imagesc(u0), colorbar, title('real solution')
subplot(222), imagesc(u5), colorbar, title('solution by BJacobi')
subplot(223), imagesc(abs(u0-u5)./u0), colorbar,
title('relative error between real solution and solution by BJacobi')
subplot(224), imagesc(reshape(absE,N,N)./reshape(b-d,N,N)), colorbar,
title('relative error of BJacobi')
saveas(gcf, ['results\N', num2str(N),'-1-BJacobi.png']);

%% Block Gauss-Seidel
omega = 1;
[x,count] = my_BSOR(AC, bdc, omega, 1e-15, 1e6);
disp(['BGS: ', num2str(count)])

absE = abs(A*cell2mat(x) - (b-d));
max(absE(:))

u6 = reshape(cell2mat(x),N,N); max(max(abs(u0-u6)./u0))
figure, set(gcf, 'outerposition', [11,11,1180,1034]);
subplot(221), imagesc(u0), colorbar, title('real solution')
subplot(222), imagesc(u6), colorbar, title('solution by BGS')
subplot(223), imagesc(abs(u0-u6)./u0), colorbar,
title('relative error between real solution and solution by BGS')
subplot(224), imagesc(reshape(absE,N,N)./reshape(b-d,N,N)), colorbar,
title('relative error of BGS')
saveas(gcf, ['results\N', num2str(N),'-1-BGS.png']);

%% BSOR, with different values of omega
omega = 0.6;
[x,count] = my_BSOR(AC, bdc, omega, 1e-15, 1e6);
disp(['BSOR: ', num2str(count)])

absE = abs(A*cell2mat(x) - (b-d));
max(absE(:))

u7 = reshape(cell2mat(x),N,N); max(max(abs(u0-u7)./u0))
figure, set(gcf, 'outerposition', [11,11,1180,1034]);
subplot(221), imagesc(u0), colorbar, title('real solution')
subplot(222), imagesc(u7), colorbar, title(['solution by BSOR (w = ',num2str(omega),')'])
subplot(223), imagesc(abs(u0-u7)./u0), colorbar,
title('relative error between real solution and solution by BSOR')
subplot(224), imagesc(reshape(absE,N,N)./reshape(b-d,N,N)), colorbar,
title('relative error of BSOR')
saveas(gcf, ['results\N', num2str(N),'-1-BSOR-w',num2str(omega),'.png']);

%% Conjugate gradient (GC)
x = pcg(A, b-d, 1e-15, 1e6); % MATLAB inner function

absE = abs(A*x - (b-d));
max(absE(:))

u8 = reshape(x,N,N); max(max(abs(u0-u8)./u0))
figure, set(gcf, 'outerposition', [11,11,1180,1034]);
subplot(221), imagesc(u0), colorbar, title('real solution')
subplot(222), imagesc(u8), colorbar, title('solution by GC')
subplot(223), imagesc(abs(u0-u8)./u0), colorbar, title('relative error')
subplot(224), imagesc(reshape(absE,N,N)./reshape(b-d,N,N)), colorbar,
title('relative error of GC')
saveas(gcf, ['results\N', num2str(N),'-1-GC.png']);
