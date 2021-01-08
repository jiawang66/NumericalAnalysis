function S = CubicSpline(X, Y, d1, d2)
%
% function S = CubicSpline(X, Y, d1, d2)
%
%   Cubic Spline Interpolation with type 1 boundary conditions (known
%   1st derivative at endpoints)
%
% Input -
% X, Y: Nx1 vectors, (x_i, y_i): interpolation nodes for any i = 1,2,...,N
% d1:= S'(x_1)
% d2:= S'(x_N)
%
% Output -
% S: struct,
%   S.x: first dimension of interpolation nodes
%   S.y: second dimension of interpolation nodes
%   S.S: function handle, cubic Spline Interpolation polynimals

%% Parameters check
if nargin ~= 4
    error('Error! 4 parameters needed.')
end

if size(X,1) ~= size(Y,1) || size(X,2) ~= 1 || size(Y,2) ~= 1
    error('The shape of vector x and y need to be both "N x 1".')
end

if numel(X) ~= numel(unique(X))
    error('The elements of x must be unique.')
end

%% Initialization
N = numel(X);
h = zeros(N-1,1);
mu = zeros(N-1,1);
lambda = zeros(N-1,1);
d = zeros(N,1);

% Sort ascending according to X
[X, indx] = sort(X);
Y = Y(indx);

S.x = X;
S.y = Y;
S.S = cell(N,1);

%% Compute mu, lambda, and d
lambda(1) = 1;
mu(N-1) = 1;

% ---- compute h ----------------------------------------------------------
for k = 1 : N-1
    h(k) = X(k+1) - X(k);
end

% ---- compute mu and lambda ----------------------------------------------
for k = 2 : N-1
    lambda(k) = h(k) / (h(k) + h(k-1));
    mu(k-1) = h(k-1) / (h(k) + h(k-1));
end

% ---- compute d ----------------------------------------------------------
% 1st divided differences
div = Y;

for k = N: -1: 2
    div(k) = (div(k) - div(k-1)) / (X(k) - X(k-1));
end
d(1) = 6 / h(1) * (div(2) - d1);
d(N) = 6 / h(N-1) * (d2 - div(N));

% 2nd divided differences
for k = N: -1: 3
    div(k) = (div(k) - div(k-1)) / (X(k) - X(k-2));
end
d(2:N-1) = 6 * div(3:N);

%% Compute M
% ---- Transform mu and lambda into a tridiagonal matrix ------------------
LE_path = '..\LinearEquation\';
addpath(LE_path);
A = diag(2*ones(N,1)) + diag(lambda,1) + diag(mu,-1);

% ---- Solve A * M = d by solving a tridiagnoal linear system -------------
M = LE_tridiag(A,d);

%% Generate cubic spline interpolation polynomial
x = sym('x');

for k = 1 : N-1
    tmp = M(k) * (X(k+1) - x)^3 / 6 / h(k)...
        + M(k+1) * (x - X(k))^3 / 6 / h(k)...
        + (Y(k) - M(k) * h(k)^2 / 6) * (X(k+1) - x) / h(k)...
        + (Y(k+1) - M(k+1) * h(k)^2 / 6) * (x - X(k)) / h(k);
    S.S{k} = matlabFunction(tmp);
end

end