function L = LagrangeInterp(x, y)
%
% function f = LagrangeInterp(x, y)
%
% Lagrange Interpolation Polynomials
%   At each interpolation node, L(x_i) = y_i
%
% Input -
% x, y: Nx1 vectors, (x_i, y_i): interpolation nodes for any i = 1,2,...,N
%
% Output -
% L: struct,
%   L.x: first dimension of interpolation nodes
%   L.y: second dimension of interpolation nodes
%   L.L: function handels, Lagrange interpolation polynomials
%   L.B: function handles, Lagrange interpolation basis function

%% Parameters check
if size(x,1) ~= size(y,1) || size(x,2) ~= 1 || size(y,2) ~= 1
    error('The shape of vector x and y need to be both "N x 1".')
end

if numel(x) ~= numel(unique(x))
    error('The elements of x must be unique.')
end

%% Initialization
L.x = x;
L.y = y;
N = numel(x);
L.B = cell(N,1);
L.L = sym('0');

%% Generate Lagrange interpolation polynomials
for k = 1 : N
    L.B{k} = GLIPBase(k, x);
    L.L = L.L + y(k) * L.B{k};
    L.B{k} = matlabFunction(L.B{k});
end

%% Return
L.L = matlabFunction(L.L);

end

%% Generate kth Lagrange interpolation basis function
function pB = GLIPBase(k,X)
%
% function pB = GLIPBase(k,X)
%
%   Generate kth base function of Lagrange interpolation polynomials

%%
x = sym('x');
pB = sym('1');
N = numel(X);
denominator = 1;

%% Compute the denominator
for m = 1 : N
    if m ~= k
        denominator = denominator * (X(k) - X(m));
    end
end

%% Compute the LIP
for m = 1 : N
    if m ~= k
        pB = pB * (x - X(m));
    end
end

%% Return
pB = pB / denominator;

end