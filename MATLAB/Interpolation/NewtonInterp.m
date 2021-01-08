function P = NewtonInterp(x, y)
%
% function P = NewtonInterp(x, y)
%
% Newton Interpolation Polynomials
%   At each interpolation node, P(x_i) = y_i
%
% Input -
% x, y: Nx1 vectors, (x_i, y_i): interpolation nodes for any i = 1,2,...,N
%
% Output -
% P: struct,
%   P.x: first dimension of interpolation nodes
%   P.y: second dimension of interpolation nodes
%   P.d: divided differences of given data (x,y) = (x,f(x))
%   P.P: function handels, Newton interpolation polynomials
%   P.B: function handles, base funciton of Newton interpolation

%% Parameters check
if size(x,1) ~= size(y,1) || size(x,2) ~= 1 || size(y,2) ~= 1
    error('The shape of vector x and y need to be both "N x 1".')
end

if numel(x) ~= numel(unique(x))
    error('The elements of x must be unique.')
end

%% Initialization
P.x = x;
P.y = y;
N = numel(x);
P.B = cell(N,1);
P.P = sym('0');

%% Compute divided differences
P.d = divdif(x,y);

%% Generate Newton Interpolation Polynomials
for k = 1 : N
    P.B{k} = GNIPBase(k, x);
    P.P = P.P + P.d(k) * P.B{k};
    P.B{k} = matlabFunction(P.B{k});
end

%% Return
P.P = matlabFunction(P.P);

end

%% Generate kth Newton interpolation basis function
function pB = GNIPBase(k,X)
%
% function pB = GNIPBase(k,X)
%
%   Generate kth Newton interpolation basis function

%%
x = sym('x');
pB = sym('1');

if k == 1
    return;
end

%% Compute the LIP
for m = 1 : k-1
    pB = pB * (x - X(m));
end

end