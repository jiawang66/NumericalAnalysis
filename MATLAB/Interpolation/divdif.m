function d = divdif(x, y)
%
% function d = divdif(x, y)
%
%   Calculate the divided differences of the given values in the vector x
%   and y. y are the values of some function f(x).
%
% Input -
% x: Nx1 vector, nodes
% y: Nx1 vector, y = f(x)
%
% Output -
% d: Nx1 vector, y[k] = f[x_1,..x_k], k = 1,2,...,N

%% Parameter check
if size(x,1) ~= size(y,1) || size(x,2) ~= 1 || size(y,2) ~= 1
    error('The shape of vector x and y need to be both "N x 1".')
end

if numel(x) ~= numel(unique(x))
    error('The elements of x must be unique.')
end

%% Compute
d = y;
N = numel(x);

for k = 2 : N
    for s = N: -1: k
        d(s) = (d(s) - d(s-1)) / (x(s) - x(s-k+1));
    end
end

end
