function [x,count] = my_Jacobi(A, b, r, n)
%
% [x,count] = my_Jacobi(A, b, r, n)
%
% Solve linear equation using iterative method Jacobi
%   A = D - L - U
%   D * x_{k+1} = (L+U) * x_{k} + b
%
% Input -
% A: NxN matrix
% b: NxS matrix
% r: maximum residule error, estimated with infinite norm
% n: natural number, the maximum number of iteration
%
% Output -
% x: approximate solution of Ax=b
% count: iteration steps
%
% Example -
% >> x = my_Jacobi(A, b, 1e-12, 1000)

%% Parameter check
if nargin < 2
    error("Error! At least two parameters needed.")
elseif nargin == 2
    r = 1e-6;
    n = 1000;
elseif nargin == 3
    n = 1000;
end

if ~ismatrix(A) || ~ismatrix(b)
    error("Error! First two parameters shoube be matrixs.")
end

if size(A,1) ~= size(b,1)
    error("Error! The first dimension of A and b should be equal.")
end

%% Initialization
[rows,~] = size(A);
[~,s] = size(b);
x0 = zeros(rows,s);
x1 = x0;
count = 0;
eps = 1e-32;

% first iteration
count = count + 1;
for i = 1:rows
    if A(i,i) ~= 0
        x1(i,:) = x0(i,:) + (-1 * A(i,:) * x0 + b(i,:)) / A(i,i);
    else
        error("Error! The pivot can not be zero.")
    end
end

%% Iteration of Jacobi
while max(abs(x1-x0))/(max(abs(x1))+eps) >= r && count < n
    x0 = x1;
    count = count + 1;
    
    for i = 1:rows
        x1(i,:) = x0(i,:) + (-1 * A(i,:) * x0 + b(i,:)) / A(i,i);
    end
end

%% Return
if count >= n
    disp("Fail to obtain a solution satisfying the given conditions.")
end
x = x1;

end