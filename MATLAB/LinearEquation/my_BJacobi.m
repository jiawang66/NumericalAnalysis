function [x,count] = my_BJacobi(A, b, r, n)
%
% function [x,count] = my_BJacobi(A, b, r, n)
%
% Block Jacobi solution method for block matrix
%
% Input -
% A: NxN cell matrix
% b: NxS cell matrix
% r: maximum residule error, estimated with infinite norm
% n: natural number, the maximum number of iteration
%
% Output -
% x: approximate solution of Ax=b
% count: iteration steps
%
% Example -
% >> x = my_BJacobi(A,b)

%% Parameter check
if nargin < 2
    error("Error! At least two parameters needed.")
elseif nargin == 2
    r = 1e-6;
    n = 1000;
elseif nargin == 3
    n = 1000;
end

if size(A,1) ~= size(b,1)
    error("Error! The first dimension of A and b should be equal.")
end

%% Initialization
[rows,~] = size(A);
[~,s] = size(b);
[p,~] = size(A{1,1});

x0 = zeros(rows*p,s);
x0 = mat2cell(x0, rows*ones(p,1));
x1 = x0;
count = 0;
eps = 1e-32;

% first iteration
count = count + 1;
for i = 1:rows
    x1{i,:} = x0{i,:} + A{i,i} \ (-1 * cellMultiply(A(i,:), x0) + b{i,:});
end

%% Iteration of Jacobi
while cellMax(cellAbs(cellPlus(x1,x0,-1))) / (cellMax(cellAbs(x1))+eps) >= r...
        && count < n
    x0 = x1;
    count = count + 1;
    
    for i = 1:rows
        x1{i,:} = x0{i,:} + A{i,i} \ (-1 * cellMultiply(A(i,:), x0) + b{i,:});
    end
end

%% Return
if count >= n
    disp("BJacobi: Fail to obtain a solution satisfying the given conditions.")
end
x = x1;

end