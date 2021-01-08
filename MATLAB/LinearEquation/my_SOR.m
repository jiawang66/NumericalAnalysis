function [x,count] = my_SOR(A,b,omega, r, n, p, xstar)
%
% [x,count] = my_SOR(A,b,omega, r, n, p, xstar)
%
% Solve linear equation using iterative method SOR
%
% Input -
% A: NxN matrix
% b: NxS matrix
% omega: real number, a weighted factor, default omega = 1 (Gauss-Seidel)
% r: maximum residule error, estimated with infinite norm
% n: natural number, the maximum number of iteration
% p: Logical, if true, display detail
% xstar: precise solution
%
% Output -
% x: approximate solution of Ax=b
% count: iteration steps
%
% Example -
% >> x = my_SOR(A,b)
% >> x = my_SOR(A,b,omega,r,n)
% >> x = my_SOR(A,b,omega,r,n,true)

%% Parameter check
if nargin < 2
    error("Error! At least two parameters needed.")
elseif nargin == 2
    omega = 1;  % Gauss-Seidel iteration method
    r = 1e-5;
    n = 500;
elseif nargin == 3
    r = 1e-5;
    n = 500;
elseif nargin == 4
    n = 500;
end

if nargin < 6
    p = false;
end

if ~ismatrix(A) || ~ismatrix(b)
    error("Error! First two parameters shoube be matrixs.")
end

if size(A,1) ~= size(b,1)
    error("Error! The first dimension of A and b should be equal.")
end

[rows,~] = size(A);
[~,s] = size(b);

if nargin < 7
    xstar = zeros(rows,s);
end

flag = all(xstar==0);

%% initialization
x0 = zeros(rows,s);
x1 = x0;
count = 0;
eps = 1e-32;

% compute x1
count = count + 1;
for i = 1:rows
    if A(i,i) ~= 0
        x1(i,:) = x1(i,:) + omega * (-1 * A(i,:) * x1 + b(i,:)) / A(i,i);
    else
        error("Error! The pivot can not be zero.")
    end
end

display_detail(count,p,x1,x0,xstar)

%% iteration of SOR
if flag
    xwave = x0;
else
    xwave = xstar;
end

while max(abs(xwave - x1))/(max(abs(x1))+eps) >= r && count < n
    x0 = x1;
    count = count + 1;
    if flag
        xwave = x0;
    end
    
    for i = 1:rows
        x1(i,:) = x1(i,:) + omega * (-1 * A(i,:) * x1 + b(i,:)) / A(i,i);
    end
    
    display_detail(count,p,x1,x0,xstar)
end

%% Return
if count >= n
    disp("Fail to obtain a solution satisfying the given conditions.")
end
x = x1;

end

%%
function display_detail(count,p,x1,x0,xstar)
if p
    fprintf('| %d', count)
    for i = 1:size(x1,1)
        fprintf(' | %0.9f', x1(i,1))
    end
    if all(xstar==0)
        fprintf(' | %0.9f |\n', max(abs(x1-x0)))
    else
        fprintf(' | %0.9f |\n', max(abs(x1-xstar)))
    end
end
end