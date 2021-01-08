function x = LE_tridiag(A,b)
%
% function x = LE_tridiag(A,b)
% 
% Solve a tridiagonal linear system A*x = b
%
% Input -
% A: NxN matrix, tridiagonal and invertible
% b: NxS matrix
%
% Output -
% x: NxS matrix, solution of A*x = b

%% Parameters check
if nargin ~= 2
   error('Error! Lack of input parameter, must be equal to 2.') 
end

if ~ismatrix(A)
   error('Error! Input parameter A must be a matrix.') 
end

if ~ismatrix(b)
    error('Error! Input parameter b must be a matrix.') 
end

[m,n] = size(A);
[m1,~] = size(b);

if n ~= m
    error('Error! Only deal with square matrix.')
end

if m ~= m1
   error('Error! The length of 1st dimension of A and b must be equal.') 
end

%% LU decomposition
[L, U] = de_tridiag(A);

%% Compute x for A*x = b
% Do forward substitution to solve lower triangular system
b(1,:) = b(1,:) / L(1,1);
for i = 2:m
    b(i,:) = (b(i,:) - L(i,i-1)*b(i-1,:)) / L(i,i);
end

% Do back substitution to solve lower triangular system
for i = (m-1):-1:1
    b(i,:) = b(i,:) - U(i,i+1) * b(i+1,:);
end

%% Set output variable
x = b;

end