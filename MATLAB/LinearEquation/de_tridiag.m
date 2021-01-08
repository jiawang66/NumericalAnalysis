function [L,U] = de_tridiag(A)
% LU decomposition of tridiagonal matrix
%
% Input -
% A: NxN matrix, a tridiagonal matrix
%
% Output -
% L: NxN matrix, a lower triangular matrix
% U: NxN matrix, an unit upper triangular matrix
% 
% Example -
% >> [L,U] = de_tridiag(A)
%% Parameter check
if nargin ~= 1
    error('Error! The number of inputed parameter should be 1.')
end

if ~ismatrix(A) || length(size(A)) ~= 2
   error('Error! Inputed parameter should be a 2-D matrix.') 
end

[m,n] = size(A);

if m~=n || rank(A)~=m
    error('Errot! Inputed matrix should be invertible.')
end

% check whether A is a tridiagonal matrix or not
for i = 1:m
    for j = 1:n
       if A(i,j)~=0 && abs(i-j)>1
          error('Error! Inputed matrix is not a tridiagonal matrix.') 
       end
    end
end

%% Calculation of L and U
A(1,2) = A(1,2) / A(1,1);
for i = 2:(m-1)
   A(i,i) = A(i,i) - A(i,i-1)*A(i-1,i); 
   if A(i,i) ~= 0
       A(i,i+1) = A(i,i+1) / A(i,i);
   else
       error('Error! Denominator cannot be zero.')
   end
end
A(m,m) = A(m,m) - A(m,m-1)*A(m-1,m);

%% Return
L = zeros(m,m);
U = eye(m);
L(1,1) = A(1,1);
L(m,(m-1):m) = A(m,(m-1):m);
U(1,2) = A(1,2);

for i = 2:(m-1)
    L(i,(i-1):i) = A(i,(i-1):i);
    U(i,i+1) = A(i,i+1);
end

end