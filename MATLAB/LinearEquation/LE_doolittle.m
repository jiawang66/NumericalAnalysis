function [x,y] = LE_doolittle(A,b)
% LU decomposition: Doolittle decompolition, A = LU, 
% A is an invertible matrix.
% Ly = b
% Ux = y
%
% Input -
% A: NxN matrix
% b: NxS matrix
%
% Output -
% x: NxS matrix, solution of Ux=y
% y: NxS matrix, solution of Ly=b
% 
% Example -
% >> x = LE_doolittle(A,b);
% >> [x,y] = LE_doolittle(A,b);

%% Parameter check
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

%% doolittle decompositon
[L,U,P] = de_doolittle(A);

%% back substitution
% solve equation Ly = b
y = substitution_forward(L,P*b);

% solve equation Ux = y
x = substitution_backward(U,y);

end