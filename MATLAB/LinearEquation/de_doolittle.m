function [L,U,P] = de_doolittle(A)
% LU decomposition: Doolittle decompolition, A = LU, 
% A is an invertible matrix.
% L is an unit lower triangular matrix.
% U is an upper triangle matrix.
%
% Input -
% A: NxN matrix
%
% Output -
% L: NxN, unit lower triangular matrix
% U: NxN, upper triangular matrix
% P: NxN, permutation matrix, record of lines swaping, PA=LU
% 
% Example -
% >> [L,U] = de_doolittle(A)
% >> [L,U,P] = de_doolittle(A)
%
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

%% Initialization
L = eye(m);
U = zeros(m,n);
P = eye(m);

%% LU matrix calculation
% check first line of A, ensure that a_{11} = max(A(:,1))
[~,index] = max(abs(A(1:m,1)));
if index ~= 1
    temp = A(1,:);
    A(1,:) = A(index,:);
    A(index,:) = temp;
    temp = P(1,:);
    P(1,:) = P(index,:);
    P(index,:) = temp;
end

% first column of A ==> first column of L
A(2:m,1) = A(2:m,1)/A(1,1);

for r = 2:m
    %% select pivot
    A(r:m,r) = A(r:m,r) - A(r:m,1:(r-1))*A(1:(r-1),r);
    [~,index] = max(abs(A(r:m,r)));
    index = index + r - 1;
    
    %% swap
    if index ~= r
        temp = A(r,:);
        A(r,:) = A(index,:);
        A(index,:) = temp;
        temp = P(r,:);
        P(r,:) = P(index,:);
        P(index,:) = temp;
    end
    
    %% rth line of A ==> rth line of U ((r+1):n)
    A(r,(r+1):n) = A(r,(r+1):n) - A(r,1:(r-1))*A(1:(r-1),(r+1):n);
    
    %% rth column of A ==> rth column of L ((r+1):m)
    A((r+1):m,r) = A((r+1):m,r) / A(r,r);
end

%% return
for r = 1:m
    U(r,r:n) = A(r,r:n);
    L((r+1):m,r) = A((r+1):m,r);
end

end
