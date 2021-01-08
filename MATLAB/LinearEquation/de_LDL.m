function [L,D] = de_LDL(A)
%
% function [L,D] = de_LDL(A)
% 
% LDL decomposition of symmetrical positive determined matrix,
%   A = LDL^{T}
% Warning! In this programm, whether A is a positive determined matrix or
% not will not be performed.
%
% Input -
% A: NxN matrix, a symmetrical positive determined matrix
%
% Output -
% L: NxN matrix, unit lower triangular matrix
% D: NxN matrix, diagonal matrix

%% Parameter check
if nargin ~= 1
    error('Error! One parameter needed.')
end

if ~ismatrix(A)
    error('Error! Inputed parameter is not a matrix.')
end

if transpose(A) ~= A
    error('Error! A symmetrical positive determined matrix needed.')
end

%% Initialization
[m,~] = size(A);
L = eye(m);
d = zeros(1,m); % pivots of D

%% LDL decomposition
if A(1,1) == 0
    error('Error! Pivot cannot be zero.')
else
    d(1) = A(1,1);
end

for i = 2:m
    for j = 1:i-1
        if j == 1
            L(i,j) = A(i,j) / d(1);
        elseif d(j) ~= 0
            L(i,j) = (A(i,j) - sum(L(i,1:j-1).*d(1:j-1).*L(j,1:j-1))) / d(j);
        else
            error('Error! denominator can not be zero.')
        end
    end
    d(i) = A(i,i) - sum(L(i,1:i-1).*L(i,1:i-1).*d(1:i-1));
end

%% Return
D = diag(d);

end
