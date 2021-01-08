function c = mycond(A, v)
%
% function c = mycond(A, v)
%
% Compute condition number of matrix A, A should be invertible.
% Cond(A)_v = norm(A^{-1})_v * norm(A)_v.
%
% Input -
% A: NxN matrix
% v: 0 -- L_infinity norm
%    1 -- L1 norm
%    2 -- L2 norm
%
% Output -
% c: number, condition number of A

%% Parameter check
if nargin ~= 2
    error('Error! Two parameter needed.')
end

if ~ismatrix(A) || length(size(A)) ~= 2
    error('Error! First parameter should be a 2-D numerical matrix.')
end

[m,n] = size(A);
if m~=n || rank(A)~=m
    error('Error! First parameter should be an invertible matrix.')
end

%% Compute condition number
if v==0 || v==1 || v==2
    c = matnorm(inv(A),v) * matnorm(A,v);
else
    error('Error! Second inputed parameter should be 0, 1, or 2.')
end

end