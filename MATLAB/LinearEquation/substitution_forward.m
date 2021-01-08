function x = substitution_forward(A,b)
% back substitution of equation Ax=b, A is an invertible lower triangular matrix.
%
% Input -
% A: NxN matrix, invertible and lower triangular
% b: NxS matrix
% 
% Ouput -
% x: NxS matrix, solution of Ax=b
% 
% Example -
% >> x = substitution_forward(A,b)

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
[m1,s] = size(b);

if n ~= m
    error('Error! Only deal with square matrix.')
end

if m ~= m1
   error('Error! The length of 1st dimension of A and b must be equal.') 
end

% TO DO: calculate rank(A) by myself
if rank(A) ~= m
   error('Error! Input parameter A must be a invertible matrix.') 
end

% check A whether is an upper triangular matrix
for i = 1:(m-1)
   if A(i,(i+1):m) ~= zeros(1,m-i)
      error('Error! First inputed parameter is not a lower triangular matrix.') 
   end
end

%% back substitution, from down to up
x = zeros(m,s);

if A(1,1) ~= 0
    x(1,:) = b(1,:) / A(1,1);
else
    error('Error! The pivot cannot be zero.')
end

for i = 2:m
    if A(i,i) ~= 0
        x(i,:) = b(i,:) - sum(x(1:(i-1),:).*repmat(A(i,1:(i-1))',1,s), 1);
        x(i,:) = x(i,:) / A(i,i);
    else
        error('Error! The pivot cannot be zero.')
    end
end

end