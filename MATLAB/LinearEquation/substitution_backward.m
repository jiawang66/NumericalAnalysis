function x = substitution_backward(A, b)
% back substitution of equation Ax=b, A is an invertible upper triangular matrix.
%
% Input -
% A: NxN matrix, invertible and upper triangular
% b: NxS matrix
% 
% Ouput -
% x: NxS matrix, solution of Ax=b
% 
% Example -
% >> x = substitution_backward(A,b)

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
for i = 2:m
   if A(i,1:(i-1)) ~= zeros(1,i-1)
      error('Error! First inputed parameter is not an upper triangular matrix.') 
   end
end

%% back substitution, from down to up
x = zeros(m,s);

if A(m,m) ~= 0
    x(m,:) = b(m,:) / A(m,m);
else
    error('Error! The pivot cannot be zero.')
end

for i = fliplr(1:(m-1))
    if A(i,i) ~= 0
        x(i,:) = b(i,:) - sum(x((i+1):m,:).*repmat(A(i,(i+1):m)',1,s), 1);
        x(i,:) = x(i,:) / A(i,i);
    else
        error('Error! The pivot cannot be zero.')
    end
end

end