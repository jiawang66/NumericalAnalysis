function [x, A1, b1] = LE_col_p_elimi_v1(A,b)
% Column pivot elimination, solve linear equation Ax=b
% A need to be a invertible matrix
% 
% Input -
% A: NxN square matrix
% b: NxS array
%
% Output -
% x: NxS array, solution
% A1: NxN matrix, eliminated A
% b1: NxS array, eliminated b
%
% Example -
% >> x = LE_col_elimi_v1(A,b)
% >> [x, A1, b1] = LE_col_elimi_v1(A,b)

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

%% Initialization
x = zeros(m,s);
Ab = [A, b];

%% Column pivot elimination
for k = 1:(m-1)
    %% select pivot
    [a_max,index] = max(abs(Ab(k:m,k)));
    index = index + k-1;
    if a_max == 0
        continue
    end
    
    %% swap
    if index ~= k
        temp = Ab(k,:);
        Ab(k,:) = Ab(index,:);
        Ab(index,:) = temp;
    end
    
    %% eliminition
    Ab((k+1):m,:) = Ab((k+1):m,:) ...
        - repmat(Ab((k+1):m,k)/Ab(k,k),1,n+s) .* repmat(Ab(k,:),m-k,1);
end

A1 = Ab(:,1:n);
b1 = Ab(:,(n+1):(n+s));

for i = 2:m
   A1(i,1:(i-1)) = zeros(1,i-1);
end

%% back substitution
if A1(m,m) ~= 0
    x(m,:) = b1(m,:) / A1(m,m);
end

for i = fliplr(1:(m-1))
    x(i,:) = b1(i,:) - sum(x((i+1):m,:).*repmat(A1(i,(i+1):m)',1,s), 1);
    if A1(i,i) ~= 0
        x(i,:) = x(i,:) / A1(i,i);
    end
end

end