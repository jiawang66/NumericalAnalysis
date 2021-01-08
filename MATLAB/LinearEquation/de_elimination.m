function [A,b,P] = de_elimination(A,b)
% 
% function [A,b,P] = de_elimination(A,b)
% 
% Column pivot elimination with swaping lines, but without back substitution
%
% Input -
% A: NxN matrix
% b: NxS matrix
%
% Output -
% A: NxN, matrix after elimination of A
% b: NxS, matrix after elimination of b
% P: NxN, permutation matrix, record of lines swaping, A_{out}=PA_{in}
%
% Example -
% >> A1 = de_elimination(A);
% >> A1 = de_elimination(A,b);
% >> [A1,~,P] = de_elimination(A,b);
% >> [A1,b1,P] = de_elimination(A,b);

%% Parameter check
if nargin == 1
    flag = 1;
elseif nargin == 2
    flag = 2;
else
    error('Error! The number of input should be 1 or 2.')
end

if ~ismatrix(A)
   error('Error! The first parameter is not a matrix.') 
end

[m,n] = size(A);
if flag == 2
    if ~ismatrix(b)
        error('Error! The second parameter is not a matrix.')
    end
    
    [m1,s] = size(b);
    if m ~= m1
        error('Error! The length of 1st dimension of A and b must be equal.')
    end
end

%% Elimination
P = eye(m);
for k = 1:(m-1)
    %% select pivot
    [a_max,index] = max(abs(A(k:m,k)));
    index = index + k-1;
    if a_max == 0
        continue
    end
    
    %% swap
    if index ~= k
        temp = A(k,:);
        A(k,:) = A(index,:);
        A(index,:) = temp;
        temp = P(k,:);
        P(k,:) = P(index,:);
        P(index,:) = temp;
        if flag == 2
            temp = b(k,:);
            b(k,:) = b(index,:);
            b(index,:) = temp;
        end
    end
    
    %% eliminition
    multiplier = A((k+1):m,k)/A(k,k);
    A((k+1):m,k:n) = A((k+1):m,k:n) ...
        - repmat(multiplier,1,n-k+1) .* repmat(A(k,k:n),m-k,1);
    if flag == 2
        b((k+1):m,:) = b((k+1):m,:) ...
            - repmat(multiplier,1,s) .* repmat(b(k,:),m-k,1);
    end
end

for i = 2:m
   A(i,1:(i-1)) = zeros(1,i-1);
end

end