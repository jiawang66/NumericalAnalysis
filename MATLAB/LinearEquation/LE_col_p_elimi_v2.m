function [x, A1, b1, P] = LE_col_p_elimi_v2(A,b)
% Solve linear equation Ax=b using column pivot elimination 
% with line swap to obtain the maximum pivot.
% 
% A need to be a invertible matrix
% 
% Input -
% A: NxN square matrix
% b: NxS array
%
% Output -
% x: NxS array, solution
% A: NxN matrix, eliminated A
% b: NxS array, eliminated b
% P: NxN matrix, permutation matrix, record of lines swaping, A_{out}=PA_{in}
%
% Example -
% >> x = LE_col_elimi_v2(A,b)
% >> [x, A1, b1] = LE_col_elimi_v2(A,b)

%% Parameter check
if nargin ~= 2
   error('Error! Lack of input parameter, must be equal to 2.') 
end

if ~ismatrix(A) || length(size(A)) ~= 2
   error('Error! Inputed parameter should be a 2-D matrix.') 
end

if ~ismatrix(b)
   error('Error! Input parameter b must be a matrix.') 
end

[m,n] = size(A);
[m1,~] = size(b);

if m~=n || rank(A)~=m
    error('Errot! Inputed matrix should be invertible.')
end

if m ~= m1
   error('Error! The length of 1st dimension of A and b must be equal.') 
end

%% Column pivot elimination
[A,b,P] = de_elimination(A,b);

%% back substitution
x = substitution_backward(A,b);
A1 = A;
b1 = b;

end