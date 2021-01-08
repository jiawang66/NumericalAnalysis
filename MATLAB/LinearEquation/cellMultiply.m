function z = cellMultiply(x, y)
%
% function z = cellMultiply(x, y)
%
% Multiplication of cell matrix (block matrix)
%
% Input --
% x: cell matrix, MxT
% y: cell matrix, TxN
%
% Output --
% z: matrix, MxN
%
% Example --
% >> z = cellMultiply(x, y)

%% Parameters check
if size(x,2) ~= size(y,1)
    error('Error! The #cols of x must be equal to the #rows of y.')
end

%% Initialization
[m, t] = size(x);
[~, n] = size(y);
z = cell(m,n);

%% Compute
for i = 1 : m
    for j = 1 : n
        z{i,j} = zeros(size(x{i,1},1),size(y{1,j},2));
        for k = 1 : t
            z{i,j} = z{i,j} + x{i,k} * y{k,j};
        end
    end
end

%% Return
z = cell2mat(z);

end