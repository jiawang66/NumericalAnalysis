function y = cellAbs(x)
%
% function y = cellAbs(x)
%
% the abs function of the cell matrix
%
% Input -
% x: cell matrix, MxN
%
% Output -
% y: cell matrix, MxN of absoulte value

%% Parameters check
if ~iscell(x)
    error('Error! Cell vector input needed.')
end

%% Compute
[m, n] = size(x);
y = cell(m,n);

for i = 1 : m
    for j = 1 : n
        y{i,j} = abs(x{i,j});
    end
end

end