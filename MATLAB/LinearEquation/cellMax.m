function y = cellMax(x)
%
% function y = cellMax(x)
%
% Compute the maximal of cell vector x
%
% Input -
% x: NxM cell matrix
%
% Output -
% y: scale number
%
% Example -
% >> y = cellMax(x)

%% Parameter check
if ~iscell(x)
    error('Error! Cell input needed.')
end

%% Initialization
[m, n] = size(x);
y = zeros(m,n);

%% Compute
for i = 1 : m
    for j = 1 : n
        y(i,j) = max(x{i,j}(:));
    end
end

%% return
y = max(y(:));
    
end