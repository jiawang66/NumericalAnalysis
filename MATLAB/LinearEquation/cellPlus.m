function z = cellPlus(x, y, f)
%
% function z = cellPlus(x, y, f)
%
% Compute the addition or subtraction of cell matrix
%
% Input -
% x: NxM cell matrix
% y: NxM cell matrix
% f: 1-addition, (-1)-subtraction, default f = 1
%
% Output -
% z: NxM cell matrix
%
% Example -
% >> z = cellPlus(x, y)
% >> z = cellPlus(x, y, 1)

%% Parameter check
if ~iscell(x) || ~iscell(y)
    error('Error! Type of inputs needed to be cell matrix.')
end

if size(x) ~= size(y)
    error('Error! Size of x and y needed to be equal.')
end

if nargin < 3
    f = 1;
elseif f ~= 1 && f ~= -1
    error('Error! f = 1 or -1.')
end

%% Initialization
[m, n] = size(x);
z = cell(m,n);

%% Compute
for i = 1 : m
    for j = 1 : n
        z{i,j} = x{i,j} + f * y{i,j};
    end
end

end