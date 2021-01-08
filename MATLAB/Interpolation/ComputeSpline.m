function y = ComputeSpline(S, x)
%
% function y = ComputeSpline(S, x)
%
%   Compute the values based on the given cublic spline interpolation
%   function and array x.
%
% Input -
% S: struct, the return of function CubicSpline
% x: Nx1 vector, nodes to be computed
%
% Output -
% y: Nx1 vector, y_k = S.S(x_k), k = 1,2,...,N

%% Parameters check
if nargin ~= 2
    error('Error! Two parameters needed.')
end

if size(x,2) ~= 1
    error('Error! Second inputed should be a Nx1 vector')
end

%% Initialization
N = numel(x);
y = zeros(N,1);

np = numel(S.x);    % number of interpolation nodes

%% Compute
flag = false;

for k = 1 : N
    for n = 1 : np-1
        if x(k) >= S.x(n) && x(k) <= S.x(n+1)
            y(k) = S.S{n}(x(k));
            flag = true;
            break
        end
    end
    
    if ~flag
        disp(['Warning! Input node x(', num2str(k), ') out of range.'])
        y(k) = nan;
    end
end

end