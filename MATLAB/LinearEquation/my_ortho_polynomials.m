function base = my_ortho_polynomials(n,a,b)
% 
% function base = my_ortho_polynomials(n)
% 
% calculate the base of n-order orthogonal polynomials on C[a,b]
% 
% Input -
% n: natural number >= 0, the highest order of orthogonal polynomials
% a,b: the range of inner product space, a = -1 and b = 1 by default
%
% Output -
% base: cell, the unitized orthogonal base of n-order polynomials

%% Parameter check
if nargin == 1
    if n ~= ceil(n) || n < 0
        error('Error! The highest order of polynomials should positive integer.')
    end
    a = -1;
    b = 1;
elseif nargin ~= 3
    error('Error! The number of inputed parameters should be 1 or 3.')
elseif a > b
    t = a;
    a = b;
    b = t;
end

%% Initialization
x = sym('x');
base = cell(n+1,1);
L2 = ones(n+1,1);

% choose a simple base
for i = 0:n
    base{i+1} = x^i;
end

%% orthogonalization
L2(1) = int(base{1}*base{1}, x, a, b);

for i = 2:(n+1)
    for j = 1:(i-1)
        p1 = int(base{i}*base{j}, x, a, b);
        base{i} = base{i} - p1 / L2(j) * base{j};
    end
    L2(i) = int(base{i}*base{i}, x, a, b);
end

%% unitization
for i = 1:(n+1)
    base{i} = base{i} / sqrt(L2(i));
end

%% display
% for i = 1:(n+1)
%     disp(base{i})
% end

end
