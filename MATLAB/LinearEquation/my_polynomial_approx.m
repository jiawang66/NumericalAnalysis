function [approx, coeff, base] = my_polynomial_approx(f, n, a, b)
% 
% function approx = my_polynomial_approx(f, n, a, b)
% 
% n-order polynomial approximation of real function in C[a,b].
%
% Input -
% f: symbolic function, define in real number field
% n: the highest order of polynomial used to approximation
%
% Output -
% approx: the polynomial approximation of inputed function f
% base: cell, the unitized orthogonal base of n-order polynomials
% coeff: coefficient of the polynomial approximation on the base

%% Parameter check
if nargin ~= 4
    error('Error! Lack of inputed parameters. Four parameters needed.')
end

if class(f) ~= "sym"
    error('The class of inputed function should be "sym".')
end

%% calculate the base of n-order orthogonal polynomials in C[a,b]
base = my_ortho_polynomials(n,a,b);

%% calculate the polynomial approximation of function f
x = sym('x');
approx = sym('0');
coeff = zeros(n+1, 1);

for i = 1:(n+1)
    coeff(i) = int(f * base{i}, x, a, b);
    approx = approx + coeff(i) * base{i};
end

%% display
% base
disp('----- base of  n-order polynomials -----');
for i = 1:(n+1)
    disp(vpa(base{i},4))
end

% plot function
figure;
fplot(approx), hold on
fplot(f, [a,b])

end