function H = hilmat(n)
%
% function H = hilmat(n)
% 
% Hilbert matrix
%
% Input -
% n: order of Hilbert matrix
%
% Output -
% H: n-order Hilbert matrix

%% Parameter check
if nargin ~= 1
    error('Error! Number of inputed parameter shoule be equal to 1.')
end

if numel(n)~=1 || ~ismatrix(n)
    error('Error! Inputed parameter should be a scalar number.')
end

%% Hilbert matrix
H = zeros(n,n);

for i = 1:n
    for j = 1:n
        H(i,j)= 1 / (i+j-1);
    end
end

end