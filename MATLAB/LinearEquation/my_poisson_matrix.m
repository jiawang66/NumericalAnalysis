function [A,d] = my_poisson_matrix(h,bd,u_0y,u_1y,u_x0,u_x1)
% 
% function [A,d] = my_possion_matrix(h,bx,by,u_0y,u_1y,u_x0,u_x1)
%
% Generating five-point difference scheme of Poisson equation
%   A * u + d = b
%
% Input -
% h: step unit, h = 1/(N+1), default be N=8
% bd: boundary of variable x, y
% u_0y,u_1y,u_x0,x_y0: function handle, the boundary conditions of
% poisson's equation, default be 0
%
% Output -
% A: matrix, five-point difference scheme of Poisson equation
% d: vector, boundary value of five-point difference equation when u_{ij} \in \partial D
%
% Example -
% >> A = my_possion_matrix(h,bx,by)
% >> [A,d] = my_possion_matrix(h,bx,by,u_0y,u_1y,u_x0,u_x1)

%% Parameter check
if h <= 0
    error('Error! Parameter h must be laeger than zero.')
end

if nargin < 3
    error('Error! Lack of inputed parameters.')
elseif nargin == 3
    d = zeros(round((bd(2)-bd(1))/h-1)^2,1);
    flag_d = true;
else
    flag_d = false;
end

%% Initialization
N = round((bd(2)-bd(1))/h - 1);
N2 = N^2;
A = zeros(N2,N2);
E = -eye(N);
dA = diag(ones(N,1)*4) + diag(-ones(N-1,1),1) + diag(-ones(N-1,1),-1);


%% Generate d
if ~flag_d
    d = zeros(N2,1);
    
    %% compute the vector d
    bd1 = ones(1,N-2) * bd(1);
    
    % y = bd(1), left border
    d(1) = - u_x0(bd(1)+h, bd(1)) - u_0y(bd(1), h+bd(1));
    d(2:(N-1)) = - u_x0((2:(N-1))*h+bd1, bd(1));
    d(N) = - u_x0(h*N+bd(1), bd(1)) - u_1y(bd(2), h+bd(1));
    
    % x = bd(1), upper border
    d((N+1):N :(N*(N-2)+1)) = - u_0y(bd(1), (2:(N-1))*h+bd1);
    % x = bd(2), bottom border
    d((2*N):N:(N*(N-1))) = - u_1y(bd(2), (2:(N-1))*h+bd1);
    
    % y = bd(2), right border
    d(N*(N-1)+1) =  - u_x1(h+bd(1), bd(2)) - u_0y(bd(1), h*N+bd(1));
    d((N*(N-1)+2):(N*(N-1)+N-1)) = - u_x1((2:(N-1))*h+bd1, bd(2));
    d(N2) = - u_x1(h*N+bd(1), bd(2)) - u_1y(bd(2), h*N+bd(1));
end

%% Generate A
AC = mat2cell(A, N*ones(N,1), N*ones(1,N));

AC{1,1} = dA;
AC{1,2} = E;
for i = 2:(N-1)
    AC{i,i-1} = E;
    AC{i,i} = dA;
    AC{i,i+1} = E;
end
AC{N,N-1} = E;
AC{N,N} = dA;
A = cell2mat(AC);

end