% For testing
N = 16;
h = 1/(N+1);
bd = [0,1];
u_0y = @(x,y)1;
u_1y = @(x,y)exp(y);
u_x0 = @(x,y)1;
u_x1 = @(x,y)exp(x);
[A,d] = my_poisson_matrix(h,bd,u_0y,u_1y,u_x0,u_x1);

% LU
