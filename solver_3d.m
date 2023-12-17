function C_np1 = solver(C,F)
% Solve differential equation for one time step
% Inputs: 
    % C is 2D matrix of concentrations at points (x,y) at time t = n and
    % fixed z value (in mm)
    % F is matrix resulting from discretization of spatial components of
    % differential equation (see buildF.m)
% Output: C_nplus1 is concentration at time t = n+1

global h
global k
global rho
global xdim
global ydim
global zdim


numPoints = xdim*ydim*zdim;

% Construct matrices
I = speye(numPoints);
v = rho*C - rho*C.^2;
B = I + k*F  + spdiags(v, 0, numPoints, numPoints );
C_np1 = B*C; %explicity get C^(n+1)


end

