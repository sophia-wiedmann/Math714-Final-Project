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
global xdim
global ydim
% global zdim

% numPoints = xdim*ydim*zdim;
numPoints = xdim*ydim;

% Construct matrices
I = speye(numPoints);
A = I - k/2*F;
B = I + k/2*F;
BC = B*C;

% Solve A*C_nplus1 = BC
C_np1 = A\BC;

end

