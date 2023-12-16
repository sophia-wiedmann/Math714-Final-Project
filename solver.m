function C_nplus1 = solver(C,F)
% Solve differential equation for one time step
% Inputs: 
    % C is 2D matrix of concentrations at points (x,y,z) at time t = n and
    % fixed z value
    % F is matrix resulting from discretization of spatial components of
    % differential equation
% Output: C_nplus1 is concentration at time t = n+1

global h
global k
global xdim
global ydim
% global zdim

% Turn C into vector
% numPoints = xdim*ydim*zdim;
% C = reshape(C,numPoints,1,1);
numPoints = xdim*ydim;
C = reshape(C,numPoints,1);

% Construct matrices
I = speye(numPoints);
A = I - k/2*F;
B = I + k/2*F;
BC = B*C;

% Solve A*C_nplus1 = BC
C_nplus1 = A\BC;

% Reshape to (x,y,z) 3D matrix
% C_nplus1 = reshape(C_nplus1,xdim,ydim,zdim);
C_nplus1 = reshape(C_nplus1,xdim,ydim);

end

