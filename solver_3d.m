function C_np1 = solver_3d(C,F, t, treat)
% Solve differential equation for one time step
% Inputs: 
    % C is 2D matrix of number of tumerous cells at points (x,y,z) at time t = n
    % F is matrix resulting from discretization of spatial components of
    % differential equation (see buildF.m)
% Output: C_np1 is concentration at time t = n+1

global h
global k
global rho
global xdim
global ydim
global zdim
global rad_kill_rate
global chem_kill_rate


numPoints = xdim*ydim*zdim;

% Construct matrices
I = speye(numPoints);

% Treatment indicatior functions
if treat
    if t >= 0 && t < 42
        R_ind = 1;
        C_ind = 1;
    elseif (t >= 42 && t <= 46) || (t >= 70 && t <= 74) || (t >= 98 && t <= 102) || (t >= 126 && t <= 130) || (t >= 154 && t <= 158) || (t >= 182 && t <= 186)
        C_ind = 2;
        R_ind = 0;
    else
        C_ind = 0;
        R_ind = 0;
    end
else
    C_ind = 0;
    R_ind = 0;
end
    
v = rho*C - rho*C.^2/max(C) - R_ind * rad_kill_rate * C - C_ind * chem_kill_rate * C;
B = I + k*F;
C_np1 = B*C + k*v; % explicity get C^(n+1)
C_np1(C_np1<0) = 0; % don't let number of cells fall below zero
end

