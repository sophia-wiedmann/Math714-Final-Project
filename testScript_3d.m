
%% Data and model parameters
% Dimensions of MRI  data
global xdim
global ydim
global zdim

% MRI data
global skullVol
global backgroundVol
global greyVol
global whiteVol

% Constants
global h
%global b
global k
global rho

xdim = 181;
ydim = 217;
zdim = 181;

% Tissue-specific diffusion coefficients from Table 11.6 in textbook
    % Units are cm^2/day
    % Tumor grading is high (HH), intermediate (HL), intermediate (LH),
    % and low (LL)
Dg = 1.3*10^(-3); % HH
%Dg = 1.3*10^(-4); % HL
%Dg = 1.3*10^(-3); % LH
%Dg = 1.3*10^(-4); % LL
Dw = 5*Dg; % max diffusion coefficient

rho = 1.2*10^(-2); % HH
%rho = 1.2*10^(-2); % HL
%rho = 1.2*10^(-3); % LH
%rho = 1.2*10^(-3); % LL

h = 0.1; % h = 1mm = 0.1cm
k = 1/ceil(1/(h^2/(6*Dw))); % choose k <= h^2/(6*Dw)

% Read in data
readData

% Constants
z = 28; % this can be changed to another slice (in mm)

numPoints = xdim*ydim*zdim;

% Matrix for spatial discretization
load('F.mat');

%% Initial Condition (normal distribution -- see 11.9 in book)
x0 = [111, 50, 111]; % center of tumor
a = 1; % max density at center of tumor
r = 3; % radius of tumor in mm
cutoff = 0.01; % density at radius r
b = -r^2/log(cutoff/a); % measure of spread so that cutoff condition is satisfied

% Initialize IC
IC = zeros(xdim,ydim,zdim);

% Compute IC at each grid point
for x = 1:xdim
    for y = 1:xdim
        for z = 1:zdim
            dist2 = (x-x0(1))^2 + (y-x0(2))^2 + (z - x0(3))^2; % squared distance to center of tumor
            IC(x,y, z) = a*exp(-dist2/b);
            
        end
    end
end

%% Simulate tumor growth
C_n = reshape(IC,numPoints,1);

% Simulate for 100 time steps
for t = 1:100
    C = C_n;
    C_n = solver_3d(C,F);
end

%% Plot results
% Grab just the data where the tumor concentration is > 10^-6
X = [];
Y = [];
Z = [];
val = [];

% Reshape to matrix
C_n = reshape(C_n,xdim,ydim, zdim);
for x = 1:xdim
    for y = 1:ydim
        for z = 1:zdim
            if C_n(x, y, z)> 10^(-6)
                X(end + 1) = x;
                Y(end + 1) = y;
                Z(end + 1) = z;
                val(end+1) = C_n(x, y, z);
            end
        end
    end
end





size(C_n)
w = (C_n + .01)*100;
% Relevant data for plotting
scatter3(X, Y, Z, val*100 ,'filled')
set(gca,'XLim',[0 181],'YLim',[0 218],'ZLim',[0 181])

