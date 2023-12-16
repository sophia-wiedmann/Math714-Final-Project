
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

xdim = 181;
ydim = 217;

% Tissue-specific diffusion coefficients from Table 11.6 in textbook
    % Units are cm^2/day
    % Tumor grading is high (HH), intermediate (HL), intermediate (LH),
    % and low (LL)
Dg = 1.3*10^(-3); % HH
%Dg = 1.3*10^(-4); % HL
%Dg = 1.3*10^(-3); % LH
%Dg = 1.3*10^(-4); % LL
Dw = 5*Dg; % max diffusion coefficient

h = 0.1; % h = 1mm = 0.1cm
k = 1/ceil(1/(h^2/(6*Dw))); % choose k <= h^2/(6*Dw)

% Read in data
readData

% Constants
z = 28; % this can be changed to another slice (in mm)
zStart = -72;
zval = z-zStart;
numPoints = xdim*ydim;

% Matrix for spatial discretization
F = buildF(z);

%% Initial Condition (normal distribution -- see 11.9 in book)
x0 = [111, 50, zval]; % center of tumor
a = 1; % max density at center of tumor
r = 3; % radius of tumor in mm
cutoff = 0.01; % density at radius r
b = -r^2/log(cutoff/a); % measure of spread so that cutoff condition is satisfied

% Initialize IC
IC = zeros(xdim,ydim);

% Compute IC at each grid point
for x = 1:xdim
    for y = 1:xdim
        dist2 = (x-x0(1))^2 + (y-x0(2))^2; % squared distance to center of tumor
        IC(x,y) = a*exp(-dist2/b);
    end
end

%% Simulate tumor growth
% Initialize concentration vector
C_n = reshape(IC,numPoints,1);

% Simulate for 100 time steps
for t = 1:100
    C = C_n;
    C_n = solver(C,F);
end


%% Plot results
% Reshape to matrix
C_n = reshape(C_n,xdim,ydim);

% Relevant data for plotting
greyData = greyVol(:,:,zval);
whiteData = whiteVol(:,:,zval);
skullData = skullVol(:,:,zval);
greyData = greyData';
whiteData = whiteData';
skullData = skullData';

% Plots
figure;
s = pcolor(greyData);
s.FaceColor = 'interp';
colorbar;
title("Grey Matter");
axis image

figure;
s = pcolor(whiteData);
s.FaceColor = 'interp';
colorbar;
title("White Matter");
axis image

figure;
s = pcolor(skullData);
s.FaceColor = 'interp';
colorbar;
title("Skull");
axis image

figure;
s = pcolor(IC');
s.FaceColor = 'interp';
colorbar;
title("Initial Condition");
axis image

figure;
s = pcolor(C_n');
s.FaceColor = 'interp';
colorbar;
title("Tumor");
axis image

figure;
s = pcolor(greyData + C_n');
s.FaceColor = 'interp';
colorbar;
title("Grey Matter and Tumor");
axis image