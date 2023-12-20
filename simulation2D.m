%% Tumor grade
% Grading of tumor (1 = HH, 2 = HL, 3 = LH, 4 = LL)
type = 1;

% Growth rate [1/day] and tissue-specific diffusion coefficients [mm^2/day]
% Tumor grading is high (HH), intermediate (HL), intermediate (LH), and low (LL)
if type == 1
    rho = 1.2*10^(-2); % HH
    Dg = 1.3*10^(-3)*10^2; % HH
elseif type == 2
    rho = 1.2*10^(-2); % HL
    Dg = 1.3*10^(-4)*10^2; % HL
elseif type == 3
    rho = 1.2*10^(-3); % LH
    Dg = 1.3*10^(-3)*10^2; % LH
elseif type == 4
    rho = 1.2*10^(-3); % LL
    Dg = 1.3*10^(-4)*10^2; % LL
else disp("Type must be an integer between 1 and 4")
end


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

% Model constants
global h
global k

% Parameter values
xdim = 181;
ydim = 217;
zdim = 181;
numPoints = xdim*ydim;
Dw = 5*Dg; % diffusion in white matter
h = 0.1; % h = 1mm = 0.1cm
k = 1/ceil(1/(h^2/(6*Dw))); % choose k <= h^2/(6*Dw)

% Read in data
readData

% Constants
z = 28; % this can be changed to another slice (in mm)
zStart = -72;
zval = z-zStart;

% % Diffusion coefficient matrix D(x)
% D = diffusionCoeff(z,Dg);
% D = reshape(D,numPoints,1);
% 
% % Discretizations of D(D(x))
% D1 = buildGradient(z);
% gradD = 1/(2*h)*D1*D;

% Matrix representing discretized differential equation
F = buildF(z,type);

%% Initial Condition (normal distribution -- see 11.9 in book)
x0 = [111, 50, zval]; % center of tumor
a = 5; % max concentration at center of tumor
r = 5; % radius of tumor in mm
cutoff = 1; % concentration at radius r
b = -r^2/log(cutoff/a); % measure of spread so that cutoff condition is satisfied

% Initialize IC
IC = zeros(xdim,ydim);

% Compute IC at each grid point
for x = 1:xdim
    for y = 1:ydim
        % squared distance to center of tumor
        dist2 = (x-x0(1))^2 + (y-x0(2))^2; 

        % Compute initial tumor distribution
        IC(x,y) = a*exp(-dist2/b);
    end
end

%% Simulate tumor growth
% Initialize concentration vector
C_n = reshape(IC,numPoints,1);

% How many days to simulate
days = 100;
timesteps = days/k;

% Simulate for given number of days
tic
for t = 1:timesteps
    C = C_n;
    %C_n = explicitSolver(z,C,D,gradD);
    C_n = solver(C,F);
    
    % Plot results every 30 days
    if mod(t*k,30) == 0
        % Reshape to matrix for plotting
        C_n = reshape(C_n,xdim,ydim);

        % Get time in days
        time = string(t*k);
        text = "Tumor after " + time + " days";
        
        % Plot
        figure;
        s = pcolor(C_n');
        s.FaceColor = 'interp';
        colorbar;
        title(text);
        axis image

        % Reshape to vector to continue solving
        C_n = reshape(C_n,numPoints,1);
    end
end
toc

%% Plot results
% Reshape final solution to matrix
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
s = pcolor(greyData+IC');
s.FaceColor = 'interp';
colorbar;
title("Grey Matter with Initial Condition");
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
s = pcolor(greyData + C_n');
s.FaceColor = 'interp';
colorbar;
title("Grey Matter and Tumor");
axis image