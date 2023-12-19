function F = buildF(z)
% Constructs matrix F representing discretized differential equation
% Currently F represents grad(D(x)*grad(C)) = grad(D)*grad(C) + D*grad^2(C)
% Input: z measurement of slice (in mm)

% Constants
global xdim
global ydim
global h
global b
numPoints = xdim*ydim;

% Growth rate [1/day] from Table 11.6 in textbook
% Tumor grading is high (HH), intermediate (HL), intermediate (LH),
    % and low (LL)
rho = 1.2*10^(-2); % HH
%rho = 1.2*10^(-2); % HL
%rho = 1.2*10^(-3); % LH
%rho = 1.2*10^(-3); % LL

% Tissue-specific diffusion coefficients from Table 11.6 in textbook
    % Units are cm^2/day
    % Tumor grading is high (HH), intermediate (HL), intermediate (LH),
    % and low (LL)
Dg = 1.3*10^(-3); % HH
%Dg = 1.3*10^(-4); % HL
%Dg = 1.3*10^(-3); % LH
%Dg = 1.3*10^(-4); % LL

% Diffusion coefficient D(x)
D = diffusionCoeff(z,Dg);
D = reshape(D,numPoints,1);
D_diag = spdiags(D,0,numPoints,numPoints);

% Discretizations of derivatives
D1 = buildGradient(z); 
D2 = buildLaplacian(z);

gradD = 1/(2*h)*D1*D; % graph looks good, range is about +/- 0.15
% Graph to test gradD
 gradD = reshape(gradD,xdim,ydim);
 figure;
 s = pcolor(gradD');
 s.FaceColor = 'interp';
 colorbar;
 title("gradD");
 axis image
size(gradD)
numPoints

gradD = reshape(gradD, numPoints, 1);
gradD_diag = spdiags(gradD,0,numPoints,numPoints);

%F = b/h^2*D2; % constant diffusion
%F = 1/h^2*D2*D_diag; % treats D like a constant, but depends on tissue type
F = 1/(2*h)*D1*gradD_diag + 1/h^2*D2*D_diag; % treats D like D(x)

end