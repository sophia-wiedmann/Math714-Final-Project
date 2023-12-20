function F = buildF(z,type)
% Constructs matrix F representing discretized differential equation
% Input: z measurement of slice (in mm) and type = 1,2,3,4 determines
% grading

% Constants
global xdim
global ydim
global h
numPoints = xdim*ydim;

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

% Diffusion coefficient D(x)
D = diffusionCoeff(z,Dg);
D = reshape(D,numPoints,1);
D_diag = spdiags(D,0,numPoints,numPoints);

% Discretizations of derivatives
D1 = buildGradient(z); 
D2 = buildLaplacian(z);

gradD = 1/(2*h)*D1*D;
gradD = reshape(gradD, numPoints, 1);
gradD_diag = spdiags(gradD,0,numPoints,numPoints);

F = 1/(2*h)*D1*gradD_diag + 1/h^2*D2*D_diag;

end