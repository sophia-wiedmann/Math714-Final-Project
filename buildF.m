function F = buildF(z)
% Constructs matrix F representing discretized differential equation
% Currently F represents grad(D(x)*grad(C)) = grad(D)*grad(C) + D*grad^2(C)
% Input: z measurement of slice (in mm)

% MRI data
global greyVol
global whiteVol
global xdim
global ydim

% Constants
zStart = -72;
zval = z-zStart;
numPoints = xdim*ydim;

% Construct matrices
D = diffusionCoeff(greyVol,whiteVol,zval);
D = reshape(D,numPoints,1);
D_diag = spdiags(D,0,numPoints,numPoints);

D1 = buildGradient(zval);
D2 = buildLaplacian(zval);

gradD = D1*D;
gradD_diag = spdiags(gradD,0,numPoints,numPoints);

%F = D2;
%F = D2*D_diag; % treats D like a constant
F = D1*gradD_diag + D2*D_diag; % treats D like D(x)

end