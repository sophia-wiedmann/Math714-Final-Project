function F = buildF(z)
% Constructs matrix F representing discretized differential equation
% Currently F represents grad(D(x)*grad(C)) = grad(D)*grad(C) + D*grad^2(C)
% Input: z measurement of slice (in mm)

% MRI data
global greyVol
global whiteVol
global xdim
global ydim
global h
global b

% Constants
zStart = -72;
zval = z-zStart;
numPoints = xdim*ydim;

% Construct matrices
D = diffusionCoeff(greyVol,whiteVol,zval);
D = reshape(D,numPoints,1);
D_diag = spdiags(D,0,numPoints,numPoints);

D1 = buildGradient(z);
D2 = buildLaplacian(z);

gradD = D1*D;
gradD_diag = spdiags(gradD,0,numPoints,numPoints);

F = b/h^2*D2;
%F = D2*D_diag; % treats D like a constant
%F = D1*gradD_diag + D2*D_diag; % treats D like D(x)

end