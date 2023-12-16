function D = diffusionCoeff(z,Dg)
% Computes diffusion coefficient matrix D (2D matrix)

% Inputs: z is measurement of slice in mm, Dg is diffusion coefficent for
% grey matter

% Output: 2D matrix D where D(x,y,z) is diffusion coefficient at point
    % (x,y,z) computed based on percentage of grey and white matter (0
    % everywhere else)

% MRI data
global greyVol
global whiteVol

% MRI data dimensions
global xdim
global ydim
% global zdim

zStart = -72;
zval = z-zStart;

% Tissue-specific diffusion coefficients from Table 11.6 in textbook
    % Units are cm^2/day
    % Tumor grading is high (HH), intermediate (HL), intermediate (LH),
    % and low (LL)
Dw = 5*Dg;

% Initialize diffusion coefficient matrix as square matrix using max
% dimension
% D = zeros(xdim, ydim, zdim);
D = zeros(xdim, ydim);

% Compute diffusion coefficent for each cell
for x = 1:xdim
    for y = 1:ydim
         % for z = 1:zdim
            greyVal = greyVol(x,y,zval);
            whiteVal = whiteVol(x,y,zval);
            D(x,y) = greyVal*Dg + whiteVal*Dw;
        % end
    end
end

end