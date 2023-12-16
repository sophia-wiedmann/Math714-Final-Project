function D = diffusionCoeff(greyVol,whiteVol,z)
% Computes diffusion coefficient matrix D (2D matrix)

% Inputs: greyVol and whiteVol are 3D matrices describing what percentage of
    % each cell is grey (white) matter in range [0,1]

% Output: 3D matrix D where D(x,y,z) is diffusion coefficient at point
    % (x,y,z) computed based on percentage of grey and white matter (0
    % everywhere else)

% MRI data dimensions
global xdim
global ydim
% global zdim

zStart = -72;
zval = z-zStart;

% Tissue-specific diffusion coefficients (need to update)
D_grey = 1;
D_white = 1;

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
            D(x,y) = greyVal*D_grey + whiteVal*D_white;
        % end
    end
end

end