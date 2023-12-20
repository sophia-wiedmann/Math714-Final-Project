function test = isBoundary(point)
% Checks if a gridpoint is in the anatomic boundary
% Inputs: point = [x, y, z] is gridpoint being checked, skullVol and
    % backgroundVol contain locations of boundary as determined by readData
% Output: test = 1 if it's a boundary gridpoint (0 if not)

% Dimensions of MRI data
global xdim
global ydim
global zdim

% MRI data
global skullVol
global backgroundVol
global CSFVol

% Set threshold levels for determining boundaries
skull_threshold = 0.1;
background_threshold = 0.1;
CSF_threshold = 0.1;

% Coordinates of gridpoint
x = point(1);
y = point(2);
z = point(3);

% Check if it's in the boundary
if (x > xdim | y > ydim | z > zdim)
    test = 1;
elseif skullVol(x,y,z) >= skull_threshold
    test = 1;
elseif backgroundVol(x,y,z) >= background_threshold
    test = 1;
elseif CSFVol(x,y,z) >= CSF_threshold
    test = 1;
else
    test = 0;
end

end