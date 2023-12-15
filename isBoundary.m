function test = isBoundary(point)
% Checks if a gridpoint is in the boundary (skull or background) -- should
% include ventricles?
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

% Set threshold levels for determining boundaries
skull_threshold = 0.1; % need to decide what these values should be
background_threshold = 0.1;

% Coordinates of gridpoint
x = point(1);
y = point(2);
z = point(3);

% Check if it's in the boundary
if (x > xdim | y > ydim | z > zdim)
    test = 1;
elseif (skullVol(x,y,z) >= skull_threshold | backgroundVol(x,y,z) >= background_threshold)
    test = 1;
else
    test = 0;
end

end