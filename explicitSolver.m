function C_np1 = explicitSolver(z,rho,C,D,gradD)
% Solve differential equation for one time step using explicit methods
% Inputs: 
    % z is measurement of slice (in mm)
    % C is 2D matrix of concentrations at points (x,y) at time t = n and
    % fixed z value (in mm)
    % D is diffusion coefficient matrix
    % gradD is matrix representing discretization of grad(D(x))
% Output: C_nplus1 is concentration at time t = n+1

global h
global k
global xdim
global ydim
% global zdim

% numPoints = xdim*ydim*zdim;
numPoints = xdim*ydim;
sz = [xdim ydim];

% Constants
zStart = -72;
zval = z-zStart;

% Initialize solution vector
C_np1 = zeros(numPoints,1);

% Update each voxel using explicit method
for ijk = 1:numPoints
    % Get (x,y,z) coordinates for point
    %[i, j, k] = ind2sub(sz, ijk);
    [i, j] = ind2sub(sz, ijk);

    % Indices of neighbors
    if i+1 <= xdim
        %i1jk = sub2ind(sz, i+1, j, k);
        i1jk = sub2ind(sz, i+1, j);
    else
        i1jk = 0;
    end

    if i-1 >= 1
        %im1jk = sub2ind(sz, i-1, j, k);
        im1jk = sub2ind(sz, i-1, j);
    else
        im1jk = 0;
    end

    if j+1 <= ydim
        %ij1k = sub2ind(sz, i, j+1, k);
        ij1k = sub2ind(sz, i, j+1);
    else
        ij1k = 0;
    end

    if j-1 >= 1
        %ijm1k = sub2ind(sz, i, j-1, k);
        ijm1k = sub2ind(sz, i, j-1);
    else
        ijm1k = 0;
    end

    % if k+1 <= zdim
    %     ijk1 = sub2ind(sz, i, j, k+1);
    % else
    %     ijk1 = 0;
    % end
    % 
    % if k-1 >= 1
    %     ijkm1 = sub2ind(sz, i, j, k-1);
    % else
    %     ijkm1 = 0;
    % end

    %neighbors = [i1jk im1jk ij1k ijm1k ijk1 ijkm1];
    neighbors = [i1jk im1jk ij1k ijm1k];

    % Check if orthogonal neighbors are in the boundary
    for index = 1:length(neighbors)
        neighbor = neighbors(index);
        if neighbor ~= 0 % Neighbor is in grid
            % [nbr_x, nbr_y, nbr_z] = ind2sub(sz,neighbor);
            % point = [nbr_x, nbr_y, nbr_z];
            [nbr_x, nbr_y] = ind2sub(sz,neighbor);
            point = [nbr_x, nbr_y, zval];
        else point = [1, 1, 1]; % just assign a random point
        end
        isGhost = (isBoundary(point) | neighbor == 0);
        isNotGhost = (isGhost == 0);

        % Adjust index for ghost nodes (either no change or set = ijk)
        neighbors(index) = isNotGhost*neighbor + isGhost*ijk; 
    end

    % Neighbor indices using ghost nodes
    i1jk = neighbors(1);
    im1jk = neighbors(2);
    ij1k = neighbors(3);
    ijm1k = neighbors(4);

    % D(x) times Laplacian D^2*C
    D2 = D(i1jk)*C(i1jk) + D(im1jk)*C(im1jk) - 4*D(ijk)*C(ijk) ...
        + D(ij1k)*C(ij1k) + D(ijm1k)*C(ijm1k);

    % Product grad(D)*grad(C)
    prod = gradD(i1jk)*C(i1jk) - gradD(im1jk)*C(im1jk) ...
    + gradD(ij1k)*C(ij1k) - gradD(ijm1k)*C(ijm1k);

    % Solution to differential equation
    C_np1(ijk) = k/(2*h)*prod + k/h^2*D2 + k*rho*C(ijk)*(1-C(ijk)) + C(ijk);
end

end