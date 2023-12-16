function D2 = buildLaplacian(z)
% Builds square matrix D2 representing discretization of Laplacian in 2D at
% given z value (in mm)

% Constants
global xdim
global ydim
% global zdim
global h

hfac = 1/(2*h^2);
zStart = -72;
zval = z-zStart;

%numPoints = xdim*ydim*zdim;
%sz = [xdim, ydim, zdim];
numPoints = xdim*ydim;
sz = [xdim, ydim];

% Initialize D2 as sparse matrix
%numNonzero = 7*numPoints;
numNonzero = 5*numPoints;
D2 = spalloc(numPoints,numPoints,numNonzero);

%for ijk = 10000:11000
for ijk = 1:numPoints
    % Get (x,y,z) coordinates for point
    %[i, j, k] = ind2sub(sz, ijk);
    [i, j] = ind2sub(sz, ijk);

    % Update diagonal entry of D2
    %D2(ijk,ijk) = -6*hfac;
    D2(ijk,ijk) = -4*hfac;

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

    % Check if orthogonal neighbors are in the boundary and update F
    for index = 1:length(neighbors)
        neighbor = neighbors(index);
        if neighbor ~= 0 % Neighbor is in grid
            % [nbr_x, nbr_y, nbr_z] = ind2sub(sz,neighbor);
            % point = [nbr_x, nbr_y, nbr_z];
            [nbr_x, nbr_y] = ind2sub(sz,neighbor);
            point = [nbr_x, nbr_y, zval];
        else % Neighbor is outside grid
            point = [1, 1, 1]; % just assign a random point
        end
        isGhost = (isBoundary(point) | neighbor == 0);
        isNotGhost = (isGhost == 0);

        % Adjust index for ghost nodes (either no change or set = ijk)
        neighbor = isNotGhost*neighbor + isGhost*ijk; 

        % Update row ijk of D2
        D2(ijk,neighbor) = D2(ijk,neighbor) + 1*hfac;
    end
end

end