function D2 = buildLaplacian_3d
% Builds square matrix D2 representing discretization of Laplacian in 2D at
% given z value (in mm)

% Constants
global xdim
global ydim
global zdim
global h


numPoints = xdim*ydim*zdim;
sz = [xdim, ydim, zdim];

% Initialize D2 as sparse matrix
numNonzero = 7*numPoints;
x_lst = []
y_lst = []
val_lst = []
t_0 = datetime;
for ijk = 1:numPoints
    % Get (x,y,z) coordinates for point
    [i, j, k] = ind2sub(sz, ijk);

    % Update diagonal entry of D2
    x_lst(end + 1) = ijk;
    y_lst(end + 1) = ijk;
    val_lst(end + 1) = -6;

    % Indices of neighbors
    if i+1 <= xdim
        i1jk = sub2ind(sz, i+1, j, k);
    else
        i1jk = 0;
    end

    if i-1 >= 1
        im1jk = sub2ind(sz, i-1, j, k);
    else
        im1jk = 0;
    end

    if j+1 <= ydim
        ij1k = sub2ind(sz, i, j+1, k);
    else
        ij1k = 0;
    end

    if j-1 >= 1
        ijm1k = sub2ind(sz, i, j-1, k);
    else
        ijm1k = 0;
    end

     if k+1 <= zdim
         ijk1 = sub2ind(sz, i, j, k+1);
     else
         ijk1 = 0;
     end
     
     if k-1 >= 1
         ijkm1 = sub2ind(sz, i, j, k-1);
     else
         ijkm1 = 0;
     end

    neighbors = [i1jk im1jk ij1k ijm1k ijk1 ijkm1];

    % Check if orthogonal neighbors are in the boundary and update D2
    for index = 1:length(neighbors)
        neighbor = neighbors(index);
        if neighbor ~= 0 % Neighbor is in grid
            [nbr_x, nbr_y, nbr_z] = ind2sub(sz,neighbor);
            point = [nbr_x, nbr_y, nbr_z];
        else point = [1, 1, 1]; % just assign a random point
        end
        isGhost = (isBoundary(point) | neighbor == 0);
        isNotGhost = (isGhost == 0);

        % Adjust index for ghost nodes (either no change or set = ijk)
        neighbor = isNotGhost*neighbor + isGhost*ijk; 

        % Update row ijk of D2
        x_lst(end + 1) = ijk;
        y_lst(end + 1) = neighbor;
        val_lst(end + 1) = 1;
    end
end
D2 = sparse(x_lst, y_lst, val_lst, numPoints,numPoints);


end