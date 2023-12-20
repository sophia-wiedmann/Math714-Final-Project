function D1 = buildGradient_3d
% Builds square matrix D1 representing discretization of gradient in 2D at
% given z value (in mm)
% Uses centered difference approximation for spatial derivative

% Constants
global xdim
global ydim
global zdim

numPoints = xdim*ydim*zdim;
sz = [xdim, ydim, zdim];

% Initialize D2 as sparse matrix
numNonzero = 7*numPoints;
t_0 = datetime
x_lst = []
y_lst = []
val_lst = []
for ijk = 1:numPoints
    % Get (x,y,z) coordinates for point
    [i, j, k] = ind2sub(sz, ijk);

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


    % Check if orthogonal neighbors are in the boundary and update D1
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
        neighbors(index) = isNotGhost*neighbor + isGhost*ijk; 
    end

    % Update row ijk of D1
    for i =1:length(neighbors)
        x_lst(end + 1) = ijk;
        y_lst(end + 1) = neighbors(i);
        val_lst(end+1) = (-1)^(i+1);
    end

end
D1 = sparse(x_lst, y_lst, val_lst, numPoints,numPoints);


end