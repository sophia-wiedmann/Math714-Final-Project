function D1 = buildGradient_3d
% Builds square matrix D1 representing discretization of gradient in 2D at
% given z value (in mm)
% Uses centered difference approximation for spatial derivative

% Constants
global xdim
global ydim
global zdim

%zStart = -72;
%zval = z-zStart;

numPoints = xdim*ydim*zdim;
sz = [xdim, ydim, zdim];
%numPoints = xdim*ydim;
%sz = [xdim, ydim];

% Initialize D2 as sparse matrix
numNonzero = 7*numPoints;
%numNonzero = 4*numPoints;
%D1 = spalloc(numPoints,numPoints,numNonzero);
t_0 = datetime
x_lst = []
y_lst = []
val_lst = []
for ijk = 1:numPoints
    if ~mod(ijk, 100000)
        ijk
        datetime - t_0
        t_1 = datetime;
    end
    % Get (x,y,z) coordinates for poibuildGradient(z)nt
    [i, j, k] = ind2sub(sz, ijk);
    %[i, j] = ind2sub(sz, ijk);

    % Indices of neighbors
    if i+1 <= xdim
        i1jk = sub2ind(sz, i+1, j, k);
        %i1jk = sub2ind(sz, i+1, j);
    else
        i1jk = 0;
    end

    if i-1 >= 1
        im1jk = sub2ind(sz, i-1, j, k);
        %im1jk = sub2ind(sz, i-1, j);
    else
        im1jk = 0;
    end

    if j+1 <= ydim
        ij1k = sub2ind(sz, i, j+1, k);
        %ij1k = sub2ind(sz, i, j+1);
    else
        ij1k = 0;
    end

    if j-1 >= 1
        ijm1k = sub2ind(sz, i, j-1, k);
        %ijm1k = sub2ind(sz, i, j-1);
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
    %neighbors = [i1jk im1jk ij1k ijm1k];
    if ~mod(ijk, 1000)
        %milliseconds(datetime - t_1)
        t_1 = datetime;

    end

    % Check if orthogonal neighbors are in the boundary and update D1
    for index = 1:length(neighbors)
        neighbor = neighbors(index);
        if neighbor ~= 0 % Neighbor is in grid
            [nbr_x, nbr_y, nbr_z] = ind2sub(sz,neighbor);
            point = [nbr_x, nbr_y, nbr_z];
            %[nbr_x, nbr_y] = ind2sub(sz,neighbor);
            %point = [nbr_x, nbr_y, zval];
        else point = [1, 1, 1]; % just assign a random point
        end
        isGhost = (isBoundary(point) | neighbor == 0);
        isNotGhost = (isGhost == 0);

        % Adjust index for ghost nodes (either no change or set = ijk)
        neighbors(index) = isNotGhost*neighbor + isGhost*ijk; 
    end
    if ~mod(ijk, 1000)
        %milliseconds(datetime - t_1)
        t_1 = datetime;

    end


    % Update row ijk of D1
    for i =1:length(neighbors)
        x_lst(end + 1) = ijk;
        y_lst(end + 1) = neighbors(i);
        val_lst(end+1) = (-1)^(i+1);
    end
    %D1(ijk,neighbors(1)) = 1;
    %D1(ijk,neighbors(2)) = -1;
    %D1(ijk,neighbors(3)) = 1;
    %D1(ijk,neighbors(4)) = -1;
    %D1(ijk,neighbors(5)) = 1;
    %D1(ijk,neighbors(6)) = -1;


end
D1 = sparse(x_lst, y_lst, val_lst, numPoints,numPoints);


end