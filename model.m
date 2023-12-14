A = [1 2 3; 4 5 6; 7 8 9];
D_test = cat(3, [1 2 3; 4 5 6; 7 8 9], [10 11 12; 13 14 15; 16 17 18], [19 20 21; 22 23 24; 25 26 27]);
disp(D_test);
[m, n, p] = size(D_test);
global k;
global h;
global n;
global mu;
k = 1/10
h = 1/181
n = 181
mu = k/(2*h^2)


for i = 1:m
    for j = 1:n
        for k = 1:p
            % Access the element at (i, j, k)
            elementValue = Dx(D_test, i, j, k, m);
            
            % Perform operations on the element if needed
            % (e.g., display the value)
            fprintf('Position (%d, %d, %d): %f\n', i, j, k, elementValue);
        end
    end
end



function diff = Dx(D, i, j, k, n)
    if i == 1 || i == n
        diff = 0;
    else 
        diff = D(j, i+1, k) - D(j, i-1, k);
    end
end
function diff = Dy(D, i, j, k, n)
    if n == 1 || j == n
        diff = 0;
    else 
        diff = D(j+1, i, k) - D(j-1, i, k);
    end
end

function diff = Dz(D, i, j, k, n)
    if k == 1 || k == n
        diff = 0;
    else 
        diff = D(j, i, k+1) - D(j, i, k-1);
    end
end

function ind = index(i, j, k, n)
    ind = i + n*(j-1) + n^2*(k-1)
end


function c = diffusion_x(D, C, n, timestep)
    D_x = zeros(1,n^3)
    for k = 1:n
        for j = 1:n
            for i = 1:n
                % Access the element at (i, j, k)
                D_x(int64(index(i, j,k, n))) = Dx(D_test, i, j, k, n);
                
                % Perform operations on the element if needed
                % (e.g., display the value)
                fprintf('Position (%d, %d, %d): %f\n', i, j, k, elementValue);
            end
        end
    end

    A = diag(((1/mu) - 2*reshape(D, 1, [])) *ones(1, n^3)) + diag(((-D_x/4)-reshape(D, 1, []))*ones(1, n^3), 1) + diag(((+D_x/4)-reshape(D, 1, []))*ones(1, n^3), -1);
    B = diag(((1/mu) - 2*reshape(D, 1, [])) *ones(1, n^3))+ diag(((D_x/4)+reshape(D, 1, []))*ones(1, n^3), 1) + diag(((-D_x/4)+reshape(D, 1, []))*ones(1, n^3), -1);
    
    c = A \ (B * C)
end

function c = diffusion_y(D, C, n, timestep)
    D_y = zeros(1,n^3)
    D_vec = reshape(D, 1, [])
    for k = 1:n
        for j = 1:n
            for i = 1:n
                % Add element to vector at (i, j, k)
                D_y(int64(index(i, j,k, n))) = Dy(D_test, i, j, k, n);
                
            end
        end
    end
    A_blk_lst = []
    for i in 1:n
        A_blk = diag(((1/mu) - 2*D_vec) *ones(1, n^2)) + diag(((-D_y((i-1)*n+1:i*n)/4)-D_vec((i-1)*n+1:i*n)/4))*ones(1, n^2 - n), n+1) + diag(((+D_y((i-1)*n+1:i*n)/4)-D_vec((i-1)*n+1:i*n)/4))*ones(1, n^2 - n), -1-n)
        B_blk = diag(((1/mu) - 2*D_vec) *ones(1, n^2)) + diag(((D_y((i-1)*n+1:i*n)/4)+D_vec((i-1)*n+1:i*n)/4))*ones(1, n^2 - n), n+1) + diag(((-D_y((i-1)*n+1:i*n)/4)+D_vec((i-1)*n+1:i*n)/4))*ones(1, n^2 - n), -1-n)

        if i == 1
            A = A_blk
            B = B_blk
        else
            A = blkdiag(A, A_blk)
            B = blkdiag(B, B_blk)

    c = A \ (B * C)
end


function c = diffusion_z(D, C, n, timestep)
    D_z = zeros(1,n^3)
    for k = 1:n
        for j = 1:n
            for i = 1:n
                % Add element to vector at (i, j, k)
                D_z(int64(index(i, j,k, n))) = Dz(D_test, i, j, k, n);
                
            end
        end
    end
    A_blk = []
    A = diag(((1/mu) - 2*reshape(D, 1, [])) *ones(1, n^3)) + diag(((-D_z/4)-reshape(D, 1, []))*ones(1, n^3), n^2+1) + diag(((+D_x/4)-reshape(D, 1, []))*ones(1, n^3), -1-n);
    B = diag(((1/mu) - 2*reshape(D, 1, [])) *ones(1, n^3))+ diag(((D_z/4)+reshape(D, 1, []))*ones(1, n^3), n^2+1) + diag(((-D_x/4)+reshape(D, 1, []))*ones(1, n^3), -1-n);
    
    c = A \ (B * C)
end


