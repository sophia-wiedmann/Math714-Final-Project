
%% Data and model parameters
% Dimensions of MRI  data
global xdim
global ydim
global zdim

% MRI data
global skullVol
global backgroundVol
global greyVol
global whiteVol

% Constants
global h
%global b
global k
global rho
global chem_kill_rate
global rad_kill_rate

plot = true;
treat = true;

xdim = 181;
ydim = 217;
zdim = 181;

% Tissue-specific diffusion coefficients from Table 11.6 in textbook
    % Units are cm^2/day
    % Tumor grading is high (HH), intermediate (HL), intermediate (LH),
    % and low (LL)
Dg = 1.3*10^(-3); % HH
%Dg = 1.3*10^(-4); % HL
%Dg = 1.3*10^(-3); % LH
%Dg = 1.3*10^(-4); % LL
Dw = 5*Dg; % max diffusion coefficient

rho = 1.2*10^(-2); % HH
%rho = 1.2*10^(-2); % HL
%rho = 1.2*10^(-3); % LH
%rho = 1.2*10^(-3); % LL

chem_kill_rate = .06;
rad_kill_rate = .9;

detection_threshold = 0.0019;

h = 0.1; % h = 1mm = 0.1cm
k = 1/ceil(1/(h^2/(6*Dw))); % choose k <= h^2/(6*Dw)
az = -37.5;
% Read in data
readData

% Constants

numPoints = xdim*ydim*zdim;

% Matrix for spatial discretization
load('Matricies/F.mat');

%% Initial Condition (normal distribution -- see 11.9 in book)
x0 = [111, 50, 111]; % center of tumor
a = 1; % max density at center of tumor
r = 8.3; % radius of tumor in mm
cutoff = 0.01; % density at radius r
b = -r^2/log(cutoff/a); % measure of spread so that cutoff condition is satisfied

% Initialize IC
IC = zeros(xdim,ydim,zdim);

% Compute IC at each grid point
for x = 1:xdim
    for y = 1:xdim
        for z = 1:zdim
            dist2 = (x-x0(1))^2 + (y-x0(2))^2 + (z - x0(3))^2; % squared distance to center of tumor
            IC(x,y, z) = a*exp(-dist2/b);
        end
    end
end


%% Plot IC
        % Grab just the data where the tumor concentration is > threshold
        X = [];
        Y = [];
        Z = [];
        val = [];
        
        % Reshape to matrix
        C_n = IC;
        for x = 1:xdim
            for y = 1:ydim
                for z = 1:zdim
                    if C_n(x, y, z)>= detection_threshold
                        X(end + 1) = x;
                        Y(end + 1) = y;
                        Z(end + 1) = z;
                        val(end+1) = C_n(x, y, z);
                    end
                end
            end
        end

        r = findRadius(X, Y, Z);
        tot_concent = sum(C_n, "all");
        if treat
            rad_arr_treat = [double(0.0), r, tot_concent]
        else
            rad_arr = [double(0.0), r, tot_concent];
        end
        
        
        
        
        az = -37.5;
        w = (C_n + .01)*100;
        if plot
            if ~treat
                % Relevant data for plotting
                scatter3(X, Y, Z, val*100 ,'filled');
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/tumorstill0000.png');
                saveas(gcf, fname);
                xlim([95, 125]);
                ylim([35, 65]);
                zlim([100, 130]);
                fname = append('tumor_images/tumorstillzoom0000.png');
                saveas(gcf, fname);
                view(az, 30);
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/tumor0000.png');
                saveas(gcf, fname);
                xlim([95, 125]);
                ylim([35, 65]);
                zlim([100, 130]);;
                fname2 = append('tumor_images/tumorzoom0000.png');
                saveas(gcf, fname2);
                az = az + 1;
            else
                % Relevant data for plotting
                scatter3(X, Y, Z, val*100 ,'filled');
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/treattumorstill0000.png');
                saveas(gcf, fname);
                xlim([95, 125]);
                ylim([35, 65]);
                zlim([100, 130]);
                fname = append('tumor_images/treattumorstillzoom0000.png');
                saveas(gcf, fname);
                view(az, 30);
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/treattumor0000.png');
                saveas(gcf, fname);
                xlim([95, 125]);
                ylim([35, 65]);
                zlim([100, 130]);
                fname2 = append('tumor_images/treattumorzoom0000.png');
                saveas(gcf, fname2);
                az = az + 1;

            end
        end

%% Simulate tumor growth
C_n = reshape(IC,numPoints,1);
t_0 = datetime
% Simulate for 365 time steps
for t = 1:365
    if mod(t, 10) == 0
    t
    datetime - t_0
    end
    C = C_n;
    C_n = solver_3d(C,F, t, treat);
        %% Plot results
        % Grab just the data where the tumor concentration is > 10^-6
        X = [];
        Y = [];
        Z = [];
        val = [];
        
        % Reshape to matrix
        C_n = reshape(C_n,xdim,ydim, zdim);
        for x = 1:xdim
            for y = 1:ydim
                for z = 1:zdim
                    if C_n(x, y, z)> detection_threshold
                        X(end + 1) = x;
                        Y(end + 1) = y;
                        Z(end + 1) = z;
                        val(end+1) = C_n(x, y, z);
                    end
                end
            end
        end
        
        
        r = findRadius(X, Y, Z);
        tot_concent = sum(C_n, "all");
        if treat
            rad_arr_treat = [rad_arr_treat; t, r, tot_concent];
        else
            rad_arr = [rad_arr; t, r, tot_concent];
        end
        
        
        w = (C_n + .01)*100;
        if plot
            if ~treat
                % Relevant data for plotting
                
                scatter3(X, Y, Z, val*100 ,'filled');
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/tumorstill', num2str(t, '%04d'), '.png');
                saveas(gcf, fname); 
                xlim([95, 125]);
                ylim([35, 65]);
                zlim([100, 130]);
                fname = append('tumor_images/tumorstillzoom', num2str(t, '%04d'), '.png');
                saveas(gcf, fname);
                view(az, 30);
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/tumor', num2str(t, '%04d'), '.png');
                saveas(gcf, fname);
                xlim([95, 125]);
                ylim([35, 65]);
                zlim([100, 130]);
                fname2 = append('tumor_images/tumorzoom', num2str(t, '%04d'), '.png');
                saveas(gcf, fname2);
                az = az + 1;
            else
                % Relevant data for plotting
                scatter3(X, Y, Z, val*100 ,'filled');
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/treattumorstill', num2str(t, '%04d'), '.png');
                saveas(gcf, fname);
                xlim([95, 125]);
                ylim([35, 65]);
                zlim([100, 130]);
                fname = append('tumor_images/treattumorstillzoom', num2str(t, '%04d'), '.png');
                saveas(gcf, fname);
                view(az, 30);
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/treattumor', num2str(t, '%04d'), '.png');
                saveas(gcf, fname);
                xlim([95, 125]);
                ylim([35, 65]);
                zlim([100, 130]);
                fname2 = append('tumor_images/treattumorzoom', num2str(t, '%04d'), '.png');
                saveas(gcf, fname2);
                az = az + 1;


            end

        end
         C_n = reshape(C_n,numPoints,1);
         max(C_n);
    
end
total_time = datetime - t_0
if treat
    rad_arr_treat
    clear plot
    xValues = rad_arr_treat(:, 1);
    yValues = rad_arr_treat(:, 2);
    plot(xValues, yValues, '-x');
    ylim([0 inf])
    xlabel('Time (days)', fontsize=16) 
    ylabel('Tumor Radius (mm)', fontsize = 16) 
    title('Plot of Radius of Tumor with Chemotherapy and Radiation Treatment over Time', fontsize = 20)

else
    rad_arr
    clear plot
    xValues1 = rad_arr(:, 1);
    yValues1 = rad_arr(:, 2);
    plot(xValues1, yValues1, '-x');
    ylim([0 inf])
    xlabel('Time (days)', fontsize=16) 
    ylabel('Tumor Radius (mm)', fontsize = 16) 
    title('Plot of Radius of Tumor over Time', fontsize = 20)

end




