global plot
global treat

%% Input Parameters
plot = false;
treat = true;

%% Tumor grade
% Grading of tumor (1 = HH, 2 = HL, 3 = LH, 4 = LL)
type = 1;

% Growth rate [1/day] and tissue-specific diffusion coefficients [mm^2/day]
% Tumor grading is high (HH), intermediate (HL), intermediate (LH), and low (LL)
if type == 1
    rho = 1.2*10^(-2); % HH
    Dg = 1.3*10^(-3)*10^2; % HH
elseif type == 2
    rho = 1.2*10^(-2); % HL
    Dg = 1.3*10^(-4)*10^2; % HL
elseif type == 3
    rho = 1.2*10^(-3); % LH
    Dg = 1.3*10^(-3)*10^2; % LH
elseif type == 4
    rho = 1.2*10^(-3); % LL
    Dg = 1.3*10^(-4)*10^2; % LL
else disp("Type must be an integer between 1 and 4")
end
Dw = 5*Dg; % max diffusion coefficient


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
global k
global rho
global chem_kill_rate
global rad_kill_rate

xdim = 181;
ydim = 217;
zdim = 181;

chem_kill_rate = .03;
rad_kill_rate = 1;

detection_threshold = 400;

h = 1; % h = 1mm = 0.1cm
k = 1/ceil(1/(h^2/(6*Dw))); % choose k <= h^2/(6*Dw)
az = -37.5;

% Read in data
readData

% Constants

numPoints = xdim*ydim*zdim;

% Matrix for spatial discretization
load('Matricies/F.mat');
kill_time = 600;
detect_time = 600;

%% Initial Condition (normal distribution -- see 11.9 in book)
x0 = [111, 51, 111]; % center of tumor
a = 1600; % max density at center of tumor
r = 15; % radius of tumor in mm
cutoff = 400; % density at radius r
b = -r^2/log(cutoff/a); % measure of spread so that cutoff condition is satisfied
%b = 300;
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


% Plot IC
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

        % Volume Data
        vol = length(X);
        if vol > 14137
            detect_time = 0;
        end

        if treat
            rad_arr_treat = [double(0.0), vol]
        else
            rad_arr = [double(0.0), vol];
        end
        
        
        
        %Plotting
        az = -37.5;
        if plot
            if ~treat
                % Relevant data for plotting
                scatter3(X, Y, Z, log(val + 1)/100 + 1 ,'filled');
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                xlabel('x (mm)')
                ylabel('y (mm)')
                zlabel('z (mm)')
                title('Concentration of Tumor','t = 0')
                fname = append('tumor_images/tumorstill0000.png');
                saveas(gcf, fname);
                view(az, 30);
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/tumor0000.png');
                saveas(gcf, fname);
                az = az + 1;
            else
                % Relevant data for plotting
                scatter3(X, Y, Z, log(val + 1)/100 + 1 ,'filled');
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                xlabel('x (mm)')
                ylabel('y (mm)')
                zlabel('z (mm)')
                title('Concentration of Tumor with Treatment','t = 0')
                fname = append('tumor_images/treattumorstill0000.png');
                saveas(gcf, fname);
                view(az, 30);
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/treattumor0000.png');
                saveas(gcf, fname);
                az = az + 1;

            end
        end

%% Simulate tumor growth
C_n = reshape(IC,numPoints,1);
t_0 = datetime
% Simulate for 365 time steps
for t = 1:365/k
    C = C_n;
    C_n = solver_3d(C,F, t*k, treat);
    if mod(k*t, 1) == 0 % Only plot or compute volume if time is an integer
        t*k
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
        
        % Volume Data
        vol = length(X);
        if vol > 113097
            kill_time = min(kill_time, t*k);
        end
        if vol > 14137
            detect_time = min(detect_time, t*k);
        end


        if treat
            rad_arr_treat = [rad_arr_treat; t*k, vol];
        else
            rad_arr = [rad_arr; t*k, vol];
        end
        
        % Create Plots
        if plot
            if ~treat
                % Relevant data for plotting
                
                scatter3(X, Y, Z, log(val + 1)/100 + 1 ,'filled');
                xlabel('x (mm)')
                ylabel('y (mm)')
                zlabel('z (mm)')
                title('Concentration of Tumor', append('t = ', int2str(t*k)))
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/tumorstill', num2str(t*k, '%04d'), '.png');
                saveas(gcf, fname); 
                view(az, 30);
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/tumor', num2str(t*k, '%04d'), '.png');
                saveas(gcf, fname);
                az = az + 1;
            else
                % Relevant data for plotting
                scatter3(X, Y, Z, log(val + 1)/100 + 1,'filled');
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                xlabel('x (mm)')
                ylabel('y (mm)')
                zlabel('z (mm)')
                title('Concentration of Tumor with Treatment', append('t = ', int2str(t*k)))
                fname = append('tumor_images/treattumorstill', num2str(t*k, '%04d'), '.png');
                saveas(gcf, fname);
                view(az, 30);
                xlim([0, 181]);
                ylim([0, 218]);
                zlim([0, 181]);
                fname = append('tumor_images/treattumor', num2str(t*k, '%04d'), '.png');
                saveas(gcf, fname);
                az = az + 1;


            end

        end
    end
         C_n = reshape(C_n,numPoints,1);    
end

%Analysis of vollume data
total_time = datetime - t_0
if treat
    rad_arr_treat
    clear plot
    xValues = rad_arr_treat(:, 1);
    yValues = rad_arr_treat(:, 2);
    plot(xValues, yValues, '-x');
    ylim([0 inf])
    xlabel('Time (days)', fontsize=16) 
    ylabel('Tumor Volume (mm^3)', fontsize = 16) 
    title('Plot of Volume of Tumor with Chemotherapy and Radiation Treatment over Time', fontsize = 20)

else
    rad_arr
    clear plot
    xValues1 = rad_arr(:, 1);
    yValues1 = rad_arr(:, 2);
    plot(xValues1, yValues1, '-x');
    ylim([0 inf])
    xlabel('Time (days)', fontsize=16) 
    ylabel('Tumor Volume (mm^3)', fontsize = 16) 
    title('Plot of Volume of Tumor over Time', fontsize = 20)

end
kill_time - detect_time
detect_time




