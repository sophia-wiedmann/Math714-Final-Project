function plotMRI(z,volumeData)
%%% Plots 2D slice of given MRI data 
% z = height of slice in mm 
    % Should get same image as in BrainWeb with given z value
% volumeData is cell array of matrices created by readData


% Volume data containers
greyVol = volumeData{1};
whiteVol = volumeData{2};
skullVol = volumeData{3};
skinVol = volumeData{4};
connectiveVol = volumeData{5};
CSFVol = volumeData{6};
fatVol = volumeData{7};
glialVol = volumeData{8};
muscleVol = volumeData{9};
backgroundVol = volumeData{10};


%% Get image data for a 2D slice

% Starting measurements in mm
zStart = -72;
yStart = -126;
xStart = -90;

% Get matching index for z
zval = z-zStart;

% Get MRI image data for 2D slice at z = zval
greyData = greyVol(:,:,zval);
whiteData = whiteVol(:,:,zval);
skullData = skullVol(:,:,zval);
skinData = skinVol(:,:,zval);
connectiveData = connectiveVol(:,:,zval);
CSFData = CSFVol(:,:,zval);
fatData = fatVol(:,:,zval);
glialData = glialVol(:,:,zval);
muscleData = muscleVol(:,:,zval);
backgroundData = backgroundVol(:,:,zval);


%% Plot 2D slices

% Flip x and y for plotting with pcolor (color map)
greyData = greyData';
whiteData = whiteData';
skullData = skullData';
skinData = skinData';
connectiveData = connectiveData';
CSFData = CSFData';
fatData = fatData';
glialData = glialData';
muscleData = muscleData';
backgroundData = backgroundData';

% Create figures
figure;
s = pcolor(greyData);
s.FaceColor = 'interp';
colorbar;
title("Grey Matter");
axis image

figure;
s = pcolor(whiteData);
s.FaceColor = 'interp';
colorbar;
title("White Matter");
axis image

figure;
s = pcolor(skullData);
s.FaceColor = 'interp';
colorbar;
title("Skull");
axis image

figure;
s = pcolor(skinData);
s.FaceColor = 'interp';
colorbar;
title("Skin");
axis image

figure;
s = pcolor(connectiveData);
s.FaceColor = 'interp';
colorbar;
title("Connective Tissue"); % What is this?
axis image

figure;
s = pcolor(CSFData);
s.FaceColor = 'interp';
colorbar;
title("Cerebrospinal Fluid (CSF)");
axis image

figure;
s = pcolor(fatData);
s.FaceColor = 'interp';
colorbar;
title("Fat");
axis image

figure;
s = pcolor(glialData);
s.FaceColor = 'interp';
colorbar;
title("Glial Matter");
axis image

figure;
s = pcolor(muscleData);
s.FaceColor = 'interp';
colorbar;
title("Muscle and Skin");
axis image

figure;
s = pcolor(backgroundData);
s.FaceColor = 'interp';
colorbar;
title("Background");
axis image

end