%%% Script to read in MRI data and plot images

%% Read in anatomical structure data

% File names
path = "./Brain_Structure/";
grey = "phantom_1.0mm_normal_gry.rawb";
white = "phantom_1.0mm_normal_wht.rawb";
skull = "phantom_1.0mm_normal_skl.rawb";
skin = "phantom_1.0mm_normal_skn.rawb";
connective = "phantom_1.0mm_normal_mit.rawb";
CSF = "phantom_1.0mm_normal_csf.rawb";
fat = "phantom_1.0mm_normal_fat.rawb";
glial = "phantom_1.0mm_normal_gli.rawb";
muscle = "phantom_1.0mm_normal_m+s.rawb";
background = "phantom_1.0mm_normal_bck.rawb";

filenames = [grey white skull skin connective CSF ...
    fat glial muscle background];

% Container for all the volume data
volumeData = {};

% Read out the data
for i = 1:length(filenames)
    % Read to an array in volumeData
    filename = path + filenames(i);
    fileID = fopen(filename);
    volumeData{i} = fread(fileID);
    fclose(fileID);

    % Reshape to x y z dimensions
    volumeData{i} = reshape(volumeData{i},181,217,181);
end

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

% Slice measurement in mm (change this to get a different slice)
% Should get same image as in BrainWeb with given z value
z = 28;

% Get matching index
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

% Flip x and y for plotting with surf
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

figure;
surf(greyData);
colorbar;
title("Grey Matter");
axis equal

figure;
surf(whiteData);
colorbar;
title("White Matter");
axis equal

figure;
surf(skullData);
colorbar;
title("Skull");
axis equal

figure;
surf(skinData);
colorbar;
title("Skin");
axis equal

figure;
surf(connectiveData);
colorbar;
title("Connective Tissue"); % What is this?
axis equal

figure;
surf(CSFData);
colorbar;
title("Cerebrospinal Fluid (CSF)");
axis equal

figure;
surf(fatData);
colorbar;
title("Fat");
axis equal

figure;
surf(glialData);
colorbar;
title("Glial Matter");
axis equal

figure;
surf(muscleData);
colorbar;
title("Muscle and Skin");
axis equal

figure;
surf(backgroundData);
colorbar;
title("Background");
axis equal