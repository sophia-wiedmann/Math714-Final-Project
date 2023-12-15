%%% Reads in MRI data from simulated normal brain (BrainWeb)

% Starting measurements in mm
% zStart = -72;
% yStart = -126;
% xStart = -90;

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

% Read in the data
for i = 1:length(filenames)
    % Read to an array in volumeData
    filename = path + filenames(i);
    fileID = fopen(filename);
    volumeData{i} = fread(fileID);
    fclose(fileID);

    % Reshape to x y z dimensions
    volumeData{i} = reshape(volumeData{i},181,217,181);

    % Rescale MRI data to [0,1]
    volumeData{i} = 1/255*volumeData{i};
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
