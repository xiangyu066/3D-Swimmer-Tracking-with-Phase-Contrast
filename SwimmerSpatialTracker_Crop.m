echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: SwimmerSpatialTracker_Crop.m
% - Author: XYZ
% - Created date: April 11, 2020
% - Modified date: April 15, 2020
% - Notes:
%       1.) 
% - Version: 1.0.1
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
close all, clear all, warning('off')
disp('Running...')

%%
global px
px = 1;

%% Define image infomation
inputfile = [pwd, '\Lib\', 'Ecoli_RP437_40X_Ph2_collar1.2_FM4-64cube_ROI10_5.tif'];
nLayers = 601;
width = 512 *(px);
height = 512 *(px);
roi_Square = 400 *(px);

outputfile = 'ROI10_5.tif'

%% Preallocating vavriables and functions
Proj_Images = zeros(height, width);

%% Do maximal projection for all chosen cells
tic
for nLayer = 1:nLayers
    origina = double(imread(inputfile, nLayer));
    Proj_Images(origina>Proj_Images) = origina(origina>Proj_Images);
end
toc

%% Crop
figure, imshow(Proj_Images, [])
h = imrect(gca,[1, 1, roi_Square, roi_Square]);
position = wait(h);
sx = round(position(1));
sy = round(position(2));
ex = sx+roi_Square-1;
ey = sy+roi_Square-1;

% save stacked image
neworigina = zeros(width,height);
for nLayer = 1:nLayers
    origina = double(imread(inputfile, nLayer));
    neworigina = origina(sy:ey,sx:ex);
    
    if (nLayer <2)
        imwrite(uint16(neworigina),outputfile,'Compression','none')
    else
        imwrite(uint16(neworigina),outputfile,'WriteMode','append','Compression','none')
    end
end

%%
disp('Done')
