echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: SwimmerSpatialTracker.m
% - Author: XYZ
% - Created date: April 9, 2020
% - Modified date: April 13, 202
% - Notes:
%       1.) 
% - Next modified:
%       1.) Automatically chosen target
% - Version: 1.0.0
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
close all, clear all, warning('off')
disp('Running...')

%% Define units
global um px sec msec
um = 1;
px = 1;
sec = 1;
msec = 1E-3 *(sec);

%% Define parameters of imaging system
% define imaging system parameters 
dz = 0.25 *(um);
pixelsize = 6.5 *(um);
Mag = 40; % magnification of objective
dt = 30 *(msec);

%
speed_upper = 100 *(um/sec);
inputfile = 'D:\20200409\Motile test\RP437_ROI038_BF.tif';
outputdir = 'D:\20200409\Analysis\';
outputname = 'RP437_ROI038_BF-3';

start_frame = 1;
end_frame = 100;

%% Preallocating 1vavriables and functions
% Load depth library
load('Library.mat')
load('Library_Pos.mat')

Lib_Square = size(Library,1);
nLayers = size(Library,3);
Library = Library - mean(mean(Library,2),1); % remove mean value to obtain variation
Library = Library./std(Library,0,[1, 2]); % normalized

% Preallocating vavriables
eff_pixelsize = pixelsize/Mag;
search_Box_xy = Lib_Square + 2*round(speed_upper*dt/eff_pixelsize);
search_Box_z = 2*round(speed_upper*dt/dz);

% padding library
Lib_padding = zeros([search_Box_xy,search_Box_xy,nLayers]);
Lib_padding(1:Lib_Square,1:Lib_Square,:) = Library;

%
idx = [];

% Make database folder
% [~,msg] = mkdir('Analysis');

%%
vSession = VideoWriter([outputdir,outputname,'.avi']);
vSession.FrameRate = 33;
open(vSession)

figure(1), set(gcf,'WindowStyle','docked')
for nFrame = start_frame:end_frame
    
    % read image
    origina = double(imread(inputfile,nFrame));
    
    % human-made define target
    if nFrame ==start_frame
        figure(2), imshow(origina,[])
        h = imrect(gca,[1, 1, search_Box_xy, search_Box_xy]);
        position = wait(h);
        sx = round(position(1));
        sy = round(position(2));
        ex = sx+search_Box_xy-1;
        ey = sy+search_Box_xy-1;
        sLayer = 1;
        eLayer = nLayers;
    else
        if ~isempty(idx)
            sx = sx+shiftX+Lib_Square/2-search_Box_xy/2;
            sy = sy+shiftY+Lib_Square/2-search_Box_xy/2;
            ex = sx+search_Box_xy-1;
            ey = sy+search_Box_xy-1;
        end
        sLayer = old_nLayer-search_Box_z/2;
        eLayer = old_nLayer+search_Box_z/2;
    end
    
    neworigina = origina(sy:ey,sx:ex);
    neworigina = neworigina-mean2(neworigina);
    neworigina = neworigina./std2(neworigina);
    
    % using FFT to reduce calculating correlation time
    tic
    Crrs = zeros(nLayers,4);
    for nLayer = sLayer:eLayer
        Crr = abs( ifft2( fft2(neworigina) .* conj(fft2(Lib_padding(:,:,nLayer))) ) );
        [shiftY,shiftX] = find(Crr==max(Crr(:)));
        Crrs(nLayer,:) = [shiftX,shiftY,nLayer,max(Crr(:)) ] ;
    end
    toc
    
    % filt out of boundary 
    if nFrame>1
        Crrs((Crrs(:,3)==0),:) = '';
        Crrs((Crrs(:,1)>size(neworigina,2)-Lib_Square+1),:) = '';
        Crrs((Crrs(:,2)>size(neworigina,1)-Lib_Square+1),:) = '';
    end
    
    % find maximal correlation coefficient
    idx = find(Crrs(:,4)==max(Crrs(:,4)));
    if isempty(idx)
        shiftX = 0;
        shiftY = 0;
    else
        shiftX = Crrs(idx,1);
        shiftY = Crrs(idx,2);
        old_shiftX = shiftX;
        old_shiftY = shiftY;
        old_nLayer = Crrs(idx,3);
    end

    % write databse
    if isempty(idx)
        data(nFrame,:) = [[sx+shiftX,sy+shiftY]*eff_pixelsize,data(nFrame-1,3)];
    else
        data(nFrame,:) = [[sx+shiftX,sy+shiftY]*eff_pixelsize,Library_Pos(Crrs(idx,3))];
    end
    
    % real-time monitoring
    figure(1), clf(gcf), 
    
    subplot 221, imshow(neworigina,[]), title(['nFrame = ', num2str(nFrame)])
    hold on, rectangle('position',[old_shiftX,old_shiftY,Lib_Square,Lib_Square],'edgecolor','r')
    
    subplot 222, plot3(data(1,1),data(1,2),data(1,3), 's','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on, plot3(data(:,1),data(:,2),data(:,3),'k-')
    hold on, plot3(data(end,1),data(end,2),data(end,3), '>','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    xlabel('X [\mum]'), ylabel('Y [\mum]'), zlabel('Z [\mum]')
    grid on
    
    subplot 223, imshow(neworigina(old_shiftY:old_shiftY+Lib_Square-1,old_shiftX:old_shiftX+Lib_Square-1), []), title('target roi')
    hold on, line([Lib_Square/2+1,Lib_Square/2+1],[1,Lib_Square],'color','r')
    line([1,Lib_Square],[Lib_Square/2+1,Lib_Square/2+1],'color','r')
    
    subplot 224, imshow(Library(:,:,old_nLayer), []), title(['nLayer = ', num2str(old_nLayer), ' (',num2str(Library_Pos(old_nLayer)),' \mum)'])
    hold on, line([Lib_Square/2+1,Lib_Square/2+1],[1,Lib_Square],'color','r')
    line([1,Lib_Square],[Lib_Square/2+1,Lib_Square/2+1],'color','r')
    
    drawnow
    
    % save tracking process
    frame = getframe(gcf);
    writeVideo(vSession, frame);
end
close(vSession)
save([outputdir,outputname,'.mat'],'data')

%%
%
figure, set(gcf,'WindowStyle','docked')
subplot 221, plot3(data(2:end,1),data(2:end,2),data(2:end,3))
xlabel('X [\mum]'), ylabel('Y [\mum]'), zlabel('Z [\mum]')
grid on, axis equal
subplot 222, plot(data(2:end,1),data(2:end,2))
xlabel('X [\mum]'), ylabel('Y [\mum]')
grid on, axis equal
subplot 223, plot(data(2:end,1),data(2:end,3))
xlabel('X [\mum]'), ylabel('Z [\mum]')
grid on, axis equal
subplot 224, plot(data(2:end,2),data(2:end,3))
xlabel('Y [\mum]'), ylabel('Z [\mum]')
grid on, axis equal

frame = getframe(gcf);
imwrite(frame.cdata, [outputdir,outputname,'.png'])

%%
disp('Done.')
