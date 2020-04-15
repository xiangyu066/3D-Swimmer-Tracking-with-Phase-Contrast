echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: SwimmerSpatialTracker.m
% - Author: XYZ
% - Created date: April 9, 2020
% - Modified date: April 16, 202
% - Notes:
%       1.) This version have to manually choose target, and the target
%       have to roughly put in the middle.
% - Next modified:
%       1.) Automatically chosen target
% - Version: 1.1.2
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
close all, clear all, warning('off')
disp('Running...')

%% Define units
global um px sec msec
um = 1;
px = 1;
msec = 1;
sec = 1E+3 *(msec);

%% Define parameters of imaging system
% define imaging system parameters 
dz = 0.25 *(um);
pixelsize = 6.5 *(um);
Obj_Mag = 40; % magnification of objective
dt = 30 *(msec);
nFrames = 100;
sFrame = 1;
eFrame = nFrames;

%
speed_upper = 100 *(um/sec);
inputfile = 'C:\Users\motorsgroup\Downloads\20200409\RP437\Motile test\RP437_ROI052_BF.tif';
outputdir = 'C:\Users\motorsgroup\Downloads\20200409\RP437';
outputname = 'RP437_ROI052_BF-2';

%% Preallocating 1vavriables and functions
% load depth library
load('C:\Users\motorsgroup\Downloads\20200409\RP437\Lib\Library.mat')
load('C:\Users\motorsgroup\Downloads\20200409\RP437\Lib\Library_Pos.mat')

% constrain library depth range
Library(:,:,Library_Pos<-40*(um)) = '';
Library_Pos(Library_Pos<-40*(um)) = '';
Library(:,:,Library_Pos>50*(um)) = '';
Library_Pos(Library_Pos>50*(um)) = '';

Lib_Square = size(Library,1);
nLayers = size(Library,3);

% normalized library
Library = Library - mean(mean(Library,2),1); % remove mean value to obtain variation
Library = Library./std(Library,0,[1, 2]); % normalized

% Preallocating vavriables
eff_pixelsize = pixelsize/Obj_Mag;
search_Box_xy = Lib_Square + 2*round(speed_upper*dt/eff_pixelsize);
search_Box_z = 2*round(speed_upper*dt/dz);

% padding library
Lib_padding = zeros([search_Box_xy,search_Box_xy,nLayers]);
Lib_padding(1:Lib_Square,1:Lib_Square,:) = Library;

%
maxCrr_nLayer = 0;

% Make database folder
[~,msg] = mkdir('Analysis');

%%
vSession = VideoWriter([outputdir,'\Analysis\',outputname,'.avi']);
vSession.FrameRate = 33;
open(vSession)

figure(1), set(gcf,'WindowStyle','docked'), tic
for nFrame = sFrame:eFrame
    
    % read image
    origina = double(imread(inputfile,nFrame));
    
    % label tracking seeds
    if (nFrame ==sFrame)
        % human-made define target
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
        if ~isempty(maxCrr_nLayer)
            sx = sx+shiftX+Lib_Square/2-search_Box_xy/2;
            sy = sy+shiftY+Lib_Square/2-search_Box_xy/2;
            ex = sx+search_Box_xy-1;
            ey = sy+search_Box_xy-1;
        end
        sLayer = old_nLayer-search_Box_z/2;
        eLayer = old_nLayer+search_Box_z/2;
    end
    
    % determinate whether out of boundary
    checksum = (sx>=1)+(sy>=1)+(ex<=size(origina,2))+(ey<=size(origina,1));
    if checksum<4
        break;
    end
    
    % normalized
    neworigina = origina(sy:ey,sx:ex);
    neworigina = neworigina-mean2(neworigina);
    neworigina = neworigina./std2(neworigina);
    
    % using FFT to boost calculating correlation efficiency
    Crrs = zeros(nLayers,4);
    for nLayer = sLayer:eLayer
        Crr = abs( ifft2( fft2(neworigina) .* conj(fft2(Lib_padding(:,:,nLayer))) ) );
        [shiftY,shiftX] = find(Crr==max(Crr(:)));
        Crrs(nLayer,:) = [shiftX,shiftY,nLayer,max(Crr(:)) ] ;
    end
    
    % filt out of boundary
    if nFrame>1
        Crrs((Crrs(:,3)==0),:) = '';
        Crrs((Crrs(:,1)>size(neworigina,2)-Lib_Square+1),:) = '';
        Crrs((Crrs(:,2)>size(neworigina,1)-Lib_Square+1),:) = '';
    end
    
    % find maximal correlation coefficient
    maxCrr_nLayer = find(Crrs(:,4)==max(Crrs(:,4)));
    if isempty(maxCrr_nLayer)
        shiftX = 0;
        shiftY = 0;
    else
        shiftX = Crrs(maxCrr_nLayer,1);
        shiftY = Crrs(maxCrr_nLayer,2);
        old_shiftX = shiftX;
        old_shiftY = shiftY;
        old_nLayer = Crrs(maxCrr_nLayer,3);
    end

    % write databse
    if isempty(maxCrr_nLayer)
        Pos(nFrame,:) = [[sx+shiftX,sy+shiftY]*eff_pixelsize,Pos(nFrame-1,3)];
    else
        Pos(nFrame,:) = [[sx+shiftX,sy+shiftY]*eff_pixelsize,Library_Pos(Crrs(maxCrr_nLayer,3))];
    end
    
    % real-time monitoring
    figure(1), clf(gcf), 
    
    subplot 221, imshow(neworigina,[]), title(['nFrame = ', num2str(nFrame), ' (', num2str((nFrame-1)*dt), ' msec)'])
    hold on, rectangle('position',[old_shiftX,old_shiftY,Lib_Square,Lib_Square],'edgecolor','r')
    
    subplot 222, plot3(Pos(1,1),Pos(1,2),Pos(1,3), '>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on, plot3(Pos(:,1),Pos(:,2),Pos(:,3),'k-')
    hold on, plot3(Pos(end,1),Pos(end,2),Pos(end,3), 's','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    xlabel('X [\mum]','fontweight','bold'), ylabel('Y [\mum]','fontweight','bold'), zlabel('Z [\mum]','fontweight','bold')
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
close(vSession), toc

%% plot swimmer's trajectory
figure(3), set(gcf,'WindowStyle','docked')
subplot 221, plot3(Pos(1,1),Pos(1,2),Pos(1,3), '>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on, plot3(Pos(:,1),Pos(:,2),Pos(:,3),'k')
hold on, plot3(Pos(end,1),Pos(end,2),Pos(end,3), 's','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
xlabel('X [\mum]','fontweight','bold'), ylabel('Y [\mum]','fontweight','bold'), zlabel('Z [\mum]','fontweight','bold')
grid on, axis equal

subplot 222, plot(Pos(1,1),Pos(1,2), '>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on, plot(Pos(end,1),Pos(end,2), 's','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on, plot(Pos(:,1),Pos(:,2),'k')
xlabel('X [\mum]','fontweight','bold'), ylabel('Y [\mum]','fontweight','bold')
legend({'Start','End'})
grid on, axis equal

subplot 223, plot(Pos(1,1),Pos(1,3), '>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on, plot(Pos(end,1),Pos(end,3), 's','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on, plot(Pos(:,1),Pos(:,3),'k')
xlabel('X [\mum]','fontweight','bold'), ylabel('Z [\mum]','fontweight','bold')
legend({'Start','End'})
grid on, axis equal

subplot 224, plot(Pos(1,2),Pos(1,3), '>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on, plot(Pos(end,2),Pos(end,3), 's','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on, plot(Pos(:,2),Pos(:,3),'k')
xlabel('Y [\mum]','fontweight','bold'), ylabel('Z [\mum]','fontweight','bold')
legend({'Start','End'})
grid on, axis equal

frame = getframe(gcf);
imwrite(frame.cdata, [outputdir,'\Analysis\',outputname,'.png'])

% save tracking data
data.unit = {'um', 'msec'};
data.dt = dt;
data.Pos = Pos;
save([outputdir,'\Analysis\',outputname,'.mat'], 'data')

%%
disp('Done.')
