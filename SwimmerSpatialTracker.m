echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: SwimmerSpatialTracker.m
% - Author: XYZ
% - Created date: April 9, 2020
% - Modified date: May 7, 2020
% - Notes:
%       1.) This version have to manually choose target, and the target
%       have to roughly put in the middle.
% - Next modified:
%       1.) Automatically chosen target
% - Version: 2.1.0
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
clear all, warning('off')
disp('Running...')

%% Define units
global um px sec msec
um = 1;
px = 1;
msec = 1;
sec = 1E+3 *(msec);

%% Define parameters of imaging system
% define imaging system parameters 
dz = 0.25 *(um);                    % the axial depth between layers
pixelsize = 6.5 *(um);
Obj_Mag = 40;                       % the magnification of objective
dt = 30 *(msec);                    % the timestep between images
nFrames = 100;                      % the number of stacked images
sFrame = 1;                         % the starting tracking frame
eFrame = nFrames;                   % the ending tracking frame

%
speed_upper = 100 *(um/sec);
sz = -40 *(um);                     % the lower working depth
ez = 60 *(um);                      % the upper working depth

%
inputdir = 'C:\Users\motorsgroup\Downloads\20200416\HM06\SOC\Motile\T3h';   % the directory path of the time-series movie
outputdir = 'C:\Users\motorsgroup\Downloads\20200416\HM06\SOC\Motile\T3h';  % save tracking results

%% Preallocating 1vavriables and functions
listing = dir([inputdir,'\*.tif']);

% load depth library
load('C:\Users\motorsgroup\Downloads\20200416\HM06\SOC\Lib\Library.mat')
load('C:\Users\motorsgroup\Downloads\20200416\HM06\SOC\Lib\Library_Pos.mat')

% constrain library depth range
Library(:,:,Library_Pos<sz) = [];
Library_Pos(Library_Pos<sz) = [];
Library(:,:,Library_Pos>ez) = [];
Library_Pos(Library_Pos>ez) = [];

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
[~,msg] = mkdir([outputdir,'\Analysis']);

%% Human-made define target
nSeeds = 0;
listSeeds = [];
for nFile = 1:2%length(listing)
    inputfile = [inputdir,'\',listing(nFile).name];
    origina = double(imread(inputfile,1));
    
    figure(2), clf(gcf), imshow(origina,[]), title([num2str(nFile),' / ', num2str(length(listing))])
    nCounts = str2double(cell2mat(inputdlg( 'Enter number:')));
    if (nCounts>0)
        for nCount = 1:nCounts
            nSeeds = nSeeds+1;
            h = drawrectangle('Label',['Target',num2str(nSeeds)],'DrawingArea',[1 1 flip(size(origina))],'Position',[1, 1, search_Box_xy, search_Box_xy]);
%             h = imrect(gca,[1, 1, search_Box_xy, search_Box_xy]);
            wait(h);
            position = h.Position;
            sx = round(position(1));
            sy = round(position(2));
            ex = sx+search_Box_xy-1;
            ey = sy+search_Box_xy-1;
            
            listSeeds(nSeeds,:) = [nFile, nCount, sx, sy, ex, ey];
        end
    end
end

%% Tracking
for nSeed = 1:nSeeds
    
    outputname = [strrep(listing(listSeeds(nSeed,1)).name,'.tif',''),'-',num2str(listSeeds(nSeed,2))];
    vSession = VideoWriter([outputdir,'\Analysis\',outputname,'.avi']);
    vSession.FrameRate = 33;
    open(vSession)
    
    Pos = [];
    figure(1), clf(gcf), set(gcf,'WindowStyle','docked'), tic
    for nFrame = sFrame:eFrame
        
        % read image
        inputfile = [inputdir,'\',listing(listSeeds(nSeed,1)).name];
        origina = double(imread(inputfile,nFrame));
        
        % label tracking seeds
        if (nFrame ==sFrame)
            sx = listSeeds(nSeed,3);
            sy = listSeeds(nSeed,4);
            ex = listSeeds(nSeed,5);
            ey = listSeeds(nSeed,6);
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
        checksum_xy = (sx>=1)+(sy>=1)+(ex<=size(origina,2))+(ey<=size(origina,1));
        checksum_z = (sLayer>=1)+(eLayer<=nLayers);
        if (checksum_xy<4) || (checksum_z<2)
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
        if (nFrame>1)
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
    
    % plot swimmer's trajectory
    figure(3), clf(gcf), set(gcf,'WindowStyle','docked')
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
    
end

%%
disp('Done.')
