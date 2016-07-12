clc;
disp('=========SL_ZED_WITH_MATLAB -- Tracking=========');
close all;
clear mex; clear functions; clear all;

% SVO playback
% path = '../MySVO.svo';
% mexZED('create', path)

% Live mode
mexZED('create', 720, 60);

% parameter struct, the same as sl::zed::InitParams
% values as enum number, defines in : include/zed/utils/GlobalDefine.hpp
% 1: true, 0: false for boolean
param.unit = 1; % in this sample we use METER
param.mode = 2;
result = mexZED('init', param)

if(strcmp(result,'SUCCESS'))
 
    %enable Tracking
    mexZED('enableTracking');
    
    Display2D = 1;
    
    %Store tracking informations
    PositionArray = [];

    % Create Figure and wait for keyboard interruption to quit
    f = figure('name','ZED SDK','keypressfcn','close');
    ok = 1;
    % loop over frames
    while ok

        % grab the current image and compute the depth
        mexZED('grab', 'STANDARD')

        % retrieve letf image
        image_l = mexZED('retrieveImage', 'left');    
        % retrieve camera Path
        position = mexZED('getPosition');
        
        % display
        subplot(3,1,1)
        imshow(image_l);
        title('Image Left')
      
        % Draw 2D/3D ZED position
        subplot(3,1,2:3)
        Tx = position(1,4);
        Ty = position(2,4);
        Tz = position(3,4);
        PositionArray = [PositionArray; Tx Ty Tz];
        if Display2D
            plot(PositionArray(:,1), PositionArray(:,3),'-b');
            xlabel('Tx (M)');
            ylabel('Tz (M)');
        else
            plot3(PositionArray(:,1), PositionArray(:,2), PositionArray(:,3), '-r');
            xlabel('Tx (M)');
            ylabel('Ty (M)');
            zlabel('Tz (M)');
        end
        grid on;
        title('Odometry')

         drawnow; %this checks for interrupts
         ok = ishandle(f); %does the figure still exist
    end
end

% Make sure to call this function to free the memory before use this again
mexZED('delete')
