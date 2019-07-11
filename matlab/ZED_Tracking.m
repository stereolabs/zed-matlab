clc;
disp('=========SL_ZED_WITH_MATLAB -- Tracking=========');
close all;
clear mex; clear functions; clear all;

mexZED('create');

% parameter struct, the same as sl::InitParameters
% values as enum number, defines in : sl/defines.hpp 
% or from https://www.stereolabs.com/developers/documentation/API/
% 1: true, 0: false for boolean

InitParameters.camera_resolution = 2; %HD720
InitParameters.camera_fps = 60;
InitParameters.coordinate_units = 2; %METER
InitParameters.depth_mode = 1; %PERFORMANCE
InitParameters.coordinate_system = 3; %COORDINATE_SYSTEM_RIGHT_HANDED_Z_UP
%InitParameters.svo_input_filename = '../mySVOfile.svo'; % Enable SVO playback
result = mexZED('open', InitParameters)

if(strcmp(result,'SUCCESS'))
    
    %enable Tracking
    TrackingParameters.enable_spatial_memory = 1;
    %TrackingParameters.initial_world_transform = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
    mexZED('enableTracking', TrackingParameters);
    
    % for tracking informations storage
    PositionArray = [];
    
    % Create Figure and wait for keyboard interruption to quit
    f = figure('name','ZED SDK : Positional Tracking','NumberTitle','off');
    %create 2 sub figure
    ha1 = axes('Position',[0.05,0.7,0.9,0.25]);
    ha2 = axes('Position',[0.05,0.05,0.9,0.6]);
    xlabel('Tx (M)');
    ylabel('Tz (M)');
    zlabel('Ty (M)');
    xlim(ha2, [-2 2]);
    ylim(ha2, [-2 2]);
    axis equal, grid on;
    hold on;   
    % init 3d display
    h = plot3(0,0,0, 'r');
    
    ok = 1;
    % loop over frames
    while ok        
        % grab the current image and compute the positional tracking
        result = mexZED('grab');
        if(strcmp(result,'SUCCESS'))
            % retrieve letf image
            image_left = mexZED('retrieveImage', 0); %left
            %displays it
            axes(ha1);
            imshow(image_left);

            % retrieve camera Path
            position = mexZED('getPosition');
            %stack positions
            PositionArray = [PositionArray; position(1,4) position(2,4) position(3,4)];

            % retrieve IMU Data
            IMUdata = mexZED('getIMUData');

            axes(ha2);
            set(h,'XData',PositionArray(:,1))
            set(h,'YData',PositionArray(:,2))
            set(h,'ZData',PositionArray(:,3))

            drawnow; %this checks for interrupts
            ok = ishandle(f); %does the figure still exist
        end
    end
end

% Make sure to call this function to free the memory before use this again
mexZED('close')
clear mex;