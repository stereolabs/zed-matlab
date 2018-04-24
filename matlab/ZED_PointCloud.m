clc;
disp('=========SL_ZED_WITH_MATLAB -- Point Cloud=========');
close all;
clear mex; clear functions; clear all;

mexZED('create');

% parameter struct, the same as sl::InitParameters
% values as enum number, defines in : sl/defines.hpp 
% or from https://www.stereolabs.com/developers/documentation/API/
% 1: true, 0: false for boolean

InitParameters.camera_resolution = 2; %HD720
InitParameters.camera_fps = 60;
InitParameters.system_units = 2; %METER
InitParameters.depth_mode =  1; %PERFORMANCE
InitParameters.coordinate_system = 3; %COORDINATE_SYSTEM_RIGHT_HANDED_Z_UP
%param.svo_input_filename = '../mySVOfile.svo'; % Enable SVO playback
result = mexZED('open', InitParameters)

% DepthClamp value
depth_max = 5;
% Step for mesh display
data_Step = 10;

if(strcmp(result,'SUCCESS'))
    mexZED('setDepthMaxRangeValue', depth_max)
    % Create Figure
    f = figure('name','ZED SDK : Point Cloud','NumberTitle','off');
    %create 2 sub figure
    ha1 = axes('Position',[0.05,0.7,0.9,0.25]);
    ha2 = axes('Position',[0.05,0.05,0.9,0.6]);
    axis([-depth_max, depth_max, 0, depth_max, -depth_max ,depth_max])
    xlabel('X');
    ylabel('Z');
    zlabel('Y');
    grid on;
    hold on;
    
    % init mesh data
    image_size = mexZED('getResolution')    
    requested_mesh_size = [128 72];
    pt_X = zeros(requested_mesh_size);
    pt_Y = zeros(requested_mesh_size);
    pt_Z = zeros(requested_mesh_size);
    
    nb_elem = requested_mesh_size(1) * requested_mesh_size(2);
    h = plot3(reshape(pt_X, 1,nb_elem), reshape( pt_Y, 1,nb_elem), reshape( pt_Z, 1,nb_elem), '.');
    
    % Setup runtime parameters
    RuntimeParameters.sensing_mode = 0; % STANDARD sensing mode
    RuntimeParameters.enable_depth = 1;
    RuntimeParameters.enable_point_cloud = 1;
    
    ok = 1;
    % loop over frames
    while ok
        % grab the current image and compute the depth
        mexZED('grab', RuntimeParameters)
        
        % retrieve letf image
        image_left = mexZED('retrieveImage', 0); %left
        %displays it
        axes(ha1);
        imshow(image_left);
        
        % retrieve the point cloud, resized
        [pt_X, pt_Y, pt_Z] = mexZED('retrieveMeasure', 3, requested_mesh_size(1), requested_mesh_size(2)); %XYZ pointcloud
                
        %displays it
        axes(ha2);        
        set(h,'XData',reshape(pt_X, 1,nb_elem))
        set(h,'YData',reshape(pt_Y, 1,nb_elem))
        set(h,'ZData',reshape(pt_Z, 1,nb_elem))
               
        drawnow; %this checks for interrupts
        ok = ishandle(f); %does the figure still exist
    end
end

% Make sure to call this function to free the memory before use this again
mexZED('close')
clear mex;