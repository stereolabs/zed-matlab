clc;
disp('========= ZED SDK PLUGIN =========');
disp('-- Get 3D Point Cloud --');
close all;
clear mex; clear functions; clear all;

% initial parameter structure, the same as sl::InitParameters
% values as enum number, defines in : sl/defines.hpp
% or from https://www.stereolabs.com/docs/api/structsl_1_1InitParameters.html

InitParameters.camera_resolution = 6; %AUTO
InitParameters.coordinate_units = 2; %METER
InitParameters.depth_mode =  1; %PERFORMANCE
InitParameters.coordinate_system = 3; %COORDINATE_SYSTEM_RIGHT_HANDED_Z_UP
%InitParameters.svo_input_filename = '../mySVOfile.svo'; % Enable SVO playback

% DepthClamp value, maximum depth (in METER)
depth_max = 5;
InitParameters.depth_maximum_distance = depth_max; 
result = mexZED('open', InitParameters);

if(strcmp(result,'SUCCESS'))
    % Create Figure
    f = figure('name','ZED SDK : Point Cloud','NumberTitle','off','keypressfcn',@(obj,evt) 0);
    %create 2 sub figure
    ha1 = axes('Position',[0.05,0.7,0.9,0.25]);
    ha2 = axes('Position',[0.05,0.05,0.9,0.6]);
    axis([-depth_max, depth_max, 0, depth_max, -depth_max ,depth_max])
    xlabel('X');
    ylabel('Z');
    zlabel('Y');
    grid on;
    hold on;
    
    % init point cloud data
    requested_size = [128 72];
    nb_elem = requested_size(1) * requested_size(2);
    pt_X = zeros(requested_size);
    pt_Y = zeros(requested_size);
    pt_Z = zeros(requested_size);
    
    % Setup runtime parameters
    RuntimeParameters.sensing_mode = 0; % STANDARD sensing mode
    
    h = plot3(reshape(pt_X, 1,nb_elem), reshape( pt_Y, 1,nb_elem), reshape( pt_Z, 1,nb_elem), '.');
    
    key = 1;
    % loop over frames, till Esc is pressed
    while (key ~= 27)
        % grab the current image and compute the depth
        result = mexZED('grab', RuntimeParameters);
        if(strcmp(result,'SUCCESS'))
            % retrieve letf image
            image_left = mexZED('retrieveImage', 0); %left
            %displays it
            axes(ha1);
            imshow(image_left);

            % retrieve the point cloud, resized
            [pt_X, pt_Y, pt_Z] = mexZED('retrieveMeasure', 3, requested_size(1), requested_size(2)); %XYZ pointcloud

            %displays it
            axes(ha2);
            set(h,'XData',reshape(pt_X, 1,nb_elem))
            set(h,'YData',reshape(pt_Y, 1,nb_elem))
            set(h,'ZData',reshape(pt_Z, 1,nb_elem))  

            drawnow; %this checks for interrupts
            key = uint8(get(f,'CurrentCharacter'));
            if(~length(key))
                key=0;
            end
        end
    end
    close(f)
end

% Make sure to call this function to free the memory before use this again
mexZED('close')
disp('========= END =========');
clear mex;