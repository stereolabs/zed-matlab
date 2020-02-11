clc;
disp('========= ZED SDK PLUGIN =========');
disp('-- Object Detection --');
close all;
clear mex; clear functions; clear all;

% initial parameter structure, the same as sl::InitParameters
% values as enum number, defines in : sl/defines.hpp
% or from https://www.stereolabs.com/docs/api/structsl_1_1InitParameters.html

InitParameters.camera_resolution = 2; %HD720
InitParameters.camera_fps = 30;
InitParameters.coordinate_units = 2; %METER
InitParameters.depth_mode = 3; %ULTRA
InitParameters.depth_maximum_distance = 10;% Define maximum depth (in METER)
InitParameters.coordinate_system = 3; %COORDINATE_SYSTEM_RIGHT_HANDED_Z_UP
result = mexZED('open', InitParameters);

if(strcmp(result,'SUCCESS')) % the Camera is open    
    % basic informations
    camInfo = mexZED('getCameraInformation');    
    if(camInfo.model ==2) % Need a ZED2
        
        mexZED('enablePositionalTracking');
        mexZED('enableObjectDetection');
        
        % Create Figure and wait for keyboard interruption to quit
        f = figure('name','ZED SDK: Detection','NumberTitle','off','keypressfcn','close');
                
        ok = 1;
        % loop over frames
        while ok        
            % grab the current image and compute the depth
            result = mexZED('grab');        
            if(strcmp(result,'SUCCESS'))
                % retrieve letf image
                image_left = mexZED('retrieveImage', 0); %left
                objs = mexZED('retrieveObjects', 30);
                
                % display
                imshow(image_left);
                title('Detection Left');
                
                for o = 1: length(objs.object_list)
                    bb2d = objs.object_list(o).bounding_box_2d;
                    rectangle('Position', [bb2d(1).u, bb2d(1).v, bb2d(3).u - bb2d(1).u, bb2d(3).v - bb2d(1).v]...
                    ,'Curvature',0.2, 'LineWidth',3);
                end
                                
                drawnow; %this checks for interrupts
                ok = ishandle(f); %does the figure still exist
            end
         end
     end
end

% Make sure to call this function to free the memory before use this again
mexZED('close')
disp('========= END =========');
clear mex;