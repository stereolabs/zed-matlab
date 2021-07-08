clc;
disp('========= ZED SDK PLUGIN =========');
disp('-- ZED Camera Spatial Mapping --');
close all;
clear mex; clear functions; clear all;

% initial parameter structure, the same as sl::InitParameters
% values as enum number, defines in : sl/defines.hpp
% or from https://www.stereolabs.com/docs/api/structsl_1_1InitParameters.html

InitParameters.camera_resolution = 2; %HD720
InitParameters.camera_fps = 60;
InitParameters.coordinate_units = 2; %METER
InitParameters.depth_mode = 1; %PERFORMANCE
InitParameters.coordinate_system = 3; %COORDINATE_SYSTEM_RIGHT_HANDED_Z_UP
%InitParameters.svo_input_filename = '../mySVOfile.svo'; % Enable SVO playback
result = mexZED('open', InitParameters);

if(strcmp(result,'SUCCESS'))
        
    %enable Tracking
    PositionalTrackingParameters.enable_spatial_memory = 1;
    mexZED('enablePositionalTracking', PositionalTrackingParameters);
    
    %enable Spatial Mapping
    SpatialMappingParameters.map_type = 0;
    SpatialMappingParameters.range_meter = 5.;
    SpatialMappingParameters.resolution_meter = 0.08;
    mexZED('enableSpatialMapping', SpatialMappingParameters);
    
    f = figure('name','ZED SDK : Spatial Mapping','NumberTitle','off','keypressfcn',@(obj,evt) 0);
     
    key = 1;
    % loop over frames, till Esc is pressed
    while (key ~= 27)       
        % grab the current image and compute the positional tracking
        result = mexZED('grab');
        if(strcmp(result,'SUCCESS'))
            % retrieve letf image
            image_left = mexZED('retrieveImage', 0); %left
            %displays it
            imshow(image_left);
            
            drawnow; %this checks for interrupts
            key = uint8(get(f,'CurrentCharacter'));
            if(~length(key))
                key=0;
            end
        end
    end
    
    if(SpatialMappingParameters.map_type == 0) % Mesh
        [vertices, faces] = mexZED('extractWholeSpatialMap');
    else % Fused Point Cloud
        [vertices, colors] = mexZED('extractWholeSpatialMap');
    end
    close(f)
end

% Make sure to call this function to free the memory before use this again
mexZED('close')
disp('========= END =========');
clear mex;