clc;
disp('========= ZED SDK PLUGIN =========');
disp('-- Retrieve Images and Depth --');
close all;
clear mex; clear functions; clear all;

% initial parameter structure, the same as sl::InitParameters
% values as enum number, defines in : sl/defines.hpp
% or from https://www.stereolabs.com/docs/api/structsl_1_1InitParameters.html

InitParameters.camera_resolution = 2; %HD720
InitParameters.camera_fps = 60;
InitParameters.coordinate_units = 2; %METER
InitParameters.depth_mode = 1; %PERFORMANCE
%InitParameters.svo_input_filename = '../mySVOfile.svo'; % Enable SVO playback
InitParameters.depth_maximum_distance = 7;% Define maximum depth (in METER)
result = mexZED('open', InitParameters);

if(strcmp(result,'SUCCESS')) % the Camera is open
    
    % basic informations
    camInfo = mexZED('getCameraInformation');
    image_size = [camInfo.left_cam.width camInfo.left_cam.height] 
    
    requested_depth_size = [720 404];
    
    % init depth histogram
    binranges = 0.5:0.25:InitParameters.depth_maximum_distance;
    
    % (optional) Get number of frames (if SVO)
    nbFrame = mexZED('getSVONumberOfFrames');    
    
    % Create Figure and wait for keyboard interruption to quit
    f = figure('name','ZED SDK : Images and Depth','NumberTitle','off','keypressfcn','close');
        
    % Setup runtime parameters
    RuntimeParameters.sensing_mode = 0; % STANDARD sensing mode
    
    ok = 1;
    % loop over frames
    while ok        
        % grab the current image and compute the depth
        result = mexZED('grab', RuntimeParameters);        
        if(strcmp(result,'SUCCESS'))
            % retrieve letf image
            image_left = mexZED('retrieveImage', 0); %left
            % retrieve right image
            image_right = mexZED('retrieveImage', 1); %right
            
            % image timestamp
            im_ts = mexZED('getTimestamp', 0) 
                        
            % retrieve depth as a normalized image
            image_depth = mexZED('retrieveImage', 9); %depth
            % retrieve the real depth, resized
            depth = mexZED('retrieveMeasure', 1, requested_depth_size(1), requested_depth_size(2)); %depth

            % display
            subplot(2,2,1)
            imshow(image_left);
            title('Image Left')
            subplot(2,2,2)
            imshow(image_right);
            title('Image Right')
            subplot(2,2,3)
            imshow(image_depth);
            title('Depth')
            subplot(2,2,4)
            % Compute the depth histogram
            val_ = find(isfinite(depth(:))); % handle wrong depth values
            depth_v = depth(val_);
            [bincounts] = histc(depth_v(:),binranges);
            bar(binranges,bincounts,'histc')
            title('Depth histogram')
            xlabel('meters')

            drawnow; %this checks for interrupts
            ok = ishandle(f); %does the figure still exist
        end
    end
end

% Make sure to call this function to free the memory before use this again
mexZED('close')
disp('========= END =========');
clear mex;