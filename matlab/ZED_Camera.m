clc;
disp('=========SL_ZED_WITH_MATLAB -- Basic=========');
close all;
clear mex; clear functions; clear all;

mexZED('create');

% parameter struct, the same as sl::InitParameters
% values as enum number, defines in : sl/GlobalDefine.hpp
% or from https://www.stereolabs.com/developers/documentation/API/
% 1: true, 0: false for boolean

InitParameters.camera_resolution = 2; %HD720
InitParameters.camera_fps = 60;
InitParameters.system_units = 2; %METER
InitParameters.depth_mode = 3; %QUALITY
%param.svo_filename = '../mySVOfile.svo'; % Enable SVO playback
result = mexZED('open', InitParameters)

if(strcmp(result,'Error code:  Success'))
    
    size = mexZED('getResolution')
    
    % Set Confidence Threshold
    mexZED('setConfidenceThreshold', 98);
    
    % Define maximum depth (in METER)
    maxDepth = 10;
    binranges = 0:0.1:maxDepth ;
    mexZED('setDepthMaxRangeValue',maxDepth);
    
    % Get number of frames (if SVO)
    nbFrame = mexZED('getSVONumberOfFrames');
    
    camInfo = mexZED('getCameraInformation');
    
    % Create Figure and wait for keyboard interruption to quit
    f = figure('name','ZED SDK : Images and Depth','NumberTitle','off','keypressfcn','close');
    ok = 1;
    % loop over frames
    while ok
        
        % grab the current image and compute the depth
        RuntimeParameters.sensing_mode = 0; %STANDARD
        RuntimeParameters.enable_depth = 1;
        RuntimeParameters.enable_point_cloud = 0;
        mexZED('grab', RuntimeParameters)
        
        % retrieve letf image
        image_left = mexZED('retrieveImage', 0); %left
        % retrieve right image
        image_right = mexZED('retrieveImage', 1); %right
        
        % retrieve depth as a normalized image
        image_depth = mexZED('retrieveImage', 9); %depth
        % retrieve the real depth
        depth = mexZED('retrieveMeasure', 1); %depth
        
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

% Make sure to call this function to free the memory before use this again
mexZED('close')
