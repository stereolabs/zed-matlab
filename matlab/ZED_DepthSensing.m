clc;
disp('========= ZED SDK PLUGIN =========');
disp('-- Retrieve Images and Depth and compute Depth Histogram --');
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
    requested_depth_size = [720 404];

    % init depth histogram
    binranges = 0.5:0.25:InitParameters.depth_maximum_distance;

    % Create Figure and wait for keyboard interruption to quit
    f = figure('name','ZED SDK : Images and Depth','NumberTitle','off','KeyPressFcn',@(obj,evt) 0);
    % Setup runtime parameters
    RuntimeParameters.sensing_mode = 0; % STANDARD sensing mode

    ok = 1;
    % loop over frames, till Esc is pressed
    while (ok ~= 27)
        % grab the current image and compute the depth
        result = mexZED('grab', RuntimeParameters);        
        if(strcmp(result,'SUCCESS'))
            % retrieve letf image
            image_left = mexZED('retrieveImage', 0); %left

            % retrieve depth as a normalized image
            image_depth = mexZED('retrieveImage', 9); %depth
            % retrieve the real depth, resized
            depth = mexZED('retrieveMeasure', 1, requested_depth_size(1), requested_depth_size(2)); %depth

            % display
            subplot(2,2,1)
            imshow(image_left);
            title('Image Left')
            subplot(2,2,2)
            imshow(image_depth);
            title('Depth')
            subplot(2,2,3:4)
            % Compute the depth histogram
            val_ = find(isfinite(depth(:))); % handle wrong depth values
            depth_v = depth(val_);
            [bincounts] = histc(depth_v(:),binranges);
            bar(binranges,bincounts,'histc')
            title('Depth histogram')
            xlabel('meters')
            
            % redraw figure
            drawnow;
            % check for interrupts
            ok = uint8(get(f,'CurrentCharacter'));
            if(~length(ok))
                ok=0;
            end
        end
    end
end

close(f)
% Make sure to call this function to free the memory before use this again
mexZED('close')
disp('========= END =========');
clear mex;
