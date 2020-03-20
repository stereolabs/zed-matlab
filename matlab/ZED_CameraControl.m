clc;
disp('========= ZED SDK PLUGIN =========');
disp('-- Retrieve Images and change Video Settings --');
close all;
clear mex; clear functions; clear all;

% initial parameter structure, the same as sl::InitParameters
% values as enum number, defines in : sl/defines.hpp
% or from https://www.stereolabs.com/docs/api/structsl_1_1InitParameters.html

InitParameters.camera_resolution = 2; %HD720
InitParameters.camera_fps = 60;
InitParameters.depth_mode = 0; %NONE % no depth computation needed here
result = mexZED('open', InitParameters);

if(strcmp(result,'SUCCESS')) % the Camera is open

    % basic informations
    camInfo = mexZED('getCameraInformation');
    image_size = [camInfo.left_cam.width camInfo.left_cam.height] 
    
    roi = [0,0,image_size(1), image_size(2)];

    % Create Figure and wait for keyboard interruption to quit
    f = figure('name','ZED SDK : Images and Depth','NumberTitle','off','KeyPressFcn',@(obj,evt) 0);
    
    key = 1;
    % loop over frames, till Esc is pressed
    while (key ~= 27)
        % grab the current image and compute the depth
        result = mexZED('grab');        
        if(strcmp(result,'SUCCESS'))
            % retrieve letf image
            image_left = mexZED('retrieveImage', 0); %left
            % retrieve right image
            image_right = mexZED('retrieveImage', 1); %right

            % image timestamp
            im_ts = mexZED('getTimestamp', 0);

            % display
            subplot(1,2,1)
            imshow(image_left);
            if(roi(3) < image_size(1)) % display Gain/Exposure ROI if neeeded
                rectangle('Position',roi,'EdgeColor', 'r','LineWidth', 1,'LineStyle','-');
            end
            title('Image Left')
            subplot(1,2,2)
            imshow(image_right); 
            if(roi(3) < image_size(1)) % display Gain/Exposure ROI if neeeded
                rectangle('Position',roi,'EdgeColor', 'r','LineWidth', 1,'LineStyle','-');
            end
            title('Image Right')
            
            % redraw figure
            drawnow;
            % check for interrupts
            key = uint8(get(f,'CurrentCharacter'));
            if(~length(key))
                key=0;
            else
                ask_plus = key == 45; % press '+'
                ask_minus = key == 43; % press '-'
                if(ask_plus | ask_minus)
                    % get current Camera brightness value
                    brightness = mexZED('getCameraSettings','brightness');
                    if(ask_plus & (brightness>0)) % decrease value
                        brightness = brightness - 1;
                    end
                    if(ask_minus & (brightness<8)) % increase value
                        brightness = brightness + 1;
                    end
                    brightness
                    % set the new value
                    mexZED('setCameraSettings', 'brightness', brightness);
                end
                
                if(key == 100) % press 'd' to reset to the default value
                    disp('reset to default');
                    mexZED('setCameraSettings', 'brightness', -1); % set auto value
                    roi = [0,0,image_size(1), image_size(2)];
                    mexZED('setCameraSettings', 'aec_agc_roi', roi, 2, 1); % set auto Gain/Exposure on full image
                end
                
                if(key == 114) % press 'r'  to use the auto Gain/Exposure on a defineded ROI
                    roi = [350, 250, 250, 125];
                    mexZED('setCameraSettings', 'aec_agc_roi', roi);
                end
            end            
            set(f,'CurrentCharacter','0'); % reset pressed key
        end
    end
end

close(f)
% Make sure to call this function to free the memory before use this again
mexZED('close')
disp('========= END =========');
clear mex;
