clc;
disp('========= ZED SDK PLUGIN =========');
disp('-- Object Detection --');
close all;
clear mex; clear functions; clear all;


% initial parameter structure, the same as sl::InitParameters
% values as enum number, defines in : sl/defines.hpp
% or from https://www.stereolabs.com/docs/api/structsl_1_1InitParameters.html

InitParameters.camera_resolution = 0; %HD2K
result = mexZED('open', InitParameters);

clrs = jet(6);

if(strcmp(result,'SUCCESS')) % the Camera is open    
	% basic informations
    camInfo = mexZED('getCameraInformation');    
    if(camInfo.model ==2) % Need a ZED2
        
        mexZED('enablePositionalTracking');
        
        ObjectDetectionParameters.detection_model = 0; % MULTI_CLASS_BOX
        ObjectDetectionParameters.enable_tracking = 1; % enable Object tracking
        mexZED('enableObjectDetection', ObjectDetectionParameters);        
        ObjectDetectionRuntimeParameters.detection_confidence_threshold = 35
        
        % Create Figure and wait for keyboard interruption to quit
        f = figure('name','ZED SDK: Detection','NumberTitle','off', 'keypressfcn',@(obj,evt) 0);
        
        display_size = [720 404];
        ratio = [display_size(1) / camInfo.left_cam.width display_size(2) / camInfo.left_cam.height];
        
        key = 1;
        % loop over frames, till Esc is pressed
        while (key ~= 27)
            % grab the current image and compute the depth
            result = mexZED('grab');
            if(strcmp(result,'SUCCESS'))
                % retrieve letf image
                image_left = mexZED('retrieveImage', 0, display_size(1), display_size(2)); %left
                
                objs = mexZED('retrieveObjects', ObjectDetectionRuntimeParameters);
                
                % display
                imshow(image_left);
                
                for o = 1: length(objs.object_list)
                    bb2d = objs.object_list(o).bounding_box_2d;
                    clr = clrs(mod(max(0, objs.object_list(o).id), 6)+1, :);
                    rectangle('Position', [bb2d(1).u * ratio(1), bb2d(1).v * ratio(2), (bb2d(3).u - bb2d(1).u) * ratio(1), (bb2d(3).v - bb2d(1).v) * ratio(2)]...
                        ,'Curvature',0.05, 'LineWidth',3, 'EdgeColor', clr);
                end
                
                drawnow;  %this checks for interrupts  
                key = uint8(get(f,'CurrentCharacter'));
                if(~length(key))
                    key=0;
                end
            end
        end
        close(f)
    end
    mexZED('close')
    disp('========= END =========');
    clear mex;
end
