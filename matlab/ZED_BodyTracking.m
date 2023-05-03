clc;
disp('========= ZED SDK PLUGIN =========');
disp('-- Body Tracking --');
close all;
clear mex; clear functions; clear all;


% initial parameter structure, the same as sl::InitParameters
% values as enum number, defines in : sl/defines.hpp
% or from https://www.stereolabs.com/docs/api/structsl_1_1InitParameters.html

InitParameters.camera_resolution = 6; %AUTO
result = mexZED('open', InitParameters);

clrs = jet(6);

if(strcmp(result,'SUCCESS')) % the Camera is open    
	% basic informations
    camInfo = mexZED('getCameraInformation');    
    if(camInfo.model ==2) % Need a ZED2
        
        mexZED('enablePositionalTracking');
        
        BodyTrackingParameters.detection_model = 0; % HUMAN_BODY_FAST
        BodyTrackingParameters.body_format = 2; % BODY_38        
        BodyTrackingParameters.enable_tracking = 1; % enable tracking
        BodyTrackingParameters.enable_body_fitting = 1; %enable fiting 
        mexZED('enableBodyTracking', BodyTrackingParameters);
        BodyTrackingRuntimeParameters.detection_confidence_threshold = 50;
        BodyTrackingRuntimeParameters.minimum_keypoints_threshold = 5;
        
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
                
                bodies = mexZED('retrieveBodies', BodyTrackingRuntimeParameters);
                
                % display
                imshow(image_left);
                
                for b = 1: length(bodies.body_list)
                    
                    clr = clrs(mod(max(0, bodies.body_list(b).id), 6)+1, :);
                    
                    keypoints2d =  bodies.body_list(b).keypoint_2d;
                    for kp = 1: length(keypoints2d)
                        kp2d = keypoints2d(kp);
                        if(kp2d.u >= 0) % discard undetected keypoints
                            rectangle('Position', [(kp2d.u * ratio(1) - 2), (kp2d.v * ratio(2) - 2), 4, 4]...
                                ,'Curvature',[1,1], 'LineWidth',3, 'EdgeColor', clr);
                        end
                    end
                    
                    bb2d = bodies.body_list(b).bounding_box_2d;
                    rectangle('Position', [bb2d(1).u * ratio(1), bb2d(1).v * ratio(2), (bb2d(3).u - bb2d(1).u) * ratio(1), (bb2d(3).v - bb2d(1).v) * ratio(2)]...
                        ,'Curvature',0.05, 'LineWidth',3, 'EdgeColor', clr);
                    
                    ['Id: ' num2str(bodies.body_list(b).id) ' Conf: ', num2str(bodies.body_list(b).confidence)]
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
