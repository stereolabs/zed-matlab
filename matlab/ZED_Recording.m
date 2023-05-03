clc;
disp('========= ZED SDK PLUGIN =========');
disp('-- Record SVO File --');
close all;
clear mex; clear functions; clear all;

% initial parameter structure, the same as sl::InitParameters
% values as enum number, defines in : sl/defines.hpp
% or from https://www.stereolabs.com/docs/api/structsl_1_1InitParameters.html

InitParameters.camera_resolution = 6; %AUTO
InitParameters.camera_fps = 60;
InitParameters.system_units = 2; %METER
InitParameters.depth_mode = 0; %DEPTH_MODE_NONE
result = mexZED('open', InitParameters);

nb_frame_to_save = 1000;
if(strcmp(result,'SUCCESS'))
    RecordingParameters.video_filename = 'MySVO.svo';
    RecordingParameters.compression_mode = 2;
    result =   mexZED('enableRecording', RecordingParameters);
    if(strcmp(result,'SUCCESS'))
        disp('Start Recording');
        f = 0;
        while f < nb_frame_to_save
            % grab the current image
            result = mexZED('grab');
            if(strcmp(result,'SUCCESS'))
                % record the grabbed image
                f = f+1;
            end
        end
        disp('End of Recording');
    end
    mexZED('disableRecording');
end
mexZED('close')
disp('========= END =========');
clear mex;