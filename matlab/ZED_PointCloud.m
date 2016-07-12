clc;
disp('=========SL_ZED_WITH_MATLAB -- Point Cloud=========');
close all;
clear mex; clear functions; clear all;

% SVO playback
% path = '../MySVO.svo';
% mexZED('create', path)

% Live mode
mexZED('create', 720, 60);

% parameter struct, the same as sl::zed::InitParams
% values as enum number, defines in : include/zed/utils/GlobalDefine.hpp
% 1: true, 0: false for boolean
param.unit = 1; % in this sample we use METER
param.mode = 2;
param.coordinate = 2;
result = mexZED('init', param)

if(strcmp(result,'SUCCESS'))

    s = mexZED('getImageSize')

    % Get cameras parameters
    params = mexZED('getCameraParameters')

    mexZED('setDepthClampValue', 5)

    % Create Figure and wait for keyboard interruption to quit
    f =figure('name','ZED SDK','keypressfcn','close', 'Color',[1 1 1]);

    az = 10;
    el = 25;

    ok = 1;
    % loop over frames
    while ok

        % grab the current image and compute the depth
        mexZED('grab', 'STANDARD')

        % retrieve the point cloud 
        [pt_X, pt_Y, pt_Z] = mexZED('retrieveMeasure', 'XYZ');            
        mesh(pt_X(:,:),  pt_Z(:,:), pt_Y(:,:))
        grid on;

        view(az, el);

        drawnow; %this checks for interrupts
        ok = ishandle(f); %does the figure still exist
    end
end

% Make sure to call this function to free the memory before use this again
mexZED('delete')
