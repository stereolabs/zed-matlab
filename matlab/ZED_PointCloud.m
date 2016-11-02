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

% DepthClamp value
Depth_max = 5;
% Step for mesh display
data_Step = 10;

if(strcmp(result,'SUCCESS'))
    mexZED('setDepthClampValue', Depth_max)
    % Create Figure
    f = figure('name','ZED SDK : Point Cloud','NumberTitle','off');
    %create 2 sub figure
    ha1 = axes('Position',[0.05,0.7,0.9,0.25]);
    ha2 = axes('Position',[0.05,0.05,0.9,0.6]);
    axis([-Depth_max, Depth_max, 0, Depth_max, -Depth_max ,Depth_max])
    xlabel('X');
    ylabel('Z');
    zlabel('Y');
    grid on;
    hold on;
    
    % init mesh data
    size = mexZED('getImageSize')
    pt_X = zeros(size(1),size(2));
    pt_Y = zeros(size(1),size(2));
    pt_Z = zeros(size(1),size(2));
    h = mesh(pt_X(:,:),  pt_Z(:,:), pt_Y(:,:))
    
    ok = 1;
    % loop over frames
    while ok
        % grab the current image and compute the depth
        mexZED('grab', 'STANDARD')
        
        % retrieve letf image
        image_l = mexZED('retrieveImage', 'left');
        %displays it
        axes(ha1);
        imshow(image_l);
        
        % retrieve the point cloud
        [pt_X, pt_Y, pt_Z] = mexZED('retrieveMeasure', 'XYZ');
        
        %displays it
        axes(ha2);
        set(h,'XData',pt_X(1:data_Step:end,1:data_Step:end))
        set(h,'YData',pt_Z(1:data_Step:end,1:data_Step:end))
        set(h,'ZData',pt_Y(1:data_Step:end,1:data_Step:end))
        
        drawnow; %this checks for interrupts
        ok = ishandle(f); %does the figure still exist
    end
end

% Make sure to call this function to free the memory before use this again
mexZED('delete')
