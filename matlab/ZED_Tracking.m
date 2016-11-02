clc;
disp('=========SL_ZED_WITH_MATLAB -- Tracking=========');
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
result = mexZED('init', param)

if(strcmp(result,'SUCCESS'))
    
    %enable Tracking
    mexZED('enableTracking');
    
    % for tracking informations storage
    PositionArray = [];
    
    % Create Figure and wait for keyboard interruption to quit
    f = figure('name','ZED SDK : Positional Tracking','NumberTitle','off');
    %create 2 sub figure
    ha1 = axes('Position',[0.05,0.7,0.9,0.25]);
    ha2 = axes('Position',[0.05,0.05,0.9,0.6]);
    xlabel('Tx (M)');
    ylabel('Tz (M)');
    zlabel('Ty (M)');
    axis equal, grid on;
    hold on;
    
    % init 3d display
    h = plot3(0,0,0, 'r');
    
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
        
        % retrieve camera Path
        position = mexZED('getPosition');
        %stack positions
        PositionArray = [PositionArray; position(1,4) position(3,4) position(2,4)];
        
        axes(ha2);
        set(h,'XData',PositionArray(:,1))
        set(h,'YData',PositionArray(:,2))
        set(h,'ZData',PositionArray(:,3))
        
        drawnow; %this checks for interrupts
        ok = ishandle(f); %does the figure still exist
    end
end

% Make sure to call this function to free the memory before use this again
mexZED('delete')
