clc;
disp('=========SL_ZED_WITH_MATLAB=========');
close all;
clear mex; clear functions; clear all;

% SVO playback
% path = '../MySVO.svo';
% result = mexZED('init', path, 'quality')

% Live mode
result = mexZED('init', 720, 'quality')

if(strcmp(result,'SUCCESS'))

    size = mexZED('getImageSize');

    % Set Confidence Threshold
    mexZED('setConfidenceThreshold', 95);

    % Define maximum depth (in mm)
    maxDepth = 10 * 1000;
    binranges = 0:100:maxDepth ;
    mexZED('setDepthClampValue',maxDepth);

    % Get number of frames (if SVO)
    nbFrame = mexZED('getSVONumberOfFrames');

    % Get cameras parameters
    params = mexZED('getCameraParameters');

    % Create Figure and wait for keyboard interruption to quit
    f = figure('keypressfcn','close','windowstyle','modal');
    ok = 1;
    % loop over frames
    while ok

        % grab the current image and compute the depth
        mexZED('grab', 'raw')

        % retrieve letf image
        image_l = mexZED('retrieveImage', 'left');
        % retrieve right image
        image_r = mexZED('retrieveImage', 'right');

        % retrieve depth as a normalized image
        depth_im = mexZED('normalizeMeasure', 'depth');
        % retrieve the real depth
        depth = mexZED('retrieveMeasure', 'depth');

        % display
        subplot(2,2,1)
        imshow(image_l);
        title('Image Left')
        subplot(2,2,2)
        imshow(image_r);
        title('Image Right')
        subplot(2,2,3)
        imshow(depth_im);
        title('Depth')
        subplot(2,2,4)

        % Compute the depth histogram
        val_ = find(depth(:) > 0); % valid depth is above 0
        depth_v = depth(val_);
        [bincounts] = histc(depth_v(:),binranges);
        bar(binranges,bincounts,'histc')
        title('Depth histogram')

        drawnow; %this checks for interrupts
        ok = ishandle(f); %does the figure still exist
    end
end

% Make sure to call this function to free the memory before use this again
mexZED('delete')
