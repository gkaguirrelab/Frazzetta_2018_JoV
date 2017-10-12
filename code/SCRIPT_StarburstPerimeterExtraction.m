% find perimeter with starburst

p.Results.nFrames  = 1000;
grayVideoName = '/Users/giulia/Desktop/ETDemo/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/tfMRI_RETINO_PA_run01_gray.avi';
videoOutFileName = '/Users/giulia/Desktop/ETDemo/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/tfMRI_RETINO_PA_run01_STARBURSTperim.avi';
ellipseStarburstFile = '/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/Frazzetta_201x_transparentTrack/ExperimentalData/3020rfMRI_REST_AP_run01_STARBURSTellipse.mat';
perimeterFile = '/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/Frazzetta_201x_transparentTrack/ExperimentalData/3020rfMRI_REST_AP_run01_STARBURSTperimeter.mat';

grayVideoName ='/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/TOME_processing/session1_restAndStructure/TOME_3018/040717/EyeTracking/rfMRI_REST_AP_run01_gray.avi';
videoOutFileName = '/Users/giulia/Desktop/TEST/3018rfMRI_REST_AP_run01_STARBURSTellipse.avi';
ellipseStarburstFile = '/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/Frazzetta_201x_transparentTrack/ExperimentalData/3018rfMRI_REST_AP_run01_STARBURSTellipse.mat';
perimeterFile = '/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/Frazzetta_201x_transparentTrack/ExperimentalData/3018rfMRI_REST_AP_run01_STARBURSTperimeter.mat';
%% Read video file into memory
% load pupilPerimeter
videoInObj = VideoReader(grayVideoName);
% get number of frames
if p.Results.nFrames == Inf
    nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
else
    nFrames = p.Results.nFrames;
end
% get video dimensions
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;
% initialize variable to hold the perimeter data
grayVideo = zeros(videoSizeY,videoSizeX,nFrames,'uint8');
% read the video into memory, adjusting gamma and local contrast
for ii = 1:nFrames
    thisFrame = readFrame(videoInObj);
    grayVideo(:,:,ii) = rgb2gray (thisFrame);
end
% close the video object
clear videoInObj


%%
% intialize perim data
perimeter.data = zeros(videoSizeY,videoSizeX,nFrames,'uint8');
ellipseStarburst = zeros(nFrames,5);
% open figure
frameFig = figure;
for ii = 1:nFrames
    
    % get the frame
    thisFrame = squeeze(grayVideo(:,:,ii));
    
    
    % starburst
    sigma = .5;                      % Standard deviation of image smoothing
    angle_delta = 1*pi/180;         % discretization step size (radians)
    cr_window_size=101;             % corneal reflection search window size (about [sx,sy] center)
    min_feature_candidates=30;      % minimum number of pupil feature candidates
    max_ransac_iterations=10000;    % maximum number of ransac iterations
    rays= 50;                        % number of rays to use to detect feature points
    
    %
    sx = videoSizeX/2;
    sy = videoSizeY/2;
    edge_thresh = 400;
    
    I = gaussian_smooth_image(thisFrame, sigma);
    [crx, cry, crar] = locate_corneal_reflection(I, sx, sy, cr_window_size);
    
    crr = fit_circle_radius_to_corneal_reflection(I, crx, cry, crar, angle_delta);
    crr = ceil(crr*2.5);
    
    I = remove_corneal_reflection(I, crx, cry, crr, angle_delta);
    cr_circle = [crx cry crr];
    
    % %     % adjust the gamma to make the glint patch less evident
    %     I = imadjust(I,[],[],1.2);
    
    [epx, epy] = starburst_pupil_contour_detection(I, sx, sy, edge_thresh, rays, min_feature_candidates);
    
    % get perimeter points in our format
    perimFrame = im2uint8(zeros(size(thisFrame)));
    perimFrame(sub2ind(size(perimFrame),epy,epx))=255;
    perimeter.data(:,:,ii) = perimFrame;
    
    % fit ellipse
    if isempty(epx) || isempty(epy)
        pupil_ellipse = [0 0 0 0 0];
        
    else
        
        [ellipse, inliers] = fit_ellipse_ransac(epx, epy, max_ransac_iterations);
        
        
        if ellipse(3) < 1 || ellipse(3) > size(I,2) || ellipse(4) < 1 || ellipse(4) > size(I,1)
            fprintf(1, 'Error! The ellipse center lies out of the image\n');
            pupil_ellipse = [0 0 0 0 0];
        else
            [pupil_ellipse] = fit_ellipse_model(I, ellipse, angle_delta);
        end
    end
    
    % save out starburst best ellipse
    ellipseStarburst(ii,:) = pupil_ellipse;
    
    
    % save out plot
    I = plot_ellipse_in_image(thisFrame, pupil_ellipse,[0 255 0]);
    
    imshow(I, 'Border', 'tight');
    hold on
    plot(epx, epy, '.')
    hold off
    
    tmp=getframe(frameFig);
    outputVideo(:,:,:,ii)=tmp.cdata;
    
end

%% save output

save(ellipseStarburstFile, 'ellipseStarburst');
save(perimeterFile,'perimeter')

% Create a color map
cmap = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
cmap(1,:)=[1 0 0];
cmap(2,:)=[0 1 0];
cmap(3,:)=[0 0 1];
cmap(4,:)=[1 1 0];
cmap(5,:)=[0 1 1];
cmap(6,:)=[1 0 1];

% write the outputVideo to file
videoOutObj = VideoWriter(videoOutFileName,'Indexed AVI');
videoOutObj.FrameRate = 60;
videoOutObj.Colormap = cmap;
open(videoOutObj);

% loop through the frames and save them
for ii=1:nFrames
    indexedFrame = rgb2ind(squeeze(outputVideo(:,:,:,ii)), cmap, 'nodither');
    writeVideo(videoOutObj,indexedFrame);
end
% close the videoObj
clear videoOutObj

