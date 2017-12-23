% DEMO generate eye movie


% Setup the video save location and path params

% This path should be defined in the eyeModelSupport local hook
codeDirectory = '~/Documents/MATLAB/toolboxes/eyemodelSupport/code';

exportsDirectory = '~/Desktop/Blender_Simulate01';
pathParams.dataOutputDirFull = fullfile(exportsDirectory);
pathParams.dataSourceDirFull = fullfile(exportsDirectory);
pathParams.runName = 'Simulate01';
videoName = fullfile(exportsDirectory, [pathParams.runName '_gray.avi']);

% check or make a directory for output
if exist(exportsDirectory,'dir')==0
    mkdir(exportsDirectory);
end

sensorSize = [36 24]; %mm
aperture = 2.0;
fieldOfViewDEG = 45;
sceneResolution = [640 480];
focalLengthPX = (sceneResolution(1)/2) / tand(45); % this is how it is defined in blender, but I am not sure it is right
focalLengthMM = 18;
focalLengthPX = (sceneResolution(2)/2) / tand(45); % this is the focal length that works

pixelSizeMM = sensorSize./sceneResolution;


% We find that objects in the image are 25% larger (in pixels) than what we
% would calculate here. For now, we apply this fudge factor while we try to
% understand how to calculate pixel size correctly.
fudgeFactor = 1;
pixelSizeMM = pixelSizeMM ./ fudgeFactor;



% eye radius (hardcoded in the python function)
eyeRadiusMM = 14;
eyeRadiusPX = eyeRadiusMM/pixelSizeMM(1);

% The Blender coordinate system defines depth zero as being in the center
% of rotation of the eye. Therefore, the camera distance is the distance of
% the camera from the center of rotation of the eye.
cameraDistanceMM = 50;

% The camera is hard-coded in the python routine to focus upon (i.e.,
% target) the front of the center of the eye. Therefore, the scene distance
% is equal to the camera distance + the eye raidus. We compute the scene
% distance in mm and then convert to pixels here.
sceneDistancePX = (cameraDistanceMM-eyeRadiusMM)/pixelSizeMM(1);

% define azimuth and elevation timeseries
eleSteps = [-15:5:15];
aziSweeps = [-25:5:25];

allPupilAzi=[];
allPupilEle=[];
for ii = 1: length(eleSteps)
    allPupilEle = [allPupilEle eleSteps(ii)*ones(1,length(aziSweeps))];
    if mod(ii,2)
        allPupilAzi = [allPupilAzi aziSweeps];
    else
        allPupilAzi = [allPupilAzi fliplr(aziSweeps)];
    end
end

allPupilAzi=[];
allPupilEle=[];
for ii = 1: length(eleSteps)
    allPupilEle = [allPupilEle eleSteps(ii)*ones(1,length(aziSweeps))];
    if mod(ii,2)
        allPupilAzi = [allPupilAzi aziSweeps];
    else
        allPupilAzi = [allPupilAzi fliplr(aziSweeps)];
    end
end


% The eye is posed in the Blender model by defining the location in XYZ
% coordinate space of an eye target. The coordinate system that we use in
% the transparentTrack routines has the convention X = left/right; Y =
% down/up; Z = nearer / farther. The subsequent python routine transposes
% some of these dimensions to be appropriate for the Blender coordinate
% convention, but this need not concern us here.
%
% We calculate the target of the gaze using equations that describe the
% Fick axis rotation of the eye. In this system, the Y (down / up) position
% of the center of the pupil is indenpendent of the azimuthal rotation. The
% gazeTargetPositions that we generate correspond to a point that is
% directly on the surface of the eye in the center of the pupil.
gazeTargetPositionX = (eyeRadiusMM).*sind(allPupilAzi).*cosd(allPupilEle);
gazeTargetPositionY = (eyeRadiusMM).*sind(allPupilEle);
gazeTargetPositionZ = -(eyeRadiusMM).*cosd(allPupilAzi).*cosd(allPupilEle);

% Define the pupil radius in pixels and the degree of eye closedness
pupilRadiusMM = 2*ones(1,length(gazeTargetPositionX));
pupilRadiusPX = pupilRadiusMM/pixelSizeMM(1);
eyeClosedness = 0*ones(1,length(gazeTargetPositionX));


%% generate eye movie
%generateEyeMovie(codeDirectory,exportsDirectory,gazeTargetPositionX, gazeTargetPositionY, gazeTargetPositionZ, pupilRadiusMM, eyeClosedness, cameraDistanceMM)

% rename file
%movefile(fullfile(exportsDirectory,'pupil_movie.avi'),fullfile(exportsDirectory,[pathParams.runName '_gray.avi']));

%% Run the processing pipeline
% runVideoPipeline( pathParams, ...
%     'verbosity', 'full', 'useParallel',false, 'catchErrors', false,...
%     'pupilFrameMask', [60 60], 'maskBox', [0.9 0.9], 'pupilGammaCorrection',0.7, 'pupilRange', [30 80], ...
%     'overwriteControlFile',true, 'glintPatchRadius', 10, ...
%     'ellipseTransparentLB',[0,0, 20, 0, 0],...
%     'ellipseTransparentUB',[sceneResolution(1),sceneResolution(2), 20000, 1.0, pi],...
%     'candidateThetas',0:pi/16:2*pi,...    
%     'badFrameErrorThreshold', 8, ...
%     'sceneGeometryLB',[0, 0, (sceneDistancePX+eyeRadiusPX)*1.2, 25],'sceneGeometryUB',[640, 480, (sceneDistancePX+eyeRadiusPX)*1.2, 500],...
%     'skipStageByNumber',1);

% Load the result files into memory
load(fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_pupil.mat']));
load(fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_sceneGeometry.mat']));

%% Report how closely we have reconstructed the actual scene geometry
fprintf('Veridical scene geometry - eye center: [%0.1f, %0.1f, %0.1f], eye radius: %0.1f \n',sceneResolution(1)/2,  sceneResolution(2)/2, sceneDistancePX+eyeRadiusPX, eyeRadiusPX);
fprintf('Estimated scene geometry - eye center: [%0.1f, %0.1f, %0.1f], eye radius: %0.1f \n',sceneGeometry.eyeCenter.X, sceneGeometry.eyeCenter.Y, sceneGeometry.eyeCenter.Z, sceneGeometry.eyeRadius);

%% use the inverse projection function and the scene geometry to recover the angles

reconstructedPupilAzi =[];
reconstructedPupilEle = [];

nPoints = size(pupilData.ellipseParamsAreaSmoothed_mean,1);

for ii=1:nPoints
    worldPoints(ii,1) = eyeRadiusPX*sind(allPupilAzi(ii))*cosd(allPupilEle(ii));
    worldPoints(ii,2) = eyeRadiusPX*sind(allPupilEle(ii));
    worldPoints(ii,3) = eyeRadiusPX*cosd(allPupilAzi(ii))*cosd(allPupilEle(ii));
    
    imagePoints(ii,1)=pupilData.ellipseParamsUnconstrained_mean(ii,1);
    imagePoints(ii,2)=pupilData.ellipseParamsUnconstrained_mean(ii,2);
end

worldPoints1=[worldPoints, ones(nPoints,1)];
imagePoints1=[imagePoints, ones(nPoints,1)];

intrinsicCameraMatrix = [focalLengthPX 0 sceneResolution(1)/2; 0 focalLengthPX sceneResolution(2)/2; 0 0 1];
[extrinsicRotationMatrix, extrinsicTranslationVector, worldPointsReconstructed]=efficient_pnp_gauss(worldPoints1, imagePoints1, intrinsicCameraMatrix);


%% NEXT STEP --
% Create a set of 5, 3D world coordinate points that are on the border of a
% pupil that is on a tilted plane for the azi / ele of the eye. Project
% these to the image plane. Fit the 5 points with an ellipse (5 points
% are needed to uniquely specify the ellipse). Identify the center of that
% fitted ellipse. This is the predicted location of the projection of that
% pupil in the image plane. Can also project the center of the circle and
% see where that lands.



% Reconstruct the points on the image plane



projectionMatrix=intrinsicCameraMatrix*[extrinsicRotationMatrix,extrinsicTranslationVector];
imagePointsReconstructedUnscaled=(projectionMatrix*worldPoints1')';

imagePointsReconstructed=zeros(nPoints,2);
imagePointsReconstructed(:,1)=imagePointsReconstructedUnscaled(:,1)./imagePointsReconstructedUnscaled(:,3);
imagePointsReconstructed(:,2)=imagePointsReconstructedUnscaled(:,2)./imagePointsReconstructedUnscaled(:,3);



reconstructionRMSE = sqrt(sum(sqrt((imagePointsReconstructed(:,1)-imagePoints(:,1)).^2 + (imagePointsReconstructed(:,2)-imagePoints(:,2)).^2)));

figure
plot(imagePoints(:,1),imagePoints(:,2),'o','MarkerEdgeColor',[.7 .7 .7],...
    'MarkerFaceColor',[.7 .7 .7])
hold on
plot(imagePointsReconstructed(:,1),imagePointsReconstructed(:,2),'.r')

figure
plot3(worldPoints(:,1),worldPoints(:,2),worldPoints(:,3),'o','MarkerEdgeColor',[.7 .7 .7],...
    'MarkerFaceColor',[.7 .7 .7])
hold on
plot3(worldPointsReconstructed(:,1),worldPointsReconstructed(:,2),worldPointsReconstructed(:,3),'.r')

