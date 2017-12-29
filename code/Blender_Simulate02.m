% DEMO generate eye movie

close all
clear all

%% Set paths and file names
% This path should be defined in the eyeModelSupport local hook
codeDirectory = '~/Documents/MATLAB/toolboxes/eyemodelSupport/code';

% Where to save the results of our simulation
exportsDirectory = '~/Desktop/Blender_Simulate01';
pathParams.dataOutputDirFull = fullfile(exportsDirectory);
pathParams.dataSourceDirFull = fullfile(exportsDirectory);
pathParams.runName = 'Simulate01';
videoName = fullfile(exportsDirectory, [pathParams.runName '_gray.avi']);

% check or make a directory for output
if exist(exportsDirectory,'dir')==0
    mkdir(exportsDirectory);
end


%% Hard coded properties of the Blender camera
% The focal length is set in the routine eyeModel/__init__.py
% All intrinsic camera values are in units of pixels
sceneResolutionPx = [640 480];
opticalCenterPx = (sceneResolutionPx./2);
focalLengthPX = (sceneResolutionPx(1)/2.0) / tan(45*pi/180 / 2);
sceneGeometry.intrinsicCameraMatrix = ...
    [focalLengthPX 0 opticalCenterPx(1);...
     0 focalLengthPX opticalCenterPx(2);...
     0 0 1];

 
%% Hard coded aspect of the scene geometry
% The eyeRadius is set in the routine eyeModel/__init__.py
sceneGeometry.eyeRadius = 12;   % in mm


%% Variable aspects of scene geometry
% Variable aspect of the scene geometry. The Blender coordinate system
% defines depth zero as being at the center of rotation of the eye.
% Therefore, the depth position of the camera is relative to the center of
% rotation of the eye. The distance of the camera from the imaged scene,
% therefore, is equal to [cameraDepthPosition - eyeRadius]. We would like
% it to be the case that the extrinsicTranslationVector is estimated to
% have this scene distance value (50-12) for the Z dimension.
cameraDepthPositionMM = 50; % in mm

sceneGeometry.extrinsicTranslationVector = [0; 0; cameraDepthPositionMM - sceneGeometry.eyeRadius];
sceneGeometry.extrinsicRotationMatrix = [1 0 0; 0 -1 0; 0 0 -1];

%% Generate a set of gaze targets for the blender model

% We define the gaze target in terms of head-fixed azimuth and elevation
% rotations.
aziSweeps = [-25:25:25];
eleSteps = [-15:15:15];

pupilCenterAzimuths=[];
pupilCenterElevations=[];
for ii = 1: length(eleSteps)
    pupilCenterElevations = [pupilCenterElevations eleSteps(ii)*ones(1,length(aziSweeps))];
    if mod(ii,2)
        pupilCenterAzimuths = [pupilCenterAzimuths aziSweeps];
    else
        pupilCenterAzimuths = [pupilCenterAzimuths fliplr(aziSweeps)];
    end
end

nFrames = length(pupilCenterAzimuths);

% The eye is posed in the Blender model by defining the location in XYZ
% coordinate space of an eye target. Our sceneWorld coordinate system has
% the convention X = left/right; Y = down/up; Z = nearer / farther, with
% the origin corresponding to the front surface of the eye and the center
% of the pupil when the line that connects the center of rotation of the
% eye and the center of the pupil are normal to the image plane.
%
% Blender, however, uses a sceneWorld coordinate system corresponding to X
% = right / left; Y = nearer/farther, Z = down / up. Also, this coordinate
% system has its origin on the center of rotation of the eye.

% To pose the Blender model eye, we specify gaze targets in the Blender X,
% Y, Z sceneWorld coordinate space. We define azimuthal and elevation
% positions of the eye with respect to extrinsic (fixed) head-centered
% coordinates.

% Define the location of the center of the pupil in Blender
% sceneCoordinates when the eye is in primary gaze

for ii=1:nFrames
    
    % Get the rotation values for this frame
    eyeAzimuth = pupilCenterAzimuths(ii);
    eyeElevation = pupilCenterElevations(ii);
    eyeTorsion = 0;
    
    % Define a rotation matrix
    R3 = [cosd(eyeAzimuth) -sind(eyeAzimuth) 0; sind(eyeAzimuth) cosd(eyeAzimuth) 0; 0 0 1];
    R2 = [cosd(eyeElevation) 0 sind(eyeElevation); 0 1 0; -sind(eyeElevation) 0 cosd(eyeElevation)];
    R1 = [1 0 0; 0 cosd(eyeTorsion) -sind(eyeTorsion); 0 sind(eyeTorsion) cosd(eyeTorsion)];
    
    % This order (1-2-3) corresponds to a head-fixed, extrinsic, rotation
    % matrix.
    eyeRotation = R1*R2*R3;
    
    % Define the location of the eye center of rotation in the head-centered
    % coordinate frame
    pupilCenter = [-sceneGeometry.eyeRadius 0 0];
    
    % Apply the head-fixed rotation to center of the pupil
    pupilCenterCoords = (eyeRotation*pupilCenter')';
    
    % Arrange the dimensions and signs to corresponds to the Blender model
    gazeTargetPositionX(ii) = pupilCenterCoords(2)*-1;
    gazeTargetPositionY(ii) = pupilCenterCoords(3);
    gazeTargetPositionZ(ii) = pupilCenterCoords(1);
    pupilRadii(ii) = 2; % fix this at 2 mm
    eyeClosedness(ii) = 0; % fix this at fully open
end


%% generate eye movie
generateEyeMovie(codeDirectory,exportsDirectory,gazeTargetPositionX, gazeTargetPositionY, gazeTargetPositionZ, pupilRadii, eyeClosedness, cameraDepthPositionMM)
 
% rename the movie file
movefile(fullfile(exportsDirectory,'pupil_movie.avi'),fullfile(exportsDirectory,[pathParams.runName '_gray.avi']));

% move the rendered frames into a sub-directory
renderedFramesDirectory = '~/Desktop/Blender_Simulate01/renderedFrames';
if exist(renderedFramesDirectory,'dir')==0
    mkdir(renderedFramesDirectory);
end
system(['mv ' exportsDirectory '/reconstructedFrame*png ' renderedFramesDirectory '/']);


%% Run the processing pipeline
runVideoPipeline( pathParams, ...
    'verbosity', 'full', 'useParallel',false, 'catchErrors', false,...
    'pupilFrameMask', [60 60], 'maskBox', [1 1], 'pupilGammaCorrection',0.5, 'pupilRange', [20 80], ...
    'overwriteControlFile',true, 'glintPatchRadius', 10,'cutErrorThreshold',0.5, ...
    'ellipseTransparentLB',[0,0, 20, 0, 0],...
    'ellipseTransparentUB',[sceneResolutionPx(1),sceneResolutionPx(2), 20000, 1.0, pi],...
    'candidateThetas',0:pi/16:2*pi,...
    'badFrameErrorThreshold', 8, ...
    'makeFitVideoByNumber',6,...
    'skipStageByNumber',1,'lastStage','fitPupilPerimeter');

% Load the result files into memory
load(fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_pupil.mat']));



%% Test the forward projection model

% Obtain the ellipses found on the image plane
imagePlaneEllipsesObserved = pupilData.ellipseParamsUnconstrained_mean;
imagePlaneEllipseCentersObserved = imagePlaneEllipsesObserved(:,1:2);

% Open a figure and plot the observed ellipse centers
figure
plot(imagePlaneEllipseCentersObserved(:,1),imagePlaneEllipseCentersObserved(:,2),'o','MarkerEdgeColor',[.7 .7 .7],...
    'MarkerFaceColor',[.7 .7 .7])
axis equal
hold on

% Can we invert the Blender model? We test if the rendered ellipses can be
% modeled with the sceneGeometry within a 1% error on matching ellipse
% theta, eccentrivity, and error.
constraintTolerance = 0.01;
for ii=1:nFrames
    [eyeParams, projectedEllipsesOnImagePlane(ii,:), centerError(ii), shapeError(ii), areaError(ii)] = ...
        pupilProjection_inv(imagePlaneEllipsesObserved(ii,:), sceneGeometry, 'constraintTolerance', constraintTolerance);
    eyeParamError(ii,:)=eyeParams-[pupilCenterAzimuths(ii), pupilCenterElevations(ii), pupilRadii(ii)];
    ellipseParamError(ii,:)=projectedEllipsesOnImagePlane(ii,:)-imagePlaneEllipsesObserved(ii,:);
end
fprintf('Attempt to model observed image plane ellipses using design-specified sceneGeometry:\n')
fprintf('\t %d of %d ellipses were fit with < %f error in shape \n',sum(shapeError<constraintTolerance),nFrames,constraintTolerance);
fprintf('\t Maximum eye azimuth error: %f \n',max(eyeParamError(:,1)));
fprintf('\t Maximum eye elevation error: %f \n',max(eyeParamError(:,2)));
fprintf('\t Maximum pupil radius error: %f \n',max(eyeParamError(:,3)));
fprintf('\t Maximum ellipse center distance error: %f \n',max(centerError));
fprintf('\t Maximum ellipse area error: %f \n',max(areaError));

% Update the plot with these centers
plot(projectedEllipsesOnImagePlane(:,1),projectedEllipsesOnImagePlane(:,2),'.r')

% Search across sceneGeometry parameters to improve the fit to Blender
[adjustedSceneGeometry, distanceError] = findExtrinsicTranslationVector(imagePlaneEllipsesObserved, sceneGeometry);

% Can we invert the Blender model? Test again, now with the adjusted scene
% geometry
constraintTolerance = 0.01;
for ii=1:nFrames
    [eyeParams, projectedEllipsesOnImagePlane(ii,:), centerError(ii), shapeError(ii), areaError(ii)] = ...
        pupilProjection_inv(imagePlaneEllipsesObserved(ii,:), adjustedSceneGeometry, 'constraintTolerance', constraintTolerance);
    eyeParamError(ii,:)=eyeParams-[pupilCenterAzimuths(ii), pupilCenterElevations(ii), pupilRadii(ii)];
    ellipseParamError(ii,:)=projectedEllipsesOnImagePlane(ii,:)-imagePlaneEllipsesObserved(ii,:);
end
fprintf('Attempt to model observed image plane ellipses using adjusted sceneGeometry:\n')
fprintf('\t %d of %d ellipses were fit with < %f error in shape \n',sum(shapeError<constraintTolerance),nFrames,constraintTolerance);
fprintf('\t Maximum eye azimuth error: %f \n',max(eyeParamError(:,1)));
fprintf('\t Maximum eye elevation error: %f \n',max(eyeParamError(:,2)));
fprintf('\t Maximum pupil radius error: %f \n',max(eyeParamError(:,3)));
fprintf('\t Maximum ellipse center distance error: %f \n',max(centerError));
fprintf('\t Maximum ellipse area error: %f \n',max(areaError));


% Update the plot with these centers
plot(projectedEllipsesOnImagePlane(:,1),projectedEllipsesOnImagePlane(:,2),'.g')
legend({'observed ellipse centers','modeled with design sceneGeometry','modeled with adjusted sceneGeometry'})
hold off



