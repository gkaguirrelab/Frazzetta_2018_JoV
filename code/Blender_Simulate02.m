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
intrinsicCameraMatrix = ...
    [focalLengthPX 0 opticalCenterPx(1);...
    0 focalLengthPX opticalCenterPx(2);...
    0 0 1];


%% Hard coded aspect of the scene geometry
% The eyeRadius is set in the routine eyeModel/__init__.py
eyeRadius = 12;   % in mm


%% Variable aspects of scene geometry
% Variable aspect of the scene geometry. The Blender coordinate system
% defines depth zero as being at the center of rotation of the eye.
% Therefore, the depth position of the camera is relative to the center of
% rotation of the eye. The distance of the camera from the imaged scene,
% therefore, is equal to [cameraDepthPosition - eyeRadius]. We would like
% it to be the case that the extrinsicTranslationVector is estimated to
% have this scene distance value (50-12) for the Z dimension.
cameraDepthPositionMM = 50; % in mm

extrinsicTranslationVector = [0; 0; cameraDepthPositionMM - eyeRadius];
extrinsicRotationMatrix = [1 0 0; 0 -1 0; 0 0 -1];

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
    pupilCenter = [-eyeRadius 0 0];
    
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
% We force the sceneGeometry model to match the design parameters used by
% Blender to render the scene.
runVideoPipeline( pathParams, ...
    'verbosity', 'full', 'useParallel',false, 'catchErrors', false,...
    'pupilFrameMask', [60 60], 'maskBox', [1 1], 'pupilGammaCorrection',0.5, 'pupilRange', [20 80], ...
    'overwriteControlFile',true, 'glintPatchRadius', 10,'cutErrorThreshold',0.5, ...
    'ellipseTransparentLB',[0,0, 20, 0, 0],...
    'ellipseTransparentUB',[sceneResolutionPx(1),sceneResolutionPx(2), 20000, 1.0, pi],...
    'candidateThetas',0:pi/16:2*pi,...
    'badFrameErrorThreshold', 8, ...
    'makeFitVideoByNumber',6,...
    'skipStageByNumber',1,...
    'intrinsicCameraMatrix',intrinsicCameraMatrix,...
    'extrinsicRotationMatrix',extrinsicRotationMatrix,...
    'extrinsicTranslationVector',extrinsicTranslationVector,...
    'extrinsicTranslationVectorLB',extrinsicTranslationVector,...
    'extrinsicTranslationVectorUB',extrinsicTranslationVector,...
    'eyeRadius',eyeRadius,...
    'eyeRadiusLB',eyeRadius,...
    'eyeRadiusUB',eyeRadius...
    );
