% Blender_Simulate01
% Create a video of a simulated eye and recover the input parameters
%
% Description:
%   This simulation makes use of the Blender eye model for the generation
%   of ground-truth data introduced by Lech ?wirski and Neil A. Dodgson:
%       http://www.cl.cam.ac.uk/research/rainbow/projects/eyerender/
%
%   We create a video with the eye posed at different azimuths and
%   elevations. We then determine how well we can recover these parameters
%   of the eye from an analysis of the resulting video. For this
%   simulation, we aim to place the camera (and the IR light source)
%   exactly on the line that connects the center of rotation of the eye
%   with the center of the pupil when the eye is in primary position
%   (Azimuth = 0; Elevation = 0).
%
%   If fully succesful, we would be able to specify the sceneGeometry
%   parameters given only the settings used for the Blender movie
%   synthesis. Using these parameters, we would be able to perfectly fit
%   the image of the pupil in the video and recover the "ground truth"
%   rotation of the eye and pupil radius. In practice, we face the
%   following limitations:
%     - Our knowledge of the parameters of the Blender movie generation is
%       incomplete. As evidence of this, the image that results for an
%       ostensible azimuth = 0; elevation = 0 is slightly off from primary
%       gaze. We find that applying small adjustments to the X and Y gaze
%       targets can be used to center the pupil, but we are uncertain as to
%       why this is necessary.
%     - Our forward model of the pupil projection is imperfect. Two 
%       examples:
%        1) The Blender eye incorporates a model of the cornea that
%           has some refractive index. This causes a change in the appearance
%           of the "entrance pupil" of the eye in the image plane. Our model
%           does not include this distortion of the size and shape of the pupil
%           image, and thus cannot perfectly model the Blender data at far
%           rotations of the eye.
%        2) The appearance of the pupil on the image plane deviates from an
%           ellipse and becomes more egg-shaped for far rotations. Our
%           model will imerfectly capture this shape and recover the
%           eyeParams.
%   

%% Housekeeping
close all
clear all


%% Set paths and file names
% This path should be defined in the eyeModelSupport local hook
codeDirectory = '~/Documents/MATLAB/toolboxes/eyemodelSupport/code';

% Where to save the results of our simulation
exportsDirectory = '~/Desktop/Simulate01_AzimuthElevationSweeps';
pathParams.dataOutputDirFull = fullfile(exportsDirectory);
pathParams.dataSourceDirFull = fullfile(exportsDirectory);
pathParams.runName = 'Simulate01_AzimuthElevationSweeps';
videoName = fullfile(exportsDirectory, [pathParams.runName '_gray.avi']);

% check or make a directory for output
if exist(exportsDirectory,'dir')==0
    mkdir(exportsDirectory);
end


%% Hard coded properties of the camera
% The focal length is set in the routine eyeModel/__init__.py
% All intrinsic camera values are in units of pixels
sceneResolutionPx = [640 480];
opticalCenterPx = (sceneResolutionPx./2);
focalLengthPX = (sceneResolutionPx(1)/2.0) / tan(45*pi/180 / 2);
radialDistortionVector = [0 0];
intrinsicCameraMatrix = ...
    [focalLengthPX 0 opticalCenterPx(1);...
    0 focalLengthPX opticalCenterPx(2);...
    0 0 1];


%% Hard coded aspect of the scene geometry
% The eyeRadius is set in the routine eyeModel/__init__.py
% While it is specified as 12.0 in the Blender scene, we find that a value
% of 12.15 best fits the data. This may be because there is magnification
% of the pupil image by the cornea that we are not modeling.
eyeRadius = 12.15;   % in mm
extrinsicRotationMatrix = [1 0 0; 0 -1 0; 0 0 -1];


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
    % We find that a shift in the X and Y position of the gaze target is
    % needed to place the center of the pupil in line with the center of
    % rotation of the eye when azimuth, elevation = 0. We have been unable
    % to determine what aspect of the Blender model produces this
    % imperfection.
    gazeTargetPositionX(ii) = pupilCenterCoords(2)*-1+0.05;
    gazeTargetPositionY(ii) = pupilCenterCoords(3)-0.15;
    gazeTargetPositionZ(ii) = pupilCenterCoords(1);
    pupilRadii(ii) = 2.5; % fix this at 2.5 mm
    eyeClosedness(ii) = 0; % fix this at fully open
end


%% generate eye movie
%generateEyeMovie(codeDirectory,exportsDirectory,gazeTargetPositionX, gazeTargetPositionY, gazeTargetPositionZ, pupilRadii, eyeClosedness, cameraDepthPositionMM)

% rename the movie file
%movefile(fullfile(exportsDirectory,'pupil_movie.avi'),fullfile(exportsDirectory,[pathParams.runName '_gray.avi']));

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
    'radialDistortionVector',radialDistortionVector,...
    'extrinsicRotationMatrix',extrinsicRotationMatrix,...
    'extrinsicTranslationVector',extrinsicTranslationVector,...
    'extrinsicTranslationVectorLB',extrinsicTranslationVector,...
    'extrinsicTranslationVectorUB',extrinsicTranslationVector,...
    'constraintTolerance',0.03,...
    'eyeRadius',eyeRadius,...
    'eyeRadiusLB',eyeRadius,...
    'eyeRadiusUB',eyeRadius ...
    );


%% Generate result plots

% Load the pupilData into memory
load(fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_pupil.mat']));

figure
subplot(1,3,1)
plot(pupilCenterAzimuths,pupilData.radiusSmoothed.eyeParams.values(:,1), '.r')
rl = refline(1,0);
rl.Color = 'k';
xlabel('Ground Truth Pupil Azimuth in degrees')
ylabel('Reconstructed Pupil Azimuth in degrees')
ylim([-40 40]);
xlim([-40 40]);
axis square

subplot(1,3,2)
plot(pupilCenterElevations,pupilData.radiusSmoothed.eyeParams.values(:,2), '.r')
rl = refline(1,0);
rl.Color = 'k';
xlabel('Ground Truth Pupil Elevation in degrees')
ylabel('Reconstructed Pupil Elevation in degrees')
ylim([-30 30]);
xlim([-30 30]);
axis square

subplot(1,3,3)
plot(pupilData.radiusSmoothed.eyeParams.values(:,3), '.r');
hold on
plot(pupilRadii, '-k');
ylabel('pupil radius reconstructed [mm]');
ylim([2 3]);
axis square
hold off

