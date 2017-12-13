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


% features of the camera for the demo (standard blender camera)
sensorSize = [36 24]; %mm
aperture = 2.0;
sceneResolution = [640 480];
focalLength = (sceneResolution(1)/aperture) / tand(45);

pixelSizeMM = sensorSize./sceneResolution;


% eye radius and camera distance (hardcoded in the python function)
eyeRadiusMM = 14;
cameraDistanceMM = 60;

eyeRadiusPX = eyeRadiusMM/pixelSizeMM(1);
sceneDistancePX = cameraDistanceMM/pixelSizeMM(1);

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

rotationArm = abs(eyeRadiusMM);
% pupil position based on the angles (eye target)
pupilXpos = (rotationArm).*sind(allPupilAzi).*cosd(allPupilEle);
pupilYpos = (rotationArm).*sind(allPupilEle);
pupilZpos =-(rotationArm).*cosd(allPupilAzi).*cosd(allPupilEle);

pupilRadiusMM = 2*ones(1,length(pupilXpos));
pupilRadiusPX = pupilRadiusMM/pixelSizeMM(1);
eyeClosedness = 0*ones(1,length(pupilXpos));


%% generate eye movie
generateEyeMovie(codeDirectory, exportsDirectory, pupilXpos, pupilYpos, pupilZpos, pupilRadiusMM, eyeClosedness, cameraDistanceMM)

% rename file
movefile(fullfile(exportsDirectory,'pupil_movie.avi'),fullfile(exportsDirectory,[pathParams.runName '_gray.avi']));

%% Run the processing pipeline
runVideoPipeline( pathParams, ...
    'verbosity', 'full', 'useParallel',false, 'catchErrors', false,...
    'maskBox', [0.9 0.9], ...
    'overwriteControlFile',true, 'glintPatchRadius', 10,  ...
    'projectionModel', 'pseudoPerspective', ...
    'ellipseTransparentLB',[0,0, 20, 0, 0],...
    'ellipseTransparentUB',[sceneResolution(1),sceneResolution(2), 20000, 1.0, pi],...
    'candidateThetas',0:pi/16:2*pi,...
    'badFrameErrorThreshold', 8, ...
    'eyeRadius',eyeRadiusPX, 'cameraDistanceInPixels',sceneDistancePX, ...
    'sceneGeometryLB',[0,  0, sceneDistancePX+eyeRadiusPX, 100],'sceneGeometryUB',[sceneResolution(1),  sceneResolution(2), sceneDistancePX+eyeRadiusPX, 500],...
    'skipStageByNumber',1);

% Load the result files into memory
load(fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_pupil.mat']));
load(fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_sceneGeometry.mat']));

%% Report how closely we have reconstructed the actual scene geometry
fprintf('Veridical scene geometry - eye center: [%0.1f, %0.1f, %0.1f], eye radius: %0.1f \n',sceneResolution(1)/2,  sceneResolution(2)/2, sceneDistancePX+eyeRadiusPX, eyeRadiusPX);
fprintf('Estimated scene geometry - eye center: [%0.1f, %0.1f, %0.1f], eye radius: %0.1f \n',sceneGeometry.eyeCenter.X, sceneGeometry.eyeCenter.Y, sceneGeometry.eyeCenter.Z, sceneGeometry.eyeRadius);

%% use the inverse projection function and the scene geometry to recover the angles

reconstructedPupilAzi =[];
reconstructedPupilEle = [];

for ii=1:size(pupilData.ellipseParamsAreaSmoothed_mean,1)
    transparentEllipse = pupilData.ellipseParamsAreaSmoothed_mean(ii,:);
    [reconstructedPupilAzi(ii), reconstructedPupilEle(ii), reconstructedPupilArea(ii)] = pupilProjection_inv(transparentEllipse, [sceneGeometry.eyeCenter.X sceneGeometry.eyeCenter.Y sceneGeometry.eyeCenter.Z], sceneGeometry.eyeRadius, sceneGeometry.meta.projectionModel);
end

% plot real Azi,  Ele, and pupil area vs reconstructed ones
figure
subplot(1,2,1)
plot(allPupilAzi,reconstructedPupilAzi, '.r')
rl = refline(1,0);
rl.Color = 'k';
xlabel('Ground Truth Pupil Azimuth in degrees')
ylabel('Reconstructed Pupil Azimuth in degrees')
ylim([min([min(allPupilAzi) min(reconstructedPupilAzi)]) max([max(allPupilAzi) max(reconstructedPupilAzi)])]);
xlim([min([min(allPupilAzi) min(reconstructedPupilAzi)]) max([max(allPupilAzi) max(reconstructedPupilAzi)])]);
axis square

subplot(1,2,2)
plot(-allPupilEle,reconstructedPupilEle, '.r')
rl = refline(1,0);
rl.Color = 'k';
xlabel('Ground Truth Pupil Elevation in degrees')
ylabel('Reconstructed Pupil Elevation in degrees')
ylim([min([min(allPupilAzi) min(reconstructedPupilAzi)]) max([max(allPupilAzi) max(reconstructedPupilAzi)])]);
xlim([min([min(allPupilAzi) min(reconstructedPupilAzi)]) max([max(allPupilAzi) max(reconstructedPupilAzi)])]);
axis square

%% plot some fits

temporalSupport = 0:1/60.:(size(pupilData.ellipseParamsSceneConstrained_mean,1)-1)/60; % seconds
temporalSupport = temporalSupport / 60; % minutes

% Make a plot of pupil area, both on the image plane and on the eye
figure
subplot(2,1,1)
plot(temporalSupport,pupilData.ellipseParamsUnconstrained_mean(:,3),'-k','LineWidth',2);
hold on
plot(temporalSupport,pupilData.ellipseParamsSceneConstrained_mean(:,3),'-b');
plot(temporalSupport,pupilData.ellipseParamsSceneConstrained_mean(:,3)-pupilData.ellipseParamsSceneConstrained_splitsSD(:,3),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.ellipseParamsSceneConstrained_mean(:,3)+pupilData.ellipseParamsSceneConstrained_splitsSD(:,3),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.ellipseParamsAreaSmoothed_mean(:,3),'-r','LineWidth',2)
xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('pupil area [pixels in plane]');
hold off

subplot(2,1,2)
plot(temporalSupport,pupilData.meta.smoothPupilArea.pupilAreaMean,'-k');
hold on
plot(temporalSupport,pupilData.meta.smoothPupilArea.pupilAreaMean-pupilData.meta.smoothPupilArea.pupilAreaSD,'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.meta.smoothPupilArea.pupilAreaMean+pupilData.meta.smoothPupilArea.pupilAreaSD,'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.meta.smoothPupilArea.pupilAreaMean,'-r','LineWidth',2)
xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('pupil area [pixels on eye]');
hold off
