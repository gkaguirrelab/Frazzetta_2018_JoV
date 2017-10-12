% DEMO

%% set paths and make directories
% create test sandbox on desktop
transparentDir = '/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/Frazzetta_201x_transparentTrack/ExperimentalData/';
if ~exist(transparentDir,'dir')
    mkdir(transparentDir)
end

%% hard coded parameters
nFrames = 1000; % number of frames to process (set to Inf to do all)
verbosity = 'full'; % Set to none to make the demo silent
TbTbProjectName = 'eyeTrackTOMEAnalysis';


%% TbTb configuration
% We will suppress the verbose output, but detect if there are deploy
% errors and if so stop execution
tbConfigResult=tbUseProject(TbTbProjectName,'reset','full','verbose',false);
if sum(cellfun(@sum,extractfield(tbConfigResult, 'isOk')))~=length(tbConfigResult)
    error('There was a tb deploy error. Check the contents of tbConfigResult');
end
tbSnapshot=tbDeploymentSnapshot(tbConfigResult,'verbose',false);
clear tbConfigResult

% identify the base for the project code directory
%  This would normally be used as the location to save the controlFiles
codeBaseDir = tbLocateProject(TbTbProjectName,'verbose',false);




%% Define some file names

% TOME_3018
grayVideoName ='/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/TOME_processing/session1_restAndStructure/TOME_3018/040717/EyeTracking/rfMRI_REST_AP_run01_gray.avi';
glintFileName = fullfile(transparentDir, '3018rfMRI_REST_AP_run01_glint.mat');
perimeterFileName = fullfile(transparentDir, '3018rfMRI_REST_AP_run01_TRANSPARENTperimeter.mat');
starburstPerimeterFileName = fullfile(transparentDir, '3018rfMRI_REST_AP_run01_STARBURSTperimeter.mat');
controlFileName = fullfile(transparentDir, '3018rfMRI_REST_AP_run01_TRANSPARENTcontrolFile.csv');
starburstControlFileName = fullfile(transparentDir, '3018rfMRI_REST_AP_run01_STARBURSTcontrolFile.csv');
correctedPerimeterFileName = fullfile(transparentDir, '3018rfMRI_REST_AP_run01_correctedTRANSPARENTPerimeter.mat');
starburstCorrectedPerimeterFileName = fullfile(transparentDir, '3018rfMRI_REST_AP_run01_correctedSTARBURSTPerimeter.mat');
pupilFileName = fullfile(transparentDir, '3018rfMRI_REST_AP_run01_pupil.mat');
starburtstPupilFileName = fullfile(transparentDir, '3018rfMRI_REST_AP_run01_STARBURSTpupil.mat');
finalFitVideoName = fullfile(transparentDir, '3018rfMRI_REST_AP_run01_finalFit.avi');
irisFileName = fullfile(transparentDir, '3018rfMRI_REST_AP_run01_iris.mat');
starburstFitVideoName = fullfile(transparentDir, '3018rfMRI_REST_AP_run01_STARBURSTfinalFit.avi');


% TOME_3020
grayVideoName ='/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/tfMRI_RETINO_PA_run01_gray.avi';
glintFileName = fullfile(transparentDir, '3020tfMRI_RETINO_PA_run01_glint.mat');
perimeterFileName = fullfile(transparentDir, '3020tfMRI_RETINO_PA_run01TRANSPARENTperimeter.mat');
starburstPerimeterFileName = fullfile(transparentDir, '3020tfMRI_RETINO_PA_run01_STARBURSTperimeter.mat');
controlFileName = fullfile(transparentDir, '3020tfMRI_RETINO_PA_run01_TRANSPARENTcontrolFile.csv');
starburstControlFileName = fullfile(transparentDir, '3020tfMRI_RETINO_PA_run01_STARBURSTcontrolFile.csv');
correctedPerimeterFileName = fullfile(transparentDir, '3020tfMRI_RETINO_PA_run01_correctedTRANSPARENTPerimeter.mat');
starburstCorrectedPerimeterFileName = fullfile(transparentDir, '3020tfMRI_RETINO_PA_run01_correctedSTARBURSTPerimeter.mat');
pupilFileName = fullfile(transparentDir, '3020tfMRI_RETINO_PA_run01_pupil.mat');
starburtstPupilFileName = fullfile(transparentDir, '3020tfMRI_RETINO_PA_run01_STARBURSTpupil.mat');
finalFitVideoName = fullfile(transparentDir, '3020tfMRI_RETINO_PA_run01_finalFit.avi');
irisFileName = fullfile(transparentDir, '3020tfMRI_RETINO_PA_run01_iris.mat');
starburstFitVideoName = fullfile(transparentDir, '3020tfMRI_RETINO_PA_run01_STARBURSTfinalFit.avi');

% find glint
findGlint(grayVideoName,glintFileName, 'nFrames', nFrames)

% find pupil perimeter with thresholding
findPupilPerimeter(grayVideoName, perimeterFileName, 'nFrames', nFrames)

% correct OUR perimeter
makeControlFile(controlFileName,perimeterFileName,glintFileName, 'nFrames', nFrames)
applyControlFile(perimeterFileName,controlFileName,correctedPerimeterFileName, 'nFrames', nFrames)

% correct STARBURST perimeter
makeControlFile(starburstControlFileName,starburstPerimeterFileName,glintFileName, 'nFrames', nFrames, 'verbosity', 'full','glintPatchRadius', 15,  'cutErrorThreshold', 100)
applyControlFile(starburstPerimeterFileName,starburstControlFileName,starburstCorrectedPerimeterFileName, 'nFrames', nFrames)
fitPupilPerimeter(starburstCorrectedPerimeterFileName, starburtstPupilFileName, ...
    'verbosity', 'full', 'cutErrorThreshold', 100,'useParallel',true, 'exponentialTauParams', [.25, .25, 10, 1, 1],'likelihoodErrorExponent',[1.25 1.25 2 2 2])

makeFitVideo(grayVideoName, starburstFitVideoName, 'glintFileName', glintFileName, 'perimeterFileName', starburstCorrectedPerimeterFileName, 'controlFileName', starburstControlFileName, 'pupilFileName', starburtstPupilFileName,  ...
    'nFrames',nFrames,'verbosity', 'full', 'useParallel',true);









% %% Create videos of intermediate stages
% exampleFitVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_Stage1.avi']);
% makeFitVideo(grayVideoName, exampleFitVideoName, 'glintFileName', glintFileName, ...
%     'nFrames',nFrames,'verbosity', verbosity,'tbSnapshot',tbSnapshot, 'useParallel',true);
% 
% exampleFitVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_Stage2.avi']);
% makeFitVideo(grayVideoName, exampleFitVideoName, 'glintFileName', glintFileName, 'perimeterFileName', perimeterFileName, ...
%     'nFrames',nFrames,'verbosity', verbosity,'tbSnapshot',tbSnapshot, 'useParallel',true);
% 
% exampleFitVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_Stage3.avi']);
% makeFitVideo(grayVideoName, exampleFitVideoName, 'glintFileName', glintFileName, 'perimeterFileName', correctedPerimeterFileName, 'controlFileName', controlFileName, ...
%     'nFrames',nFrames,'verbosity', verbosity,'tbSnapshot',tbSnapshot, 'useParallel',true);
% 
% exampleFitVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_Stage4.avi']);
% makeFitVideo(grayVideoName, exampleFitVideoName, 'glintFileName', glintFileName, 'perimeterFileName', correctedPerimeterFileName, 'controlFileName', controlFileName, 'pupilFileName', pupilFileName, ...
%     'whichFieldToPlot','pInitialFitTransparent',...
%     'nFrames',nFrames,'verbosity', verbosity,'tbSnapshot',tbSnapshot, 'useParallel', true);
% 
% %% Perform Bayesian smoothing and iris fitting
% runVideoPipeline( pathParams, ...
%     'whichFieldToPlot','pInitialFitTransparent',...
%     'nFrames',nFrames,'verbosity', verbosity,'tbSnapshot',tbSnapshot, 'useParallel',true, ...
%     'skipInitialPupilFit', true, ...
%     'skipStage', {'convertRawToGray', 'findGlint', 'findPupilPerimeter', 'makeControlFile', 'applyControlFile', 'makeFitVideo' }, ...
%     'catchErrors', false);
% 
% %% Create some more videos
% exampleFitVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_Stage5.avi']);
% makeFitVideo(grayVideoName, exampleFitVideoName, 'glintFileName', glintFileName, 'perimeterFileName', correctedPerimeterFileName, 'controlFileName', controlFileName, 'pupilFileName', pupilFileName, ...
%     'nFrames',nFrames,'verbosity', verbosity,'tbSnapshot',tbSnapshot, 'useParallel',true);
% 
% exampleFitVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_Stage6.avi']);
% makeFitVideo(grayVideoName, exampleFitVideoName, 'glintFileName', glintFileName, 'perimeterFileName', correctedPerimeterFileName, 'controlFileName', controlFileName, 'pupilFileName', pupilFileName, 'irisFileName', irisFileName, ...
%     'nFrames',nFrames,'verbosity', verbosity,'tbSnapshot',tbSnapshot, 'useParallel',true);
% 
% 
% %% Plot some fits
% pupilFileName = fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_pupil.mat']);
% dataLoad = load(pupilFileName);
% pupilData = dataLoad.pupilData;
% clear dataLoad
% 
% temporalSupport = 0:1/60.:(size(pupilData.pPosteriorMeanTransparent,1)-1)/60; % seconds
% temporalSupport = temporalSupport / 60; % minutes
% 
% figure
% plot(temporalSupport,pupilData.pInitialFitTransparent(:,3),'-.k');
% hold on
% plot(temporalSupport,pupilData.pPosteriorMeanTransparent(:,3),'-r','LineWidth',2)
% plot(temporalSupport,pupilData.pPosteriorMeanTransparent(:,3)-pupilData.pPosteriorSDTransparent(:,3),'-b')
% plot(temporalSupport,pupilData.pPosteriorMeanTransparent(:,3)+pupilData.pPosteriorSDTransparent(:,3),'-b')
% xlim([0 max(temporalSupport)]);
% xlabel('time [mins]');
% ylabel('area [pixels]');
% hold off
% 
% figure
% plot(temporalSupport,pupilData.pInitialFitTransparent(:,1),'-.k');
% hold on
% plot(temporalSupport,pupilData.pPosteriorMeanTransparent(:,1),'-r','LineWidth',2)
% plot(temporalSupport,pupilData.pPosteriorMeanTransparent(:,1)-pupilData.pPosteriorSDTransparent(:,1),'-b')
% plot(temporalSupport,pupilData.pPosteriorMeanTransparent(:,1)+pupilData.pPosteriorSDTransparent(:,1),'-b')
% xlim([0 max(temporalSupport)]);
% xlabel('time [mins]');
% ylabel('position [pixels]');
% hold off
