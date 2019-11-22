function [  ] = mrecon_arrangeKT( rawDataFilePath, senseRefFilePath, coilSurveyFilePath, outputDirPath, outFilePrefix )
%MRECON_ARRANGEKT   MRecon reconstruction of Philips KT data
%
%   [ ] = MRECON_ARRANGEKT( rawDataFilePath, senseRefFilePath, coilSurveyFilePath )
%
%   output:
%       matFile containing MRecon processed k-space
%
%   input: 
%       rawDataFilePath             path to raw data file
%       senseRefFilePath            path to sense ref file
%       coilSurveyFilePath          path to coil survey file
%       outputDirPath               path to output directory
%       outFilePrefix               output file name prefix
%
%   see also: mrecon_sc_slw, arrangeKT

% tar (t.roberts@kcl.ac.uk)
% wrapper for code by Lucilio Cordero-Grande


%% File rename to match Lucilio's convention
datFile = rawDataFilePath;
refFile = senseRefFilePath;
surFile = coilSurveyFilePath;

outFileSuffix = '_kspace_sc_slw.mat';
kspaceMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, outFileSuffix ) );

cd(outputDirPath);

fprintf( '\n============ %s ============\n\n', outFilePrefix );


%% Setup
gpu=gpuDeviceCount;
if gpu;gpuF=2;else gpuF=0;end
if gpu;dev=gpuDevice;end

comp = 1; %To compress the data when storing

dyn = []; % To select a series of dynamics
sl  = []; % To select a series of slices


%% ARRANGE DATA
% arrangeKT uses mRecon

fprintf( 'Running ArrangeKT on %s ... \n\n', outFilePrefix );

inst=1;
if inst == 1
    if gpu;wait(dev);end;tsta=tic;
    [y,S,A,~,R,fps,voxsiz]=arrangeKT(datFile,refFile,surFile,dyn,comp);%Calibration data seems acquired in a different time, so we don't use it although there's some support inside, last parameter serves to set a series of dynamics                   
    if gpu;wait(dev);end;tend=toc(tsta);
    fprintf('Time reading RF: %.2fs\n',tend);
    if gpu;wait(dev);end;tsta=tic; 
    S=gather(S);y=gather(y);
    
    fprintf( 'Saving %s ... \n', strcat( outFilePrefix, outFileSuffix ) );
    save( kspaceMatFilePath, 'S','y','A','R','fps','voxsiz','-v7.3');
    fprintf( 'kspace saved: %s \n', kspaceMatFilePath );
    
    if gpu;wait(dev);end;tend=toc(tsta);
    fprintf('Time writing MAT: %.2fs\n',tend);
else
    if gpu;wait(dev);end;tsta=tic;        
    
    load( kspaceMatFilePath );
    
    if gpu;wait(dev);end;tend=toc(tsta);
    fprintf('Time reading MAT: %.2fs\n',tend);
end

fprintf( '\nArrangeKT complete. \n\n' );


end  % mrecon_arrangeKT(...)

