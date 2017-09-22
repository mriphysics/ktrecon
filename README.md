# ktrecon

Reconstruct dynamic k-t under-sampled MR data using ReconFrame.


## directories

* lib - external libraries
* mrecon - functions related to ReconFrame
* recon - k-t SENSE and sliding window reconstruction
* util - utility functions


## example

```matlab
% User-Specified
fcmrNo = 191;
seriesNos = [16:20];
rawDataDirPath = '/home/jva13/mnt/pnraw01-ingenia/2017_08_14/MA_362635/';
patchVersion = 'PIH1';

% Dependencies
cd ~/ktrecon  
addpath( '~/ktrecon/mrecon' )                         % required for id_pnraw_data
addpath( genpath( '~/reconframe/MRecon-3.0.535' ) ),  % required for ReconFrame

% Processing
outputDirPath  = fullfile( '/scratch/jva13', sprintf( 'fcmr%03i', fcmrNo ) );
if ~exist(outputDirPath,'dir')
    mkdir(outputDirPath)
end
diary( fullfile( outputDirPath, sprintf( 'log_mrecon_kt_%s.txt', datestr(now,'yyyymmdd_HHMMss') ) ) )
for seriesNo = seriesNos,
    idStr = sprintf( 'fcmr%03is%02i', fcmrNo, seriesNo );
    fprintf( '\n============ %s ============\n\n', idStr )
    [ rawDataFilePath, coilSurveyFilePath, senseRefFilePath ] = id_pnraw_data( rawDataDirPath, seriesNo );
    mrecon_kt( rawDataFilePath, 'senseref', senseRefFilePath, 'coilsurvey', coilSurveyFilePath, 'outputdir', outputDirPath, 'outputname', idStr, 'patchversion', patchVersion ) 
end
diary off
```
