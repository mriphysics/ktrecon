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
fcmrNo = 168;
seriesNos = [22:28,30:34,36:40,42:46,47:50];
rawDataDirPath = '/home/jva13/mnt/pnraw01-ingenia/2017_05_18/OF_295686';

% Dependencies
cd ~/ktrecon  
addpath( '~/ktrecon/mrecon' )                         % required for id_pnraw_data
addpath( genpath( '~/reconframe/MRecon-3.0.535' ) ),  % required for ReconFrame

% Processing
outputDirPath  = fullfile( '/scratch/jva13', sprintf( 'fcmr%03i', fcmrNo ) );
for seriesNo = seriesNos,
  idStr = sprintf( 'fcmr%03is%02i', fcmrNo, seriesNo );
  fprintf( '\n============ %s ============\n\n', idStr )
  [ rawDataFilePath, coilSurveyFilePath, senseRefFilePath ] = id_pnraw_data( rawDataDirPath, seriesNo );
  mrecon_kt( rawDataFilePath, 'senseref', senseRefFilePath, 'coilsurvey', coilSurveyFilePath, 'outputdir', outputDirPath, 'outputname', idStr ) 
end
```