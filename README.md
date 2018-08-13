# ktrecon

Reconstruct dynamic k-t under-sampled MR data using [ReconFrame](https://www.gyrotools.com/gt/index.php/products/reconframe).


## directories

* lib - external libraries
* mrecon - functions related to ReconFrame
* recon - k-t SENSE and sliding window reconstruction
* util - utility functions


## example

```matlab
% User-Specified
seriesNos = [16:20];
rawDataDirPath = '/path/to/raw/data/';
patchVersion = 'PIH1';

% Dependencies
cd ~/ktrecon  
addpath( '~/ktrecon/mrecon' )                         % required for id_pnraw_data
addpath( genpath( '~/reconframe/MRecon-3.0.535' ) ),  % required for ReconFrame

% Processing
outputDirPath  = '/path/for/output';
if ~exist(outputDirPath,'dir')
    mkdir(outputDirPath)
end
diary( fullfile( outputDirPath, sprintf( 'log_mrecon_kt_%s.txt', datestr(now,'yyyymmdd_HHMMss') ) ) )
for seriesNo = seriesNos,
    idStr = sprintf( 's%02i', seriesNo );
    fprintf( '\n============ %s ============\n\n', idStr )
    [ rawDataFilePath, coilSurveyFilePath, senseRefFilePath ] = id_pnraw_data( rawDataDirPath, seriesNo );
    mrecon_kt( rawDataFilePath, 'senseref', senseRefFilePath, 'coilsurvey', coilSurveyFilePath, 'outputdir', outputDirPath, 'outputname', idStr, 'patchversion', patchVersion ) 
end
diary off
```
