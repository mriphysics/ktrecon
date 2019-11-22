function recon_exam( fcmrNo, seriesNos, rawDataDirPath, patchVersion, isGeoCorrn, maskDirPath, cusTrnDirPath, isHarmonicFilter, isSelfCaliPreProc, outputDirPath )
%RECON_EXAM  reconstruct multiple series in one MRI exam

% jfpva (joshua.vanamerom@kcl.ac.uk)
% tar (t.roberts@kcl.ac.uk)


%% Dependencies

origPath  = path;
resetPath = onCleanup( @() path(origPath) );

addpath( fullfile( fileparts( which( 'mrecon_kt' ) ), 'mrecon' ) )  % required for id_pnraw_data
addpath( fullfile( fileparts( which( 'mrecon_kt' ) ), 'lib', 'nifti' ) )  % required to load mask file


%% Geometry Correction

suffixStr = '';
if ( isGeoCorrn )
    suffixStr = strcat( suffixStr, '_geocorrn' );
    reconOpts = { 'GeometryCorrection', 'Yes' };
end

isMask = false;
if exist('maskDirPath','var')
    if exist(maskDirPath,'dir')
        suffixStr = strcat( suffixStr, '_mask' );
        isMask    = true;
    end
end


%% Custom Training Data
isCusTrnData = false;
if exist('cusTrnDirPath','var')
    if exist(cusTrnDirPath,'dir')        
        suffixStr    = strcat( suffixStr, '_cusTrnData' );
        isCusTrnData = true;
    end
end


%% Make x-f Harmonic + Low Pass Filter
makeHarmonicFilter = false;
if ( isHarmonicFilter )
    suffixStr          = strcat( suffixStr, '_xfHrmFlt' );
    makeHarmonicFilter = true;
end


%% Self-calibrated Sliding Window k-t Preprocessing
if ~isSelfCaliPreProc
    isSelfCaliPreProc = false;
end


%% Output Directory
if nargin < 10
    outputDirPath  = fullfile( '/scratch/tr17/ktrecon', sprintf( 'fcmr%03i%s', fcmrNo, suffixStr ) );
end
    
if ~exist(outputDirPath,'dir')
    mkdir(outputDirPath)
end


%% Recon

diary( fullfile( outputDirPath, sprintf( 'log_mrecon_kt_%s.txt', datestr(now,'yyyymmdd_HHMMss') ) ) )

for seriesNo = seriesNos

    idStr = sprintf( 's%02i', seriesNo );
    fprintf( '\n============ %s ============\n\n', idStr )
    [ rawDataFilePath, coilSurveyFilePath, senseRefFilePath ] = id_pnraw_data( rawDataDirPath, seriesNo );
    
    if isHarmonicFilter
        fprintf( 'Performing reconstruction using x-f harmonic filter. \n' );
    end
    
    
    % Adaptive Regularization
    if ( isMask )
        
        % Get Mask
        maskFilePath  = fullfile( maskDirPath, ['s' num2str(seriesNo) '_mask_heart.nii.gz'] );
        fprintf( 'mask file:       %s\n', maskFilePath );
        niiMask       = load_untouch_nii( maskFilePath );
        mask          = logical(niiMask.img);
        
        % Custom Training Data
        if ~( isCusTrnData )
            
            mrecon_kt( rawDataFilePath, 'senseref', senseRefFilePath, 'coilsurvey', coilSurveyFilePath, 'outputdir', outputDirPath, 'outputname', idStr, 'patchversion', patchVersion, 'reconoptionpairs', reconOpts, 'mask', mask, 'isSelfCaliPreProc', isSelfCaliPreProc, 'makeHarmonicFilter', makeHarmonicFilter )
        
        elseif ( isCusTrnData )
            
            cusTrnMatFilePath = fullfile( cusTrnDirPath, strcat( ['s' num2str(seriesNo), '_sc_slw_recon.mat'] ) );  
            fprintf( 'custom training data file:       %s\n', cusTrnMatFilePath );
            mrecon_kt( rawDataFilePath, 'senseref', senseRefFilePath, 'coilsurvey', coilSurveyFilePath, 'outputdir', outputDirPath, 'outputname', idStr, 'patchversion', patchVersion, 'reconoptionpairs', reconOpts, 'mask', mask, 'isSelfCaliPreProc', isSelfCaliPreProc, 'cusTrnDirPath', cusTrnDirPath, 'makeHarmonicFilter', makeHarmonicFilter )
        
        end
        
        
    % Uniform Regularization
    else
        
        mrecon_kt( rawDataFilePath, 'senseref', senseRefFilePath, 'coilsurvey', coilSurveyFilePath, 'outputdir', outputDirPath, 'outputname', idStr, 'patchversion', patchVersion, 'reconoptionpairs', reconOpts )
    
    end
    
end


diary off


end  % recon_exam(...)
