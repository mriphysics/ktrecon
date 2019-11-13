function recon_exam_hrh( fcmrNo, seriesNos, rawDataDirPath, patchVersion, isGeoCorrn, maskDirPath, outputDirPath )
%RECON_EXAM  reconstruct multiple series in one MRI exam

% jfpva (joshua.vanamerom@kcl.ac.uk)


%% Dependencies

origPath  = path;
resetPath = onCleanup( @() path(origPath) );

addpath( fullfile( fileparts( which( 'mrecon_kt' ) ), 'mrecon' ) )  % required for id_pnraw_data
addpath( fullfile( fileparts( which( 'mrecon_kt' ) ), 'lib', 'nifti' ) )  % required to load mask file


%% Geomtry Correction

suffixStr = '';
if ( isGeoCorrn )
    suffixStr = strcat( suffixStr, '_geocorrn' );
    reconOpts = { 'GeometryCorrection', 'Yes' };
end

isMask = false;
if exist('maskDirPath','var')
    if exist(maskDirPath,'dir')
        suffixStr     = strcat( suffixStr, '_mask' );
        isMask        = true;
    end
end


%% Output Directory

% TAR: allow to change outputDirPath
if nargin < 7
    outputDirPath  = fullfile( '/scratch/tr17/ktrecon', sprintf( 'fcmr%03i%s', fcmrNo, suffixStr ) );
end
    
if ~exist(outputDirPath,'dir')
    mkdir(outputDirPath)
end


%% Recon

diary( fullfile( outputDirPath, sprintf( 'log_mrecon_kt_hrh_%s.txt', datestr(now,'yyyymmdd_HHMMss') ) ) )

for seriesNo = seriesNos
    idStr = sprintf( 's%02i', seriesNo );
    fprintf( '\n============ %s ============\n\n', idStr )
    [ rawDataFilePath, coilSurveyFilePath, senseRefFilePath ] = id_pnraw_data( rawDataDirPath, seriesNo );
    if ( isMask )
        %maskFilePath  = fullfile( maskDirPath, [idStr '_mask_heart_dil03.nii.gz'] ); %TAR --- bugfix: changed to idStr
        maskFilePath  = fullfile( maskDirPath, ['s' num2str(seriesNo) '_mask_heart.nii.gz'] );
        fprintf( 'mask file:       %s\n', maskFilePath );
        niiMask       = load_untouch_nii( maskFilePath );
        mask          = logical(niiMask.img); %TAR --- ensure logical for mrecon_kt
        mrecon_kt_hrh( rawDataFilePath, 'senseref', senseRefFilePath, 'coilsurvey', coilSurveyFilePath, 'outputdir', outputDirPath, 'outputname', idStr, 'patchversion', patchVersion, 'reconoptionpairs', reconOpts, 'mask', mask )
    else
        mrecon_kt_hrh( rawDataFilePath, 'senseref', senseRefFilePath, 'coilsurvey', coilSurveyFilePath, 'outputdir', outputDirPath, 'outputname', idStr, 'patchversion', patchVersion, 'reconoptionpairs', reconOpts )
  end
end

diary off


end  % recon_exam(...)
