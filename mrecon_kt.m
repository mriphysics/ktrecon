function mrecon_kt( rawDataFilePath, varargin )
%MRECON_KT  k-t SENSE reconstruction using ReconFrame
%
%   MRECON_KT( rawDataFilePath, 'senseref', senseRefFilePath, 'coilsurvey', coilSurveyFilePath )
%
%   MRECON_KT( ..., 'outputdir', outputDirPath )
% 
%   MRECON_KT( ..., 'outptname', outFilePrefix )
% 
%   MRECON_KT( ..., 'ktregstrength', 0.014 )
%
%   MRECON_KT( ..., 'ktregstrengthroi', 0.00014 )
%   
%   MRECON_KT( ..., 'mask', mask )        
% 
%   MRECON_KT( ..., 'verbose', true )
%
%   NOTE: senseref and coilsurvey currently required to estimate
%   sensitivity maps; requirement may be removed in future with 
%   implementation of autocalibrated sensitivity maps

% jfpva (joshua.vanamerom@kcl.ac.uk)


%% Setup on Beastie02

%{

  run on beastie02 with pnraw/raw-ingenia mounted,

  e.g., 

      sshfs jva13@10.0.1.150:/export/pnraw/raw-ingenia ~/mnt/pnraw01-ingenia
      
      cd ~/dev_recon_from_raw
      matlab -nosplash -nodisplay -nojvm -singleCompThread
        OR
      /usr/local/MATLAB/R2012b/bin/matlab -singleCompThread &
      
      fusermount -u ~/mnt/pnraw01-ingenia/

%}


%% TODO
%
% - [ ] remove non-steady-state frames
% - [ ] optional recon baseline and/or sliding window only
% - [ ] optional save k-space data
% - [ ] optional save 2D with time offset as magnitude


%% Optional Input Argument Default Values

default.senseRefFilePath    = '';
default.coilSurveyFilePath  = '';
default.outputDirPath       = pwd;
default.outFilePrefix       = '';
default.ktRegStrength       = 0.0014;
default.ktRegStrengthROI    = default.ktRegStrength / 100;
default.mask                = [];
default.patchVersion        = '';

default.isVerbose           = false;


%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(    p, 'rawDataFilePath', ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'senseref', default.senseRefFilePath, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'coilsurvey', default.coilSurveyFilePath, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'outputdir', default.outputDirPath, ...
    @(x) validateattributes( x, {'char'}, {'nonempty','vector'}, mfilename) );

add_param_fn(   p, 'outputname', default.outFilePrefix, ...
    @(x) validateattributes( x, {'char'}, {'nonempty','vector'}, mfilename) );

add_param_fn(   p, 'ktregstrength', default.ktRegStrength, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','scalar'}, mfilename) );

add_param_fn(   p, 'ktregstrengthroi', default.ktRegStrengthROI, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','scalar'}, mfilename) );

add_param_fn(   p, 'mask', default.mask, ...
    @(x) validateattributes( x, {'logical'}, {}, mfilename) );

add_param_fn(   p, 'patchversion', default.patchVersion, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, rawDataFilePath, varargin{:} );

senseRefFilePath    = p.Results.senseref;
coilSurveyFilePath  = p.Results.coilsurvey;
outputDirPath       = p.Results.outputdir;
outFilePrefix       = p.Results.outputname;
ktRegStrength       = p.Results.ktregstrength;
ktRegStrengthROI    = p.Results.ktregstrengthroi;
mask                = p.Results.mask;
patchVersion        = p.Results.patchversion;
isVerbose           = p.Results.verbose;


%% Add Dependencies

origDepPath = path;
resetDepPath = onCleanup( @()path(origDepPath) );

codeDirPath = fileparts( which( mfilename ) );

addpath( genpath( fullfile( codeDirPath, 'lib' ) ) )
addpath( genpath( fullfile( codeDirPath, 'mrecon' ) ) )
addpath( genpath( fullfile( codeDirPath, 'recon' ) ) )
addpath( genpath( fullfile( codeDirPath, 'util' ) ) )

if isempty( which( 'MRecon' ) )
    error( 'Require ReconFrame library in Matlab path.' )
end

    
%% Anonymous Functions

% Swap dimensions to/from ReconFrame
%   ReconFrame: 
%       1-2-3-4----5---6----7----8----9---10----11----12
%       x?y-z?chan?dyn-card?echo?loca?mix?extr1?extr2?aver
%   Otherwise:
%       x-y-dyn-chan-loca-z-card-echo-mix-extr1-extr12-aver
dim.x       = 1;
dim.y       = 2;
dim.z       = 3;
dim.chan    = 4;
dim.dyn     = 5;
dim.card    = 6;
dim.echo    = 7;
dim.loca    = 8;
dim.mix     = 9; 
dim.extr1   = 10;
dim.extr2   = 11;
dim.aver    = 12;
ind_reconframe_to_xydcl = [ dim.x dim.y dim.dyn dim.chan dim.loca dim.z dim.card dim.echo dim.mix dim.extr1 dim.extr2 dim.aver ];
[~,ind_xydcl_to_reconframe] = sort( ind_reconframe_to_xydcl );
swap_dim_reconframe_to_xydcl = @( data ) permute( data, ind_reconframe_to_xydcl ); 
swap_dim_xydcl_to_reconframe = @( data ) permute( data, ind_xydcl_to_reconframe ); 

% Verbose Messages
disp_start_step_msg = @( descStr ) fprintf( '%-70s ', strcat( descStr, ' ...' ) );
disp_time_elapsed_msg = @( t ) fprintf( '... %6gs elapsed\n', round(t) );
disp_write_file_msg = @( fileStr ) fprintf( '    saved to %s\n', fileStr );


%% Setup

% Output Filename Prefix
if isempty( outFilePrefix )
    [~,outFilePrefix] = fileparts( rawDataFilePath );
end

% Coil Sensitivity Map Calculation Method
csmCalcMethod = '';
if ~isempty( senseRefFilePath ) && ~isempty( coilSurveyFilePath )
    if exist( senseRefFilePath, 'file' ) && exist( coilSurveyFilePath, 'file' )
        csmCalcMethod = 'prescan';
    else
        if ~exist( senseRefFilePath, 'file' )
            warning( 'SENSE reference file does not exist: %s', senseRefFilePath )
        end
        if ~exist( coilSurveyFilePath, 'file' )
            warning( 'Coil survey file does not exist: %s', coilSurveyFilePath )
        end 
    end
end
if ~strcmp( csmCalcMethod, 'prescan' )
    error( 'Only CSM calculation method ''prescan'' currently valid' )
end


%% Start

tic,  % start timer

fprintf( '\n%s()  started  %s\n\n', mfilename, datestr(now) );


%% Load Data

% Undersampled Data
disp_start_step_msg( 'Loading undersampled data' )
ACQ = MRecon( rawDataFilePath );
ACQ.Parameter.Parameter2Read.typ = 1;
ACQ.Parameter.Parameter2Read.mix = 0;
ACQ.Parameter.Parameter2Read.Update;
ACQ.ReadData;
mrecon_setreconparam( ACQ );
mrecon_preprocess( ACQ );
disp_time_elapsed_msg( toc )

% Training Data
disp_start_step_msg( 'Loading and preprocessing training data' )
TRN = MRecon( rawDataFilePath );
TRN.Parameter.Parameter2Read.typ = 1;
TRN.Parameter.Parameter2Read.mix = 1;
TRN.Parameter.Parameter2Read.Update;
TRN.ReadData;
mrecon_setreconparam( TRN );
mrecon_preprocess( TRN );
disp_time_elapsed_msg( toc )

% Noise
disp_start_step_msg( 'Loading and preprocessing noise data' )
NOISE = MRecon( rawDataFilePath );
NOISE.Parameter.Parameter2Read.typ = 5;
NOISE.Parameter.Parameter2Read.mix = 0;
NOISE.Parameter.Parameter2Read.Update;
NOISE.ReadData;
mrecon_setreconparam( NOISE );
mrecon_preprocess( NOISE );

% Save Data
ktAcq   = swap_dim_reconframe_to_xydcl( ACQ.Data );
ktTrn   = swap_dim_reconframe_to_xydcl( TRN.Data );
ktNoise = swap_dim_reconframe_to_xydcl( NOISE.Data );
kspaceMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_kspace.mat' ) );
save( kspaceMatFilePath, 'ktAcq', 'ktTrn', 'ktNoise' );
clear ktAcq ktTrn ktNoise

disp_time_elapsed_msg( toc )

disp_write_file_msg( kspaceMatFilePath )


%% Calculate Coil Sensitivity Maps

switch csmCalcMethod
    
    case 'prescan'
        
        % Sense Reference
        disp_start_step_msg( 'Loading and preprocessing SENSE reference data' ),
        SREF = MRecon( senseRefFilePath );
        %         SREF.Parameter.Parameter2Read.typ = 1;
        %         SREF.Parameter.Parameter2Read.Update;
        %         SREF.ReadData;
        %         mrecon_preprocess( SREF );
        disp_time_elapsed_msg( toc ),
        
        % Coil Survey    
        disp_start_step_msg( 'Loading and preprocessing coil survey data' ),
        COIL = MRecon( coilSurveyFilePath );
        %         COIL.Parameter.Parameter2Read.typ = 1;
        %         COIL.Parameter.Parameter2Read.Update;
        %         COIL.ReadData;
        %         mrecon_preprocess( COIL );
        disp_time_elapsed_msg( toc ), 
            
        % Target
        TGT = MRecon( rawDataFilePath );
        TGT.Parameter.Parameter2Read.typ = 1;
        TGT.Parameter.Parameter2Read.mix = 0;
        TGT.Parameter.Recon.RemoveMOversampling = 'No';       
        
        % Calculate Coil Sensitivity Maps
        disp_start_step_msg( 'Calculating coil sensitivity maps' )
        
        % Create SENSE object
        SENS = MRsense( SREF, TGT, COIL );

        % Calculate sensitivity maps
        SENS.CalculateSensitivity   = 1;
        SENS.Mask                   = 1;
        SENS.Smooth                 = 1;
        SENS.Extrapolate            = 1;
        SENS.MatchTargetSize        = 1;
        SENS.RemoveMOversampling    = 0;
        SENS.OutputSizeSensitivity  = [ size(ACQ.Data,dim.x), size(ACQ.Data,dim.y), size(ACQ.Data,dim.loca) ];
        SENS.OutputSizeReformated   = SENS.OutputSizeSensitivity;
        
        SENS.Perform;

        % Extract data from MRsense object
        csm     = swap_dim_reconframe_to_xydcl( SENS.Sensitivity );             % coil sensitivity maps
        imBody  = swap_dim_reconframe_to_xydcl( SENS.ReformatedBodycoilData );  % body coil image
        imCoil  = swap_dim_reconframe_to_xydcl( SENS.ReformatedCoilData );      % array coil images
        psiSens = SENS.Psi;                                 % array coil noise covariance
        csmMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_csm.mat' ) );
        save( csmMatFilePath, 'csm', 'imBody', 'imCoil', 'psiSens', '-v7.3' );
        
        disp_time_elapsed_msg( toc )
        
        disp_write_file_msg( csmMatFilePath )
        
    otherwise

        error( 'Unrecognised CSM calculation method: %s', csmCalcMethod ),
        
end

% Add Sensitivity Maps to Undersampled and Training Data MRecon Objects
ACQ.Parameter.Recon.Sensitivities = SENS;
TRN.Parameter.Recon.Sensitivities = SENS;


%% TODO: Drop Non-Steady State Frames

%{
if isempty( numFrame ),
    numFrame = size(ktAcq,3);  % e.g., s = squeeze(sum(sum(bsxfun(@times,mask,abs(xtSlw)))));  figure,  plot(s),  grid on,  
end                            % or,   frameNo = find_steadystate_framenos( permute(abs(xtSlw),[1,2,4,3]), 'showFig' ); 
                               %       ktFactor = 8; 
                               %       numFrame = floor(length(frameNo)/ktFactor)*ktFactor;

ktAcq = ktAcq(:,:,(end-numFrame+1):end,:);
ktTrn = ktTrn(:,:,1:numFrame,:);
%}


%% Baseline Recon

% Display Start Message
disp_start_step_msg( 'Reconstructing baseline images' ),

% Copy MRecon Object
BLN = ACQ.Copy;
BLN.Parameter.Recon.SENSE   = 'Yes';  % CLEAR coil combination

% Get Data from MRecon Object
ktAcq = swap_dim_reconframe_to_xydcl( BLN.Data );

% Sampling Pattern
ktSmp = single( sum( sum( ktAcq, 4 ), 1 ) ~= 0 );

% Get Baseline
ktBln = bsxfun( @rdivide, sum( ktAcq, 3 ), sum( ktSmp, 3 ) );

% Put Data in Back in MRecon Object
BLN.Data = swap_dim_xydcl_to_reconframe( ktBln );

% Transform to Image Space
mrecon_k2i( BLN );

% Image Space Postprocessing
mrecon_postprocess( BLN );

% Write to NIfTI
blnNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_bln_mag' ) );
blnNiiFilePath = mrecon_writenifti( BLN, blnNiiFilePath );

% Display Time Elapsed Message
disp_time_elapsed_msg( toc )
disp_write_file_msg( blnNiiFilePath )


%% Sliding Window Recon

% Display Start Message
disp_start_step_msg( 'Computing sliding window reconstruction of undersampled data' ),

% Sliding Window Reconstruction MRecon Object
SLW = ACQ.Copy;
SLW.Parameter.Recon.SENSE   = 'Yes';  % CLEAR coil combination

% Get Data from MRecon Object
ktSlw = swap_dim_reconframe_to_xydcl( SLW.Data );

% Fill k-Space Using Sliding Window for Each Slice
for iSlice = 1:size(SLW.Data,dim.loca)
    ktSlw(:,:,:,:,iSlice) = kt_sliding_window( ktSlw(:,:,:,:,iSlice) );
end

% Put Data in Back in MRecon Object
SLW.Data = swap_dim_xydcl_to_reconframe( ktSlw );

% Recon Images
mrecon_k2i( SLW )

% Postprocessing
mrecon_postprocess( SLW );

% Get Timing
frameDuration = mrecon_calc_frame_duration( SLW );

% Save as NIfTI
slwNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_slw_mag' ) );
slwNiiFilePath = mrecon_writenifti( SLW, slwNiiFilePath, 'frameduration', frameDuration );

% Display Time Elapsed Message
disp_time_elapsed_msg( toc ),

disp_write_file_msg( slwNiiFilePath )


%% k-t SENSE Reconstruction

% Display Start Message
fprintf( 'Performing k-t SENSE reconstruction  \n' ),

% k-t SENSE Reconstruction MRecon Object
RCN = ACQ.Copy;
RCN.Data = sum (RCN.Data, dim.chan );

% k-t SENSE Prior MRecon Object
PRI = TRN.Copy;
PRI.Parameter.Recon.SENSE = 'Yes';  % CLEAR coil combination

% DC/Baseline MRecon Object
DC = ACQ.Copy;
DC.Data = sum( sum( DC.Data, dim.chan ), dim.dyn );

% Process Slice-by-Slice
numSlice = size( RCN.Data, dim.loca );
for iSlice = 1:numSlice
    
    % Display Start Message
    disp_start_step_msg( sprintf( '    slice %3i', iSlice ) ) 
    
    % Get Data from MRecon Objects
    ktAcq       = swap_dim_reconframe_to_xydcl( ACQ.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:) );
    ktTrn       = swap_dim_reconframe_to_xydcl( TRN.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:) );
    csm         = swap_dim_reconframe_to_xydcl( ACQ.Parameter.Recon.Sensitivities.Sensitivity(:,:,:,:,:,:,:,iSlice,:,:,:,:) );
    noiseCov    = SENS.Psi; 
        %kNoise      = swap_dim_reconframe_to_xydcl( NOISE.Data(:,:,:,:,:,:,:,iLoc,:,:,:,:) ); 
        
    % k-t SENSE Reconstruction
    if isempty( mask )
        [ ktRcn, ktDC ] = recon_ktsense( ktAcq, ktTrn, csm, 'noisecov', noiseCov, 'lambda0', ktRegStrength );
    else
        [ ktRcn, ktDC ] = recon_ktsense( ktAcq, ktTrn, csm, 'noisecov', noiseCov, 'lambda0', ktRegStrength,  'mask', mask, 'lambdaroi', ktRegStrengthROI ); 
    end
    % TODO: compare recon using noise covariance   
    % 1) estimated from k-t undersampled data, 
    % 2) calculated using noise samples in MRecon object NOISE, and
    % 3) calculated in sensitivity maps MRecon object SENS
    
    % Zero-Pad k-t Data to Fit MRecon Object Prior to RemoveOversampling()
    % FIXME: zero-pad in image space
    %{
        ktRcn = padarray( ktRcn, [round((size(RCN.Data,dim.x)-size(ktRcn,1))/2),0,0,0,0,0,0,0,0,0,0,0] );
        ktDC  = padarray( ktDC,  [round((size(DC.Data,dim.x)-size(ktDC,1))/2),  0,0,0,0,0,0,0,0,0,0,0] );
    %}
    
    % Put Data in Back in MRecon Objects
    RCN.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:) = swap_dim_xydcl_to_reconframe( ktRcn );
    DC.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:)  = swap_dim_xydcl_to_reconframe( ktDC );

    % Display Time Elapsed Message
    disp_time_elapsed_msg( toc ),
    
end

% Display Start Message
disp_start_step_msg( 'Finalising k-t SENSE reconstruction' ),

% Recon Images
mrecon_k2i( RCN )
mrecon_k2i( PRI )
mrecon_k2i( DC )

% Postprocessing
mrecon_postprocess( RCN );
mrecon_postprocess( PRI );
mrecon_postprocess( DC );

% Get Timing
frameDuration = mrecon_calc_frame_duration( RCN );

% Save as NIfTI
rltNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_mag' ) );
rltNiiFilePath = mrecon_writenifti( RCN, rltNiiFilePath, 'frameduration', frameDuration );
priNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_pri_mag' ) );
priNiiFilePath = mrecon_writenifti( PRI, priNiiFilePath, 'frameduration', frameDuration );
dcNiiFilePath  = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_mag' ) );
dcNiiFilePath  = mrecon_writenifti( DC, dcNiiFilePath );
disp_write_file_msg( rltNiiFilePath )
disp_write_file_msg( priNiiFilePath )
disp_write_file_msg( dcNiiFilePath )


%% Separate Slices

% Save m2d Stack as NIfTI
rltCpxNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_cpx' ) );
rltCpxNiiFilePath = mrecon_writenifti( RCN, rltCpxNiiFilePath, 'complex', true, 'frameduration', frameDuration );
rltCpxNiiFileName = filename( rltCpxNiiFilePath );

% Get Frame Duration
frameDuration = mrecon_calc_frame_duration( RCN );

% Get Study Start Time
tStudyStr = RCN.Parameter.GetValue('RFR_ECSERIES_DICOM_SERIES_TIME');
ind = strfind( tStudyStr, '.' );
hr  = str2double( tStudyStr(1:(ind-5)) );
min = str2double( tStudyStr((ind-4):(ind-3)) );
sec = str2double( tStudyStr((ind-2):end) );
tStudy = 3600*hr + 60*min + sec;
% Get Series Start Time
tSeriesStr = RCN.Parameter.GetValue('RFR_SERIES_DICOM_SERIES_TIME');
ind = strfind( tSeriesStr, '.' );
hr  = str2double( tSeriesStr(1:(ind-5)) );
min = str2double( tSeriesStr((ind-4):(ind-3)) );
sec = str2double( tSeriesStr((ind-2):end) );
tSeries = 3600*hr + 60*min + sec;

% Get Slice Duration
seriesDuration  = RCN.Parameter.GetValue('AC_total_scan_time');
numSlice        = size( RCN.Data, dim.loca );
numDynDummy     = RCN.Parameter.GetValue( 'MP_nr_dummy_dynamic_scans' );
switch patchVersion
    % NOTE:
    case 'c51c8c1'
        % acquition order: dummy, acq_1, dummy, acq_2, ..., dummy, acq_n
        %                  dummy, trn_1, dummy, trn_2, ..., dummy, trn_n
        sliceDuration    = seriesDuration / numSlice;
        sliceStartOffset = numDynDummy * frameDuration;
    otherwise
        fprintf( '    Warning: patch version (%s) unspecified or unrecognised; slice timing may be incorrect.', patchVersion )
        sliceDuration    = seriesDuration / numSlice;
        sliceStartOffset = numDynDummy * frameDuration;
end

%{
% NOTE: datetime isn't available in earlier Matlab versions
timestr2datetime = @(str) datetime( str, 'InputFormat', 'HHmmss.SSSSS' );
datetime2seconds = @(t) 360*t.Hour + 60*t.Minute + t.Second;
tStudy  = datetime2seconds( timestr2datetime( RCN.Parameter.GetValue('RFR_ECSERIES_DICOM_SERIES_TIME') ) );
tSeries = datetime2seconds( timestr2datetime( RCN.Parameter.GetValue('RFR_SERIES_DICOM_SERIES_TIME') ) );
%}

% Generate Script to Extract Slices from Stack
scriptFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_extract_slices.sh' ) );
fid = fopen( scriptFilePath, 'w' );

fwrite( fid, sprintf( '#!/bin/sh\n\n' ) );
fwrite( fid, sprintf( '# %-25s %-15s %s\n', 'study start time', '(hhmmss.sss)', tStudyStr ) );
fwrite( fid, sprintf( '# %-25s %-15s %g\n', 'study start time', '(s)', tStudy ) );
fwrite( fid, sprintf( '# %-25s %-15s %s\n', 'series start time', '(hhmmss.sss)', tSeriesStr ) );
fwrite( fid, sprintf( '# %-25s %-15s %g\n', 'series start time', '(s)', tSeries ) );
fwrite( fid, sprintf( '# %-25s %-15s %g\n', 'series duration', '(s)', seriesDuration ) );
fwrite( fid, sprintf( '# %-25s %-15s %g\n', 'slice durationn', '(s)', sliceDuration ) );
fwrite( fid, sprintf( '# %-25s %-15s %g\n', 'dummy delay', '(s)', sliceStartOffset ) );
fwrite( fid, sprintf( '# %-25s %-15s %g\n', 'frame duration', '(s)', frameDuration ) );

% Process Slice-by-Slice
for iSlice = 1:numSlice

    fwrite( fid, sprintf( '\n# Slice %i\n\n', iSlice ) );
    
    % Extract Slice from M2D Stack
    rltCpxSliceNiiFileName = sprintf( '%s_rlt_slice%02i_cpx.nii.gz', outFilePrefix, iSlice );
    cmd = sprintf( 'region %s %s -Rz1 %i -Rz2 %i\n\n', rltCpxNiiFileName, rltCpxSliceNiiFileName, iSlice-1, iSlice );
    fwrite( fid, cmd );
    
    % Set Slice Thickness and Start Time Relative to Start of Study
    sliceThickness = RCN.Parameter.Scan.RecVoxelSize(3);
    tSlice = tSeries + (iSlice-1) * sliceDuration + sliceStartOffset - tStudy;
    cmd = sprintf( 'headertool %s %s -size 0 0 %f -timeOrigin %f\n\n', rltCpxSliceNiiFileName, rltCpxSliceNiiFileName, sliceThickness, tSlice );
    fwrite( fid, cmd );

end

% Close Script File
fclose( fid );

% Display Time Elapsed Message
disp_time_elapsed_msg( toc ),

disp_write_file_msg( rltCpxNiiFilePath )

disp_write_file_msg( scriptFilePath )


%% End

fprintf( '\n%s()  finished  %s\n\n', mfilename, datestr(now) );


end  % mrecon_kt(...)
