function mrecon_kt( rawDataFilePath, varargin )
%MRECON_KT  k-t SENSE reconstruction using ReconFrame
%
%   MRECON_KT( rawDataFilePath, 'senseref', senseRefFilePath, 'coilsurvey', coilSurveyFilePath )
%
%   MRECON_KT( ..., 'outputdir', outputDirPath )
% 
%   MRECON_KT( ..., 'outptname', outFilePrefix )
%
%   MRECON_KT( ..., 'patchversion', 'PIH1' )
%
%   MRECON_KT( ..., 'reconoptionpairs', opts )
%
%   MRECON_KT( ..., 'reconktsense', false )
%
%   MRECON_KT( ..., 'ktregstrength', 0.014 )
%
%   MRECON_KT( ..., 'ktregstrengthroi', 0.00014 )
%   
%   MRECON_KT( ..., 'mask', mask )
%
%   MRECON_KT( ..., 'exportpdf', true )
%
%   MRECON_KT( ..., 'exportjson', true )
%
%   MRECON_KT( ..., 'exportgoalc', true )
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

      ```bash
      # Mount Raw Data Drive
      sshfs jva13@10.0.1.150:/export/pnraw/raw-ingenia ~/mnt/pnraw01-ingenia
      # Open Matlab
      matlab -nosplash -nodisplay -nojvm -singleCompThread
        OR
      /usr/local/MATLAB/R2012b/bin/matlab -singleCompThread &
      ```
        
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
      % Exit Matlab
      exit
      ```

      ```bash
      # Dismount Raw Data Drive
      fusermount -u ~/mnt/pnraw01-ingenia/
      ```

%}


%% Optional Input Argument Default Values

default.senseRefFilePath    = '';
default.coilSurveyFilePath  = '';
default.outputDirPath       = pwd;
default.outFilePrefix       = '';
default.reconOpts           = {};
default.isReconKtSense      = true;
default.ktRegStrength       = 0.0014;
default.ktRegStrengthROI    = default.ktRegStrength / 100;
default.mask                = []; 
default.patchVersion        = '';
default.isExportPdf         = true;
default.isExportJson        = true;
default.isExportGoalc       = true;
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

add_param_fn(   p, 'reconoptionpairs', default.reconOpts, ...
    @(x) validateattributes( x, {'cell'}, {}, mfilename) );

add_param_fn(   p, 'reconktsense', default.isReconKtSense, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

add_param_fn(   p, 'ktregstrength', default.ktRegStrength, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','scalar'}, mfilename) );

add_param_fn(   p, 'ktregstrengthroi', default.ktRegStrengthROI, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','scalar'}, mfilename) );

add_param_fn(   p, 'mask', default.mask, ...
    @(x) validateattributes( x, {'logical'}, {}, mfilename) );

add_param_fn(   p, 'patchversion', default.patchVersion, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'exportpdf', default.isExportPdf, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

add_param_fn(   p, 'exportjson', default.isExportJson, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

add_param_fn(   p, 'exportgoalc', default.isExportGoalc, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, rawDataFilePath, varargin{:} );

senseRefFilePath    = p.Results.senseref;
coilSurveyFilePath  = p.Results.coilsurvey;
outputDirPath       = p.Results.outputdir;
outFilePrefix       = p.Results.outputname;
reconOpts           = p.Results.reconoptionpairs;
isReconKtSense      = p.Results.reconktsense;
ktRegStrength       = p.Results.ktregstrength;
ktRegStrengthROI    = p.Results.ktregstrengthroi;
mask                = p.Results.mask;
patchVersion        = p.Results.patchversion;
isExportPdf         = p.Results.exportpdf;
isExportJson        = p.Results.exportjson;
isExportGoalc       = p.Results.exportgoalc;
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
disp_start_step_msg( 'Loading and preprocessing undersampled data' )
ACQ = MRecon( rawDataFilePath );
ACQ.Parameter.Parameter2Read.typ = 1;
ACQ.Parameter.Parameter2Read.mix = 0;
ACQ.Parameter.Parameter2Read.Update;
ACQ.ReadData;
mrecon_setreconparam( ACQ, 'optionpairs', reconOpts );
mrecon_preprocess( ACQ );
disp_time_elapsed_msg( toc )

% Training Data
disp_start_step_msg( 'Loading and preprocessing training data' )
TRN = MRecon( rawDataFilePath );
TRN.Parameter.Parameter2Read.typ = 1;
TRN.Parameter.Parameter2Read.mix = 1;
TRN.Parameter.Parameter2Read.Update;
TRN.ReadData;
mrecon_setreconparam( TRN, 'optionpairs', reconOpts );
mrecon_preprocess( TRN );
disp_time_elapsed_msg( toc )

% Noise
disp_start_step_msg( 'Loading and preprocessing noise data' )
NOISE = MRecon( rawDataFilePath );
NOISE.Parameter.Parameter2Read.typ = 5;
NOISE.Parameter.Parameter2Read.mix = 0;
NOISE.Parameter.Parameter2Read.Update;
NOISE.ReadData;
mrecon_setreconparam( NOISE, 'optionpairs', reconOpts );
mrecon_preprocess( NOISE );
disp_time_elapsed_msg( toc )

% Save Data
disp_start_step_msg( 'Saving k-space data' )
ktAcq   = swap_dim_reconframe_to_xydcl( ACQ.Data );
ktTrn   = swap_dim_reconframe_to_xydcl( TRN.Data );
ktNoise = swap_dim_reconframe_to_xydcl( NOISE.Data );
kspaceMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_kspace.mat' ) );
save( kspaceMatFilePath, 'ktAcq', 'ktTrn', 'ktNoise', '-v7.3' );
clear ktAcq ktTrn ktNoise
disp_time_elapsed_msg( toc )
disp_write_file_msg( kspaceMatFilePath )


%% Calculate Coil Sensitivity Maps

switch csmCalcMethod
    
    case 'prescan'
        
        % Sense Reference
        disp_start_step_msg( 'Loading SENSE reference data' ),
        SREF = MRecon( senseRefFilePath );
        disp_time_elapsed_msg( toc ),
        
        % Coil Survey    
        disp_start_step_msg( 'Loading coil survey data' ),
        COIL = MRecon( coilSurveyFilePath );
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
        psiSens = SENS.Psi;                                                     % array coil noise covariance
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
blnAbNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_bln_ab' ) );
blnAbNiiFilePath = mrecon_writenifti( BLN, blnAbNiiFilePath, 'datatype', 'magnitude' );
blnReNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_bln_re' ) );
blnReNiiFilePath = mrecon_writenifti( BLN, blnReNiiFilePath, 'datatype', 'real' );
blnImNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_bln_im' ) );
blnImNiiFilePath = mrecon_writenifti( BLN, blnImNiiFilePath, 'datatype', 'imaginary' );

% Save as .mat
xtBln = BLN.Data;
blnMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_bln_recon.mat' ) );
save( blnMatFilePath, 'xtBln', '-v7.3' );
clear xtBln 

% Display Time Elapsed Message
disp_time_elapsed_msg( toc )
disp_write_file_msg( blnAbNiiFilePath )
disp_write_file_msg( blnReNiiFilePath )
disp_write_file_msg( blnImNiiFilePath )
disp_write_file_msg( blnMatFilePath )


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
slwNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_slw_ab' ) );
slwNiiFilePath = mrecon_writenifti( SLW, slwNiiFilePath, 'frameduration', frameDuration );

% Display Time Elapsed Message
disp_time_elapsed_msg( toc ),

disp_write_file_msg( slwNiiFilePath )


%% k-t SENSE Reconstruction

if ( isReconKtSense )
    
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
        disp_start_step_msg( sprintf( '  slice %3i', iSlice ) )
        
        % Get Data from MRecon Objects
        ktAcq       = swap_dim_reconframe_to_xydcl( ACQ.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:) );
        ktTrn       = swap_dim_reconframe_to_xydcl( TRN.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:) );
        csm         = swap_dim_reconframe_to_xydcl( ACQ.Parameter.Recon.Sensitivities.Sensitivity(:,:,:,:,:,:,:,iSlice,:,:,:,:) );
        noiseCov    = SENS.Psi;
        
        % k-t SENSE Reconstruction
        if isempty( mask )
            [ ktRcn, ktDc ] = recon_ktsense( ktAcq, ktTrn, csm, 'noisecov', noiseCov, 'lambda0', ktRegStrength );
        else
            [ ktRcn, ktDc ] = recon_ktsense( ktAcq, ktTrn, csm, 'noisecov', noiseCov, 'lambda0', ktRegStrength,  'mask', mask(:,:,iSlice), 'lambdaroi', ktRegStrengthROI );
        end
        % TODO: compare recon using available noise covariance estimates:
        %   1) estimated from k-t undersampled data,
        %   2) calculated using noise samples in MRecon object NOISE, and
        %   3) calculated in sensitivity maps MRecon object SENS
        
        % Put Data in Back in MRecon Objects
        RCN.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:) = swap_dim_xydcl_to_reconframe( ktRcn );
        DC.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:)  = swap_dim_xydcl_to_reconframe( ktDc );
        
        % Save Reconstructed Data
        reconMatFilePath = fullfile( outputDirPath, sprintf( '%s_slice%02i_recon.mat', outFilePrefix, iSlice ) );
        save( reconMatFilePath, 'ktRcn', 'ktDc', '-v7.3' );
        clear ktRcn ktDc
        
        % Display Time Elapsed Message
        disp_time_elapsed_msg( toc ),
        disp_write_file_msg( reconMatFilePath ),
        
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
    
    % Write Real-Time as NIfTI
    rltAbNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_ab' ) );
    rltAbNiiFilePath = mrecon_writenifti( RCN, rltAbNiiFilePath, 'frameduration', frameDuration, 'datatype', 'magnitude' );
    rltReNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_re' ) );
    rltReNiiFilePath = mrecon_writenifti( RCN, rltReNiiFilePath, 'frameduration', frameDuration, 'datatype', 'real' );
    rltImNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_im' ) );
    rltImNiiFilePath = mrecon_writenifti( RCN, rltImNiiFilePath, 'frameduration', frameDuration, 'datatype', 'imaginary' );
    
    % Write Training as NIfTI
    priAbNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_trn_ab' ) );
    priAbNiiFilePath = mrecon_writenifti( PRI, priAbNiiFilePath, 'frameduration', frameDuration, 'datatype', 'magnitude' );
    priReNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_trn_re' ) );
    priReNiiFilePath = mrecon_writenifti( PRI, priReNiiFilePath, 'frameduration', frameDuration, 'datatype', 'real' );
    priImNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_trn_im' ) );
    priImNiiFilePath = mrecon_writenifti( PRI, priImNiiFilePath, 'frameduration', frameDuration, 'datatype', 'imaginary' );
    
    % Write DC as NIfTI
    dcAbNiiFilePath  = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_ab' ) );
    dcAbNiiFilePath  = mrecon_writenifti( DC, dcAbNiiFilePath, 'datatype', 'magnitude' );
    dcReNiiFilePath  = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_re' ) );
    dcReNiiFilePath  = mrecon_writenifti( DC, dcReNiiFilePath, 'datatype', 'real' );
    dcImNiiFilePath  = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_im' ) );
    dcImNiiFilePath  = mrecon_writenifti( DC, dcImNiiFilePath, 'datatype', 'imaginary' );
    
    % Save as .mat
    xtRcn = RCN.Data;
    rltMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_recon.mat' ) );
    save( rltMatFilePath, 'xtRcn', '-v7.3' );
    clear xtRcn
    xtPri = PRI.Data;
    priMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_pri_recon.mat' ) );
    save( priMatFilePath, 'xtPri', '-v7.3' );
    clear xtPri
    xtDc = DC.Data;
    dcMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_recon.mat' ) );
    save( dcMatFilePath, 'xtDc', '-v7.3' );
    clear xtDc
    
    % Save Parameters
    PARAM = mrecon_getparameters( RCN );
    PARAM.Timing = mrecon_getslicetiming( RCN, 'patchversion', patchVersion );
    rltParamMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_parameters.mat' ) );
    save( rltParamMatFilePath, 'PARAM', '-v7.3' );
    
    % Display Time Elapsed Message
    disp_time_elapsed_msg( toc ),
    disp_write_file_msg( rltAbNiiFilePath )
    disp_write_file_msg( rltReNiiFilePath )
    disp_write_file_msg( rltImNiiFilePath )
    disp_write_file_msg( rltParamMatFilePath )
    disp_write_file_msg( rltMatFilePath )
    disp_write_file_msg( priAbNiiFilePath )
    disp_write_file_msg( priReNiiFilePath )
    disp_write_file_msg( priImNiiFilePath )
    disp_write_file_msg( priMatFilePath )
    disp_write_file_msg( dcAbNiiFilePath )
    disp_write_file_msg( dcReNiiFilePath )
    disp_write_file_msg( dcImNiiFilePath )
    disp_write_file_msg( dcMatFilePath )
    
else
    
    fprintf( 'Skipping k-t SENSE reconstruction  \n' ),
    
end


%% Save GVE PDF

if ( isExportPdf )
    
    % Display Start Message
    disp_start_step_msg( 'Generating protocol definition file' ),
    
    % Save GVE
    gvepdfFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_pdf.gve' ) );
    ACQ.Parameter.ExtractPDFFile( gvepdfFilePath );
    
    % Display Time Elapsed Message
    disp_time_elapsed_msg( toc ),
    disp_write_file_msg( gvepdfFilePath )

end


%% Export JSON

if ( isExportJson )
    
    % Display Start Message
    disp_start_step_msg( 'Generating JSON parameter file' ),
    
    % Export JSON
    jsonFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '.json' ) );
    if ( isReconKtSense )
        RCN.Parameter.Export2Json( jsonFilePath );
    else
        ACQ.Parameter.Export2Json( jsonFilePath );
    end
    
    % Display Time Elapsed Message
    disp_time_elapsed_msg( toc ),
    disp_write_file_msg( jsonFilePath )
    
end


%% Export GoalC

if ( isExportGoalc )
    
    % Display Start Message
    disp_start_step_msg( 'Exporting GoalC information file' ),
    
    % Export GoalC
    goalcFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_goalc.txt' ) );
    mrecon_writegoalc2txt( ACQ, goalcFilePath );
    
    % Display Time Elapsed Message
    disp_time_elapsed_msg( toc ),
    disp_write_file_msg( goalcFilePath )
    
end


%% End

fprintf( '\n%s()  finished  %s\n\n', mfilename, datestr(now) );


end  % mrecon_kt(...)
