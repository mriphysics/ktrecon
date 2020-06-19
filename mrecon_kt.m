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


%% Optional Input Argument Default Values

default.senseRefFilePath    = '';
default.coilSurveyFilePath  = '';
default.outputDirPath       = pwd;
default.outFilePrefix       = '';
default.cusTrnDirPath       = ''; %TAR
default.isSelfCaliPreProc   = false; %TAR - outputDirPath must contain s*_recon_sc_slw.mat
default.reconOpts           = {};
default.isReconKtSense      = true;
default.isSweepAcq          = false; % TAR
default.swpKrnFullWidth     = 96; % TAR
default.swpKrnSpacing       = []; % TAR
default.swpKrnLoca          = []; % TAR
default.ktRegStrength       = 0.0014;
default.ktRegStrengthROI    = default.ktRegStrength / 100;
default.mask                = [];
default.makeHarmonicFilter  = false; % TAR
default.patchVersion        = 'PIH1';
default.isExportPdf         = true;
default.isExportJson        = false; % TAR - Note: not available in MRecon 515
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

add_param_fn(   p, 'cusTrnDirPath', default.cusTrnDirPath, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) ); %TAR

add_param_fn(   p, 'isSelfCaliPreProc', default.isSelfCaliPreProc, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) ); %TAR

add_param_fn(   p, 'reconoptionpairs', default.reconOpts, ...
    @(x) validateattributes( x, {'cell'}, {}, mfilename) );

add_param_fn(   p, 'reconktsense', default.isReconKtSense, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

add_param_fn(   p, 'reconktsweep', default.isSweepAcq, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) ); % TAR

add_param_fn(   p, 'sweepkernelwidth', default.swpKrnFullWidth, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','scalar'}, mfilename) ); % TAR

add_param_fn(   p, 'sweepkernelspacing', default.swpKrnSpacing, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','scalar'}, mfilename) ); % TAR

add_param_fn(   p, 'sweepkernellocations', default.swpKrnLoca, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','vector'}, mfilename) ); % TAR

add_param_fn(   p, 'ktregstrength', default.ktRegStrength, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','scalar'}, mfilename) );

add_param_fn(   p, 'ktregstrengthroi', default.ktRegStrengthROI, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','scalar'}, mfilename) );

add_param_fn(   p, 'mask', default.mask, ...
    @(x) validateattributes( x, {'logical'}, {}, mfilename) );

add_param_fn(   p, 'makeHarmonicFilter', default.makeHarmonicFilter, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

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
cusTrnDirPath       = p.Results.cusTrnDirPath; %TAR
isSelfCaliPreProc   = p.Results.isSelfCaliPreProc; %TAR
reconOpts           = p.Results.reconoptionpairs;
isReconKtSense      = p.Results.reconktsense;
isSweepAcq          = p.Results.reconktsweep; % TAR
swpKrnFullWidth     = p.Results.sweepkernelwidth; % TAR
swpKrnSpacing       = p.Results.sweepkernelspacing; % TAR
swpKrnLoca          = p.Results.sweepkernellocations; % TAR
ktRegStrength       = p.Results.ktregstrength;
ktRegStrengthROI    = p.Results.ktregstrengthroi;
mask                = p.Results.mask;
makeHarmonicFilter  = p.Results.makeHarmonicFilter; %TAR
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

% warning('kspace.mat NOT saved.');
% save( kspaceMatFilePath, 'ktAcq', 'ktTrn', 'ktNoise', '-v7.3' );

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

% % Display Start Message
% disp_start_step_msg( 'Reconstructing baseline images' ),
% 
% % Copy MRecon Object
% BLN = ACQ.Copy;
% BLN.Parameter.Recon.SENSE   = 'Yes';  % CLEAR coil combination
% 
% % Get Data from MRecon Object
% ktAcq = swap_dim_reconframe_to_xydcl( BLN.Data );
% 
% % Sampling Pattern
% ktSmp = single( sum( sum( ktAcq, 4 ), 1 ) ~= 0 );
% 
% % Get Baseline
% ktBln = bsxfun( @rdivide, sum( ktAcq, 3 ), sum( ktSmp, 3 ) );
% 
% % Put Data in Back in MRecon Object
% BLN.Data = swap_dim_xydcl_to_reconframe( ktBln );
% 
% % Transform to Image Space
% mrecon_k2i( BLN );
% 
% % Image Space Postprocessing
% mrecon_postprocess( BLN );
% 
% % Write to NIfTI
% blnAbNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_bln_ab' ) );
% blnAbNiiFilePath = mrecon_writenifti( BLN, blnAbNiiFilePath, 'datatype', 'magnitude' );
% blnReNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_bln_re' ) );
% blnReNiiFilePath = mrecon_writenifti( BLN, blnReNiiFilePath, 'datatype', 'real' );
% blnImNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_bln_im' ) );
% blnImNiiFilePath = mrecon_writenifti( BLN, blnImNiiFilePath, 'datatype', 'imaginary' );
% 
% % Save as .mat
% xtBln = BLN.Data;
% blnMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_bln_recon.mat' ) );
% save( blnMatFilePath, 'xtBln', '-v7.3' );
% clear xtBln 
% 
% % Display Time Elapsed Message
% disp_time_elapsed_msg( toc )
% disp_write_file_msg( blnAbNiiFilePath )
% disp_write_file_msg( blnReNiiFilePath )
% disp_write_file_msg( blnImNiiFilePath )
% disp_write_file_msg( blnMatFilePath )


%% Sliding Window Recon

% % Display Start Message
% disp_start_step_msg( 'Computing sliding window reconstruction of undersampled data' ),
% 
% % Sliding Window Reconstruction MRecon Object
% SLW = ACQ.Copy;
% SLW.Parameter.Recon.SENSE   = 'Yes';  % CLEAR coil combination
% 
% % Get Data from MRecon Object
% ktSlw = swap_dim_reconframe_to_xydcl( SLW.Data );
% 
% % Fill k-Space Using Sliding Window for Each Slice
% for iSlice = 1:size(SLW.Data,dim.loca)
%     ktSlw(:,:,:,:,iSlice) = kt_sliding_window( ktSlw(:,:,:,:,iSlice) );
% end
% 
% % Put Data in Back in MRecon Object
% SLW.Data = swap_dim_xydcl_to_reconframe( ktSlw );
% 
% % Recon Images
% mrecon_k2i( SLW )
% 
% % Postprocessing
% mrecon_postprocess( SLW );
% 
% % Get Timing
% frameDuration = mrecon_calc_frame_duration( SLW );
% 
% % Save as NIfTI
% slwNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_slw_ab' ) );
% slwNiiFilePath = mrecon_writenifti( SLW, slwNiiFilePath, 'frameduration', frameDuration );
% 
% % Display Time Elapsed Message
% disp_time_elapsed_msg( toc ),
% disp_write_file_msg( slwNiiFilePath )


%% Self-calibrated Sliding Window Preprocessing

if isSelfCaliPreProc
   
    % Display Start Message
    disp_start_step_msg( 'Computing sliding window reconstruction of undersampled data' ),
    
    % Self-Calibrated Sliding Window Reconstruction MRecon Object
    SCSLW = ACQ.Copy;
    SCSLW.Parameter.Recon.Sensitivities = SENS;
    
    % Perform preprocessing
    kt_sc_slw_preproc( outputDirPath, outFilePrefix, SCSLW );
    
    % Display Time Elapsed Message
    disp_time_elapsed_msg( toc ),
    
end


%% Use Custom Training Data
xtTrnCus = [];
if ~isempty( cusTrnDirPath )
    
    disp_start_step_msg( 'Loading Custom Training data' ),
    cusTrnMatFilePath = fullfile( cusTrnDirPath, strcat( outFilePrefix, '_rlt_sc_slw_recon.mat' ) );            
    load( cusTrnMatFilePath ); 
    xtTrnCus = xtRcn; clear xtRcn
    disp_time_elapsed_msg( toc ),
    
end


%% Use x-f Harmonic Filter + Low Pass Filter
if makeHarmonicFilter == true
    
    stackFrameDuration = mrecon_calc_frame_duration( ACQ );
    ktRegStrengthROI = ktRegStrength / 100; % denominator = alpha
    fprintf( 'Using x-f Harmonic Filter with ktRegStrengthROI = %02i \n', ktRegStrengthROI ),
    fprintf( 'Frame duration = %02i \n', stackFrameDuration ),
    
elseif makeHarmonicFilter == false
    
    stackFrameDuration = [];

end


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
    
    numSlice = size( RCN.Data, dim.loca );
    
    % Partition Sweep Data
    if ( isSweepAcq )
        
        fprintf( 'Reconstructing Sweep Acquisition ... \n' ),
        
        nX   = size( ACQ.Data, dim.x );
        nY   = size( ACQ.Data, dim.y );
        nZ   = size( ACQ.Data, dim.loca );
        nT   = size( ACQ.Data, dim.dyn );
        nC   = size( ACQ.Data, dim.chan );        
        nZnT = nT * nZ; % = No. 'Frames' in Sweep volume
        
        % Configure Sweep k-space Kernels
        clear swpKernels ktAcqSwpKrn ktTrnSwpKrn
        
        swpKrnHalfWidth = ceil( swpKrnFullWidth / 2 );
        
        if swpKrnLoca
            swpKrnSpacing = unique( diff( swpKrnLoca ) );
        end
        
        if isempty(swpKrnSpacing)
            swpKrnSpacing = swpKrnFullWidth; % M2D equivalent
        end
        
        if isempty(swpKrnLoca)
            swpKrnLoca = swpKrnHalfWidth:swpKrnSpacing:nZnT;
        end
        
        nSwpKrn = numel( swpKrnLoca );
        
        % Array Sweep Kernels
        for iKrn = 1:nSwpKrn
            swpKernels(:,iKrn) = ...
                swpKrnLoca(iKrn)-swpKrnHalfWidth+1:swpKrnLoca(iKrn)+swpKrnHalfWidth;
        end
        
        % Ensure Kernels within Acquisition Bounds
        [~,swpKrnOutOfBounds,~] = find(swpKernels > nZnT);
        swpKernels( :, unique(swpKrnOutOfBounds) ) = [];
        nSwpKrn = size( swpKernels, 2 );

        % Plot Sweep Kernels
        if isVerbose

            % Sweep Kernels
            figure; hold on;
            plot(swpKernels','.b');
            
            % M2D Frames
            m2dKrnLoca = 0:nT:nZnT;
            for iKrn = 1:numel(m2dKrnLoca)
                plot(1:nSwpKrn, repmat(m2dKrnLoca(iKrn),1,nSwpKrn),'k--');
            end
            
            xlabel('Sweep Kernel No.'); ylabel('Frame Index');
            legend('Sweep Kernels','Location','NorthWest');
            axis([1 nSwpKrn 1 nZnT]);
            
            % Save
            hFig = gcf; hFig.Name = strcat( outFilePrefix, '_swp_kernels' );
            saveas( hFig, [outputDirPath '/' hFig.Name, '.fig' ] );
            saveas( hFig, [outputDirPath '/' hFig.Name, '.png' ] ); clear hFig;
            
        end
        
        
        % Reshape ACQ & TRN Data
        ktAcqSwp    = swap_dim_reconframe_to_xydcl( ACQ.Data );
        ktAcqSwp    = reshape( permute( ktAcqSwp, [1,2,4,3,5] ), nX, nY, nC, nT*nZ );
        ktAcqSwpKrn = single( zeros ( nX, nY, nC, nSwpKrn, swpKrnFullWidth ) );
        for iKrn = 1:nSwpKrn
            ktAcqSwpKrn(:,:,:,iKrn,:) = ktAcqSwp(:,:,:,swpKernels(:,iKrn));
        end
        ktAcqSwpKrn = permute( ktAcqSwpKrn, [1,2,5,3,4] );
        clear ktAcqSwp
        
        ktTrnSwp    = swap_dim_reconframe_to_xydcl( TRN.Data );
        ktTrnSwp    = reshape( permute( ktTrnSwp, [1,2,4,3,5] ), nX, nY, nC, nT*nZ );
        ktTrnSwpKrn = single( zeros ( nX, nY, nC, nSwpKrn, swpKrnFullWidth ) );
        for iKrn = 1:nSwpKrn
            ktTrnSwpKrn(:,:,:,iKrn,:) = ktTrnSwp(:,:,:,swpKernels(:,iKrn));
        end
        ktTrnSwpKrn = permute( ktTrnSwpKrn, [1,2,5,3,4] );
        clear ktTrnSwp
        
        
        % View Sweep Kernel Acq K-Space
        if isVerbose
            
            figure;
            if nSwpKrn > 12
                numKrnToPlot = 12;
            else
                numKrnToPlot = nSwpKrn;
            end
            krnToPlot = round(linspace(1,nSwpKrn,numKrnToPlot));
            for ii = 1:numKrnToPlot
                subplot(4,3,ii);
                imagesc(squeeze(abs( ...
                    ktAcqSwpKrn(ceil(nX/2),:,:,ceil(nC/2),krnToPlot(ii)) )), [0, 500] );
                title(['Kernel No. ' num2str(krnToPlot(ii))]);
                colormap('gray');
            end
            
            % Save
            hFig = gcf; hFig.Name = strcat( outFilePrefix, '_swp_kernel_kspace' );
            saveas( hFig, [outputDirPath '/' hFig.Name, '.fig' ] );
            saveas( hFig, [outputDirPath '/' hFig.Name, '.png' ] ); clear hFig;
             
        end
        
        fprintf( 'No. Sweep Kernels    = %3i \n', nSwpKrn ),
        fprintf( 'Sweep Kernel Width   = %3i \n', swpKrnFullWidth ),
        fprintf( 'Sweep Kernel Spacing = %3i \n', swpKrnSpacing ),

        % Save Sweep Kernel Parameters
        swpParamMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_swp_parameters.mat' ) );
        save( swpParamMatFilePath, 'swpKernels', 'nSwpKrn', 'swpKrnFullWidth', 'swpKrnHalfWidth', 'swpKrnLoca', 'swpKrnSpacing', '-v7.3' );   
        
        
        % Weighted CSM Message
        fprintf( 'Creating Weighted CSM ... \n' ),
        
        % CSM Original Locations
        csmLoca = 1:nZnT;
        csmLoca = reshape(csmLoca, [], nZ);

        % CSM Identifier Matrix
        csmIDs = 1:nZ;
        csmIDs = repmat(csmIDs,nT,1);

        
        % Calculate CSM Weighting Matrices
        csmOverlapFractions = [];
        uWgt = {}; % Note: need to be arrays to account for kernels overlapping singular or multiple csmLoca
        W    = {}; % Big matrix array of weights for simple .* with original CSM
        
        for iKrn = 1:nSwpKrn
            
            currKrn    = swpKernels(:,iKrn);
            currWgt    = csmIDs(csmLoca(currKrn));
            uWgt{iKrn} = unique(currWgt);
            
            % Overlap Between Kernels and Original CSM
            for iWgt = 1:numel(uWgt{iKrn})
                csmOverlapFractions(iWgt,iKrn) = sum(currWgt == uWgt{iKrn}(iWgt));
                Wgts(iWgt,iKrn)          = csmOverlapFractions(iWgt,iKrn) / swpKrnFullWidth;
                
                % Create Weights Matrices
                W{iKrn}(:,:,:,:,iWgt) = Wgts(iWgt,iKrn) * single( ones( [nX, nY, 1, nC] ) );
            end
            
        end
        
        % Create Weighted CSM by Linear Combination
        csmWgt = single( zeros( [nX, nY, 1, nC, nSwpKrn] ) );
        
        for iKrn = 1:nSwpKrn
            csmWgt(:,:,:,:,iKrn) = sum( W{iKrn} .* csm(:,:,:,:,uWgt{iKrn}), 5 );
        end
        
        % Save Weighted CSM
        csmWgtMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_csm_wgt.mat' ) );
        save( csmWgtMatFilePath, 'csmWgt', 'Wgts', 'uWgt', 'csmLoca', '-v7.3' );
                
        % Memory Management
        clear W Wgts uWgt csmLoca csmIDs csmOverlapFractions currKrn currWgt

        % TODO: Do I need to update ACQ/TRN objects here?:
        % Maybe I should perform the csmSwp processing in the csm block
        %     ACQ.Parameter.Recon.Sensitivities = SENS;
        %     TRN.Parameter.Recon.Sensitivities = SENS;
        
        numSlice = nSwpKrn;
        
        RCN.Data = single( zeros( nX,nY,1,1,swpKrnFullWidth,1,1,nSwpKrn ) );
        DC.Data  = single( zeros( nX,nY,1,1,1,1,1,nSwpKrn ) );
        
    end
    
    
    % Process Slice-by-Slice
%     keyboard; 
%     W=74; RCN.Data = RCN.Data(:,:,:,:,1:W,:,:,:); % vary M2D width
    for iSlice = 1:numSlice
        
        % Get Data from MRecon Objects
        if isSweepAcq
            ktAcq       = ktAcqSwpKrn(:,:,:,:,iSlice);
            ktTrn       = ktTrnSwpKrn(:,:,:,:,iSlice);
            csm         = csmWgt(:,:,:,:,iSlice);
%             warning('REMINDER TOMO: CSM original in use.');
%             csm         = swap_dim_reconframe_to_xydcl( ACQ.Parameter.Recon.Sensitivities.Sensitivity(:,:,:,:,:,:,:,iSlice,:,:,:,:) );
            noiseCov    = SENS.Psi;
            
            % Display Start Message
            disp_start_step_msg( sprintf( '  kernel %3i of %1i  \n', iSlice, numSlice ) )
            
        else
            ktAcq       = swap_dim_reconframe_to_xydcl( ACQ.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:) );
            ktTrn       = swap_dim_reconframe_to_xydcl( TRN.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:) );
            csm         = swap_dim_reconframe_to_xydcl( ACQ.Parameter.Recon.Sensitivities.Sensitivity(:,:,:,:,:,:,:,iSlice,:,:,:,:) );
            noiseCov    = SENS.Psi;
            
%             warning('Hard-code REMINDER TOMO: Reduced nFrames in k-t recon');
%             ktAcq = ktAcq(:,:,75-(W/2):74+(W/2),:);
%             ktTrn = ktTrn(:,:,75-(W/2):74+(W/2),:);
            
            % Display Start Message
            disp_start_step_msg( sprintf( '  slice %3i of %1i  \n', iSlice, numSlice ) )

        end
        
        % Get Custom Training Data
        if ~isempty( xtTrnCus )
            xtTrnCusSlice = squeeze(xtTrnCus(:,:,:,:,:,:,:,iSlice));
        end
        
        % k-t SENSE Reconstruction
        if isempty( mask )
            % Uniform Regularization
            [ ktRcn, ktDc, ~, ~, xfMask, xfPri, ~, ~, xtRcn ] = recon_ktsense( ktAcq, ktTrn, csm, 'noisecov', noiseCov, 'lambda0', ktRegStrength );
        else
            if isempty( xtTrnCus )
                % Adaptive Regularization
                [ ktRcn, ktDc, ~, ~, xfMask, xfPri, ~, xtTrnCusSlice ] = recon_ktsense( ktAcq, ktTrn, csm, 'noisecov', noiseCov, 'lambda0', ktRegStrength,  'mask', mask(:,:,iSlice), 'lambdaroi', ktRegStrengthROI, 'makeHarmonicFilter', makeHarmonicFilter, 'frameDuration', stackFrameDuration );
            elseif ~isempty( xtTrnCus )
                % Adaptive Regularization with Custom Training Data
                [ ktRcn, ktDc, ~, ~, xfMask, xfPri, ~, xtTrnCusSlice ] = recon_ktsense( ktAcq, ktTrn, csm, 'noisecov', noiseCov, 'lambda0', ktRegStrength,  'mask', mask(:,:,iSlice), 'lambdaroi', ktRegStrengthROI, 'removeoversampling', true, 'xtTrnCus', xtTrnCusSlice, 'makeHarmonicFilter', makeHarmonicFilter, 'frameDuration', stackFrameDuration );
            end
        end
        % TODO: compare recon using available noise covariance estimates:
        %   1) estimated from k-t undersampled data,
        %   2) calculated using noise samples in MRecon object NOISE, and
        %   3) calculated in sensitivity maps MRecon object SENS
        
        % Put Data in Back in MRecon Objects
        RCN.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:) = single( swap_dim_xydcl_to_reconframe( ktRcn ) );
        DC.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:)  = single( swap_dim_xydcl_to_reconframe( ktDc ) );
        
        % Save x-f Mask and Training Data
        if isSweepAcq
            xfmaskMatFilePath = fullfile( outputDirPath, sprintf( '%s_slice%02i_pri_swp.mat', outFilePrefix, iSlice ) );
            save( xfmaskMatFilePath, 'xfMask', 'xfPri', '-v7.3' );
            clear xfMask xfPri
        end
        
        if makeHarmonicFilter
            xfmaskMatFilePath = fullfile( outputDirPath, sprintf( '%s_slice%02i_pri_hrm.mat', outFilePrefix, iSlice ) );
            save( xfmaskMatFilePath, 'xfMask', 'xfPri', 'xtTrnCusSlice', '-v7.3' );
            clear xfMask xfPri xtTrn xtTrnCusSlice
        end
        
        % Save Reconstructed Data
        reconMatFilePath = fullfile( outputDirPath, sprintf( '%s_slice%02i_recon.mat', outFilePrefix, iSlice ) );
        save( reconMatFilePath, 'ktRcn', 'ktDc', '-v7.3' );
        clear ktRcn ktDc
        
        % Display Time Elapsed Message
        disp_time_elapsed_msg( toc ),
        disp_write_file_msg( reconMatFilePath ),
        
    end
    
    % Memory Management
    clear xfMask xfPri
    clear ktAcqSwpKrn ktTrnSwpRcn
    clear ktAcqSwp ktTrnSwp
    clear ktAcq ktTrn
    clear csmSwp imCoilSwp imBodySwp
    
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
    rltAbNiiFilePath   = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_ab' ) );
    rltAbNiiFilePath   = mrecon_writenifti( RCN, rltAbNiiFilePath, 'frameduration', frameDuration, 'datatype', 'magnitude' );
    rltReNiiFilePath   = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_re' ) );
    rltReNiiFilePath   = mrecon_writenifti( RCN, rltReNiiFilePath, 'frameduration', frameDuration, 'datatype', 'real' );
    rltImNiiFilePath   = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_im' ) );
    rltImNiiFilePath   = mrecon_writenifti( RCN, rltImNiiFilePath, 'frameduration', frameDuration, 'datatype', 'imaginary' );
    rltPhNiiFilePath   = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_ph' ) );
    rltPhNiiFilePath   = mrecon_writenifti( RCN, rltPhNiiFilePath, 'frameduration', frameDuration, 'datatype', 'phase' );
    rltCplxNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_cplx' ) );
    rltCplxNiiFilePath = mrecon_writenifti( RCN, rltCplxNiiFilePath, 'frameduration', frameDuration, 'datatype', 'complex' );
   
    
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
    if isSweepAcq
        PARAM.Timing = mrecon_getslicetiming_sweep( RCN, 'patchversion', patchVersion );
        % TODO?: 
        % Update NrDyn for compatibility with preproc.m?
        % PARAM.Encoding.NrDyn = [swpKrnFullWidth, swpKrnFullWidth];
    else
        PARAM.Timing = mrecon_getslicetiming( RCN, 'patchversion', patchVersion );
    end
    rltParamMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_parameters.mat' ) );
    save( rltParamMatFilePath, 'PARAM', '-v7.3' );
    
    % Display Time Elapsed Message
    disp_time_elapsed_msg( toc ),
    disp_write_file_msg( rltAbNiiFilePath )
    disp_write_file_msg( rltReNiiFilePath )
    disp_write_file_msg( rltImNiiFilePath )
    disp_write_file_msg( rltPhNiiFilePath )
    disp_write_file_msg( rltCplxNiiFilePath )
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

% Note: not available in MRecon 515

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
