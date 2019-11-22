function mrecon_kt_hrh( rawDataFilePath, varargin )
%MRECON_KT_HRH  Hierarchical k-t SENSE processing using ReconFrame
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

% tar   (t.roberts@kcl.ac.uk)
% jfpva (joshua.vanamerom@kcl.ac.uk)


%% Optional Input Argument Default Values

default.senseRefFilePath    = '';
default.coilSurveyFilePath  = '';
default.outputDirPath       = pwd;
default.outFilePrefix       = '';
default.reconOpts           = {};
default.patchVersion        = 'PIH1';
default.isExportPdf         = false;
default.isExportJson        = false;
default.isExportGoalc       = false;
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


%% Calculate Coil Sensitivity Maps

% TODO: Replace with Lucilio's Sensitivities

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

%         % Extract data from MRsense object
%         csm     = swap_dim_reconframe_to_xydcl( SENS.Sensitivity );             % coil sensitivity maps
%         imBody  = swap_dim_reconframe_to_xydcl( SENS.ReformatedBodycoilData );  % body coil image
%         imCoil  = swap_dim_reconframe_to_xydcl( SENS.ReformatedCoilData );      % array coil images
%         psiSens = SENS.Psi;                                                     % array coil noise covariance
%         csmMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_csm.mat' ) );
%         save( csmMatFilePath, 'csm', 'imBody', 'imCoil', 'psiSens', '-v7.3' );
%         
%         disp_time_elapsed_msg( toc )
%         
%         disp_write_file_msg( csmMatFilePath )
        
    otherwise

        error( 'Unrecognised CSM calculation method: %s', csmCalcMethod ),
        
end

% Add Sensitivity Maps to Undersampled MRecon Object
ACQ.Parameter.Recon.Sensitivities = SENS;


%% Hierarchical k-t SENSE Processing
% - Required to ensure data consistent between Lucilio/Josh k-t frameworks

% Display Start Message
disp_start_step_msg( 'Making Hierarchical k-t SENSE data consistent' ),

% Hierarchical k-t SENSE MRecon Objects
HRH = ACQ.Copy;
HRH.Data = sum (HRH.Data, dim.chan );

DC = ACQ.Copy;
DC.Data = sum( sum( DC.Data, dim.chan ), dim.dyn );

% Load Hierarchical k-Space Data
disp_start_step_msg( 'Loading Hierarchical k-Space Data' ),
hrhKSpaceFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_kspace_hrh.mat' ) );
load( hrhKSpaceFilePath, 'y', 'A' );
disp_time_elapsed_msg( toc ), 

% Uncompress
NY=size(y);NY(end+1:5)=1; %TODO: replace 5 with size()
NA=size(A);NA(end+1:5)=1;
NR=NY./NA;NR(2)=1;
A=repmat(A,NR);
z=A;z(:)=0;
z(A==1)=y;
y=z;z=[];A=[];

% Put k-Space Data Back in MRecon Object
ktHrh   = permute( y, [1 2 5 4 3] );
ktHrhDc = mean( ktHrh, 3 );
clear y A
HRH.Data = swap_dim_xydcl_to_reconframe( ktHrh );
DC.Data  = swap_dim_xydcl_to_reconframe( ktHrhDc );

% Perform Dummy k2i to Update MRecon Object 
% - i.e.: mrecon_k2i( RCN )
mrecon_k2i( HRH );
mrecon_k2i( DC );

% Perform Dummy Postprocessing to Update MRecon Object 
% - i.e.: partial mrecon_postprocess( RCN );
HRH.SENSEUnfold;
HRH.PartialFourier;
HRH.ConcomitantFieldCorrection;
HRH.DivideFlowSegments;
HRH.CombineCoils;

DC.SENSEUnfold;
DC.PartialFourier;
DC.ConcomitantFieldCorrection;
DC.DivideFlowSegments;
DC.CombineCoils;

% Replace MRecon Object Data with Hrh Image Data
disp_start_step_msg( 'Loading Hierarchical x-t Data' ),
hrhReconMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_hrh_recon.mat' ) );
load( hrhReconMatFilePath );
disp_time_elapsed_msg( toc ),

HRH.Data = xtRcn;
DC.Data  = mean( xtRcn , dim.dyn );

% Perform Image Postprocessing Missing from Lucilio Framework
HRH.Average;
HRH.GeometryCorrection;
HRH.RemoveOversampling;
HRH.FlowPhaseCorrection;
HRH.ReconTKE;
HRH.ZeroFill;
HRH.RotateImage;

DC.Average;
DC.GeometryCorrection;
DC.RemoveOversampling;
DC.FlowPhaseCorrection;
DC.ReconTKE;
DC.ZeroFill;
DC.RotateImage;

% Get Timing
frameDuration = mrecon_calc_frame_duration( HRH );

% Create Postprocessed xt variables
xtRcn = HRH.Data;
xtDc  = DC.Data;

% Write to NIfTI
hrhAbNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_ab' ) );
hrhAbNiiFilePath = mrecon_writenifti( HRH, hrhAbNiiFilePath, 'frameduration', frameDuration, 'datatype', 'magnitude' );
hrhReNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_re' ) );
hrhReNiiFilePath = mrecon_writenifti( HRH, hrhReNiiFilePath, 'frameduration', frameDuration, 'datatype', 'real' );
hrhImNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_im' ) );
hrhImNiiFilePath = mrecon_writenifti( HRH, hrhImNiiFilePath, 'frameduration', frameDuration, 'datatype', 'imaginary' );

hrhDcAbNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_ab' ) );
hrhDcAbNiiFilePath = mrecon_writenifti( DC, hrhDcAbNiiFilePath, 'datatype', 'magnitude' );
hrhDcReNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_re' ) );
hrhDcReNiiFilePath = mrecon_writenifti( DC, hrhDcReNiiFilePath, 'datatype', 'real' );
hrhDcImNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_im' ) );
hrhDcImNiiFilePath = mrecon_writenifti( DC, hrhDcImNiiFilePath, 'datatype', 'imaginary' );

% Save as .mat
hrhMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_recon.mat' ) );
save( hrhMatFilePath, 'xtRcn', '-v7.3' );

hrhDcMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_recon.mat' ) );
save( hrhDcMatFilePath, 'xtDc', '-v7.3' );

clear xtHrh xtHrhDc xtRcn

% Display Time Elapsed Message
disp_time_elapsed_msg( toc )
disp_write_file_msg( hrhAbNiiFilePath )
disp_write_file_msg( hrhReNiiFilePath )
disp_write_file_msg( hrhImNiiFilePath )
disp_write_file_msg( hrhMatFilePath )
disp_write_file_msg( hrhDcAbNiiFilePath )
disp_write_file_msg( hrhDcReNiiFilePath )
disp_write_file_msg( hrhDcImNiiFilePath )
disp_write_file_msg( hrhDcMatFilePath )


%%  Save Parameters
PARAM = mrecon_getparameters( HRH );
PARAM.Timing = mrecon_getslicetiming( HRH, 'patchversion', patchVersion );
hrhParamMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_parameters.mat' ) );
save( hrhParamMatFilePath, 'PARAM', '-v7.3' );

% Display Time Elapsed Message
disp_time_elapsed_msg( toc ),
disp_write_file_msg( hrhParamMatFilePath )


%% Perform Adaptive Regularisation k-t SENSE Reconstruction with HRH Images as Priors

% TODO: think I want to implement here.
% Essentially, a call to recon_ktsense.m


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
    HRH.Parameter.Export2Json( jsonFilePath );
    
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


end  % mrecon_kt_hrh(...)
