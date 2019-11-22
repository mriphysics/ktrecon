function [  ] = kt_sc_slw_preproc( outputDirPath, outFilePrefix, SCSLW )
% KT_SC_SLW_PREPROC         k-t self-calibrated sliding window preprocessing
%
%   [  ] = kt_sc_slw_preproc( ACQ, outFilePrefix )
%
%   performs missing MRecon steps to make Lucilio and Josh reconstructions
%   consistent, such as image geometry correction
%
%   output:
%       Various matfiles
%
%   input:
%       outDirPath                  output directory path
%       outFilePrefix               output file name prefix
%       SCSLW                       MRecon object
%
%

%   tar (t.roberts@kcl.ac.uk)


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


%% Display Start Message
disp_start_step_msg( 'Performing k-t self-calibrated sliding window preprocessing' ),


%% MRecon Objects
SCSLW.Data = sum ( SCSLW.Data, dim.chan );

DC = SCSLW.Copy;
DC.Data = sum( sum( DC.Data, dim.chan ), dim.dyn );


%% Load k-Space Data
disp_start_step_msg( 'Loading Self-calibrated Sliding Window k-Space Data' ),
scSlwMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_kspace_sc_slw.mat' ) );
load( scSlwMatFilePath, 'y', 'A' );
disp_time_elapsed_msg( toc ), 


%% Uncompress Data
NY=size(y);NY(end+1:5)=1;
NA=size(A);NA(end+1:5)=1;
NR=NY./NA;NR(2)=1;
A=repmat(A,NR);
z=A;z(:)=0;
z(A==1)=y;
y=z;z=[];A=[];


%% Put k-Space Data Back in MRecon Object
ktScSlw   = permute( y, [1 2 5 4 3] );
ktScSlwDc = mean( ktScSlw, 3 );
clear y A
SCSLW.Data = swap_dim_xydcl_to_reconframe( ktScSlw );
DC.Data  = swap_dim_xydcl_to_reconframe( ktScSlwDc );


%% Perform Dummy k2i to Update MRecon Object 
mrecon_k2i( SCSLW );
mrecon_k2i( DC );


%% Perform Dummy Postprocessing to Update MRecon Object 
% - i.e.: partial mrecon_postprocess( RCN );
SCSLW.SENSEUnfold;
SCSLW.PartialFourier;
SCSLW.ConcomitantFieldCorrection;
SCSLW.DivideFlowSegments;
SCSLW.CombineCoils;

DC.SENSEUnfold;
DC.PartialFourier;
DC.ConcomitantFieldCorrection;
DC.DivideFlowSegments;
DC.CombineCoils;


%% Replace MRecon Object Data with ScSlw Image Data
disp_start_step_msg( 'Loading Self-calibrated Sliding Window x-t Data' ),
scSlwReconMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_sc_slw_recon.mat' ) );
load( scSlwReconMatFilePath );
disp_time_elapsed_msg( toc ),
    
SCSLW.Data = xtRcn;
DC.Data    = mean( xtRcn , dim.dyn );


%% Perform Image Postprocessing Missing from Lucilio Framework

% Preprocessing check
if isPreprocessed
    fprintf('\nData has already been preprocessed. Exiting function ...');
    return;    
end

SCSLW.Average;
SCSLW.GeometryCorrection;
% SCSLW.RemoveOversampling;
SCSLW.FlowPhaseCorrection;
SCSLW.ReconTKE;
SCSLW.ZeroFill;
SCSLW.RotateImage;

DC.Average;
DC.GeometryCorrection;
% DC.RemoveOversampling;
DC.FlowPhaseCorrection;
DC.ReconTKE;
DC.ZeroFill;
DC.RotateImage;

% Get Timing
frameDuration = mrecon_calc_frame_duration( SCSLW );

% Create Postprocessed xt variables
xtRcn = SCSLW.Data;
xtDc  = DC.Data;

% Change isPreprocessed flag
isPreprocessed = true;


%% Write to NIfTI
scSlwAbNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_sc_slw_ab' ) );
scSlwAbNiiFilePath = mrecon_writenifti( SCSLW, scSlwAbNiiFilePath, 'frameduration', frameDuration, 'datatype', 'magnitude' );
scSlwReNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_sc_slw_re' ) );
scSlwReNiiFilePath = mrecon_writenifti( SCSLW, scSlwReNiiFilePath, 'frameduration', frameDuration, 'datatype', 'real' );
scSlwImNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_sc_slw_im' ) );
scSlwImNiiFilePath = mrecon_writenifti( SCSLW, scSlwImNiiFilePath, 'frameduration', frameDuration, 'datatype', 'imaginary' );

scSlwDcAbNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_sc_slw_ab' ) );
scSlwDcAbNiiFilePath = mrecon_writenifti( DC, scSlwDcAbNiiFilePath, 'datatype', 'magnitude' );
scSlwDcReNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_sc_slw_re' ) );
scSlwDcReNiiFilePath = mrecon_writenifti( DC, scSlwDcReNiiFilePath, 'datatype', 'real' );
scSlwDcImNiiFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_sc_slw_im' ) );
scSlwDcImNiiFilePath = mrecon_writenifti( DC, scSlwDcImNiiFilePath, 'datatype', 'imaginary' );


%% Save as .mat
scSlwMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_sc_slw_recon.mat' ) );
save( scSlwMatFilePath, 'xtRcn', 'isPreprocessed', '-v7.3' );

scSlwDcMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_dc_sc_slw_recon.mat' ) );
save( scSlwDcMatFilePath, 'xtDc', 'isPreprocessed', '-v7.3' );

clear xtRcn


%% Display Time Elapsed Message
disp_time_elapsed_msg( toc )
disp_write_file_msg( scSlwAbNiiFilePath )
disp_write_file_msg( scSlwReNiiFilePath )
disp_write_file_msg( scSlwImNiiFilePath )
disp_write_file_msg( scSlwMatFilePath )
disp_write_file_msg( scSlwDcAbNiiFilePath )
disp_write_file_msg( scSlwDcReNiiFilePath )
disp_write_file_msg( scSlwDcImNiiFilePath )
disp_write_file_msg( scSlwDcMatFilePath )


% TODO: remove below
% %%  Save Parameters File
% PARAM = mrecon_getparameters( SCSLW );
% PARAM.Timing = mrecon_getslicetiming( SCSLW, 'patchversion', 'PIH1' ); %NB: patchversion hard-coded
% scSlwParamMatFilePath = fullfile( outputDirPath, strcat( outFilePrefix, '_rlt_parameters.mat' ) );
% save( scSlwParamMatFilePath, 'PARAM', '-v7.3' );
% 
% % Display Time Elapsed Message
% disp_time_elapsed_msg( toc ),
% disp_write_file_msg( scSlwParamMatFilePath )


end  % kt_sc_slw_preproc(...)