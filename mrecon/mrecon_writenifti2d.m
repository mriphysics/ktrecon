function niiFilePaths = mrecon_writenifti2d( MR, niiFileDir, niiFileNamePrefix, varargin )
%MRECON_WRITENIFTI2D  write MRecon object image-space data to NIfTI as 2D+time
%
%   MRECON_WRITENIFTI2D( MR, niiFileDir, niiFileNamePrefix )
%
%   niiFilePaths = MRECON_WRITENIFTI2D( MR, niiFilePath, niiFileNamePrefix )
%
%   MRECON_WRITENIFTI2D( ... , 'param', val) 
%
%   Output:
%       niiFilePaths        - cell array of paths to saved NIfTI files
%   
%   Input:
%       MR                  - ReconFrame MRecon object
%       niiFileDir          - path to directory to write to
%       niiFileNamePrefix   - NIfTI file name prefix
%
%   Parameter-Value Options:
%       frameduration       - duration, in seconds, of one dynamic frame
%       complex             - save complex-valued images
%       displayrange        - NIfTI header [min,max] values
%       datascaling         - [interept,slope], 
%                             i.e., im = slope * MR.Data + intercept
%       patchVersion        - version string of patch used to acquire data,
%                             used to determine slice timing in m2d stack

% JFPvA (joshua.vanamerom@kcl.ac.uk)


%% Default Values for Unspecified Input Arguments

default.frameDuration  = [];
default.saveComplex    = false;
default.displayRange   = [];
default.dataScaling    = [0 1];
default.patchVersion   = '';


%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(    p, 'MR', ...
    @(x) validateattributes( x, {'MRecon'}, {'scalar'}, mfilename) );
addRequired(   p, 'niiFileDir', ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );
addRequired(   p, 'niiFileNamePrefix', ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );
add_param_fn(   p, 'frameduration', default.frameDuration, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','real','scalar'}, mfilename) );    
add_param_fn(   p, 'complex', default.saveComplex, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );
add_param_fn(   p, 'displayrange', default.displayRange, ...
    @(x) validateattributes( x, {'numeric'}, {'vector','numel',2}, mfilename) );
add_param_fn(   p, 'datascaling', default.dataScaling, ...
    @(x) validateattributes( x, {'numeric'}, {'vector','numel',2}, mfilename) );
add_param_fn(   p, 'patchversion', default.patchVersion, ...
    @(x) validateattributes( x, {'char'}, {}, mfilename) );

parse( p, MR, niiFileDir, niiFileNamePrefix, varargin{:} );

frameDuration  = p.Results.frameduration;
saveComplex    = p.Results.complex;
displayRange   = p.Results.displayrange;
dataScaling    = p.Results.datascaling;
patchVersion   = p.Results.patchversion;

if saveComplex && isreal( MR.Data )
    saveComplex = false;
end


%% Anonymous Function

% mrecon: x-y-z-chan-dyn-card-echo-loca-mix-extr1-extr2-aver
    %         1-2-3-4----5---6----7----8----9---10----11----12--
    %
    % nii:    x-y-loca-dyn-other
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
reshape_mrecon_to_nii = @( d ) reshape( permute( d, [ dim.x dim.y dim.loca dim.dyn dim.z dim.chan dim.card dim.echo dim.mix dim.extr1 dim.extr2 dim.aver ] ), size(d,dim.x), size(d,dim.y), size(d,dim.loca), size(d,dim.dyn), [] ); 


%% Calculate Timing from MRecon Object

% Get Study Start Time
tStudyStr = MR.Parameter.GetValue('RFR_ECSERIES_DICOM_SERIES_TIME');
ind = strfind( tStudyStr, '.' );
hr  = str2double( tStudyStr(1:(ind-5)) );
min = str2double( tStudyStr((ind-4):(ind-3)) );
sec = str2double( tStudyStr((ind-2):end) );
tStudy = 3600*hr + 60*min + sec;

% Get Series Start Time
tSeriesStr = MR.Parameter.GetValue('RFR_SERIES_DICOM_SERIES_TIME');
ind = strfind( tSeriesStr, '.' );
hr  = str2double( tSeriesStr(1:(ind-5)) );
min = str2double( tSeriesStr((ind-4):(ind-3)) );
sec = str2double( tSeriesStr((ind-2):end) );
tSeries = 3600*hr + 60*min + sec;

%{
% NOTE: datetime isn't available in earlier Matlab versions,
%       otherwise, could do something like the following
timestr2datetime = @(str) datetime( str, 'InputFormat', 'HHmmss.SSSSS' );
datetime2seconds = @(t) 360*t.Hour + 60*t.Minute + t.Second;
tStudy  = datetime2seconds( timestr2datetime( MR.Parameter.GetValue('RFR_ECSERIES_DICOM_SERIES_TIME') ) );
tSeries = datetime2seconds( timestr2datetime( MR.Parameter.GetValue('RFR_SERIES_DICOM_SERIES_TIME') ) );
%}

% Get Frame Duration
if isempty( frameDuration )
    frameDuration = mrecon_calc_frame_duration( MR );
end

% Get Slice Duration
seriesDuration  = MR.Parameter.GetValue('AC_total_scan_time');
numLoc        = size( MR.Data, dim.loca );
numDynDummy     = MR.Parameter.GetValue( 'MP_nr_dummy_dynamic_scans' );
switch patchVersion
    case 'c51c8c1'
        % acquition order: dummy, acq_1, dummy, acq_2, ..., dummy, acq_n
        %                  dummy, trn_1, dummy, trn_2, ..., dummy, trn_n
        sliceDuration    = seriesDuration / numLoc;
        sliceStartOffset = numDynDummy * frameDuration;
    otherwise
        fprintf( '    Warning: patch version (%s) unspecified or unrecognised; slice timing may be incorrect.  \n', patchVersion )
        sliceDuration    = seriesDuration / numLoc;
        sliceStartOffset = numDynDummy * frameDuration;
end

fprintf( '    %-25s %-15s %15.3f\n', 'study start time', '(hhmmss.sss)', str2num( tStudyStr ) );
fprintf( '    %-25s %-15s %15.3f\n', 'series start time', '(hhmmss.sss)', str2num( tSeriesStr ) );
fprintf( '    %-25s %-15s %15.3f\n', 'study start time', '(s)', tStudy );
fprintf( '    %-25s %-15s %15.3f\n', 'series start time', '(s)', tSeries );
fprintf( '    %-25s %-15s %15.3f\n', 'series time offset', '(s)', tSeries - tStudy );
fprintf( '    %-25s %-15s %15.3f\n', 'series duration', '(s)', seriesDuration );
fprintf( '    %-25s %-15s %15.3f\n', 'slice durationn', '(s)', sliceDuration );
fprintf( '    %-25s %-15s %15.3f\n', 'dummy delay', '(s)', sliceStartOffset );
fprintf( '    %-25s %-15s %15.3f\n', 'frame duration', '(s)', frameDuration );


%% Get Image Information

% Compute M2D Image to World Transformation
A  = mrecon_image_to_nii_world_transform( MR );

% Calculate Scaling of m2d Data (Voxel Size)
Sm2d  = sqrt( sum( A(1:3,1:3).^2, 1 ) );

% Get Slice Thickness
sliceThickness = MR.Parameter.Scan.RecVoxelSize(3);

% Scaling of 2d Data (Voxel Size)
S2d = [ Sm2d(1) Sm2d(2) sliceThickness ];


%% Extract 2D+Time Slices from M2D Stack

niiFilePaths = cell(numLoc,1);

locOrder = mrecon_locorder( MR ) + 1;  % NOTE: change from zero index to index 1 

fprintf( '    slice order: ' );
fprintf( '%i, ', locOrder-1 );
fprintf( '\b\b  \n' );

for iLoc = locOrder(:).'  % NOTE: for loop index values should be row vector

    % Full Path of NIfTI File
    niiFilePaths{iLoc} = fullfile( niiFileDir, sprintf( '%s_slice%02i.nii', niiFileNamePrefix, iLoc ) );

    % Format Image Data
    im = reshape_mrecon_to_nii( MR.Data(:,:,:,:,:,:,:,iLoc,:,:,:,:) );
    if ~saveComplex
        im = abs( im );
    end
    im = dataScaling(2) * im + dataScaling(1);

    % Slice Time Offset
    sliceTimeOffset = ( tSeries + double(iLoc-1) * sliceDuration + sliceStartOffset - tStudy );  % seconds
    fprintf( '    %-25s %-15s %15.3f\n', sprintf( 'slice%02i time offset', iLoc ), '(s)', sliceTimeOffset );
    
    % Calculate Affine Transformation for Slice
    Ashift2slice = eye(4); Ashift2slice(3,4) = iLoc-1;
    Aunscalevox = inv( [Sm2d(1) 0 0 0; 0 Sm2d(2) 0 0; 0 0 Sm2d(3) 0; 0 0 0 1] );
    Arescalevox = [S2d(1) 0 0 0; 0 S2d(2) 0 0; 0 0 S2d(3) 0; 0 0 0 1];
    Aslice = A * Arescalevox * Aunscalevox * Ashift2slice;
    
    % Make NIfTI Stucture
    N = make_nii( im, S2d );

    % Adjust NIfTI Header Values
    N.hdr.dime.pidxim(4) = sliceThickness;
    N.hdr.dime.pixdim(5) = frameDuration;
    N.hdr.dime.toffset   = sliceTimeOffset;
    N.hdr.hist.srow_x = Aslice(1,:);
    N.hdr.hist.srow_y = Aslice(2,:);
    N.hdr.hist.srow_z = Aslice(3,:);
    N.hdr.hist.sform_code = 1;
    if ~isempty( dataScaling )
        N.hdr.dime.scl_inter = dataScaling(1);
        N.hdr.dime.scl_slope = dataScaling(2);
    end
    if ~isempty( displayRange )
        N.hdr.dime.cal_min  = displayRange(1);
        N.hdr.dime.cal_max  = displayRange(2);
        N.hdr.dime.glmin    = displayRange(1);
        N.hdr.dime.glmax    = displayRange(2);
    end

    % Save NIfTI
    save_nii( N, niiFilePaths{iLoc} );
    system( sprintf( 'gzip -f %s', niiFilePaths{iLoc} ) );
    niiFilePaths{iLoc} = strcat( niiFilePaths{iLoc}, '.gz' );

end  


end  % mrecon_writenifti2d(...)
