function T = mrecon_getslicetiming_sweep( MR, varargin )
%MRECON_GETSLICETIMING  get timing of Sweep slices from MRecon object
%
%   MRECON_GETSLICETIMING_SWEEP( MR )
%
%   T = MRECON_GETSLICETIMING_SWEEP( MR )
%
%   MRECON_GETSLICETIMING_SWEEP( ... , 'param', val) 
%
%   Output:
%       T                   - struct of timing parameter values
%   
%   Input:
%       MR                  - ReconFrame MRecon object
%
%   Parameter-Value Options:
%       patchVersion        - version string of patch used to acquire data,
%                             used to determine slice timing in sweep stack
%       verbose             - if true, display timing information

% TAR   (t.roberts@kcl.ac.uk)
% JFPvA (joshua.vanamerom@kcl.ac.uk)


%% TODO:

% Calculate timings based on arbitrary sweep kernel locations and widths
% - currently only works for full acquisition volume



%% Default Values for Unspecified Input Arguments

default.patchVersion   = '';
default.isVerbose      = false;


%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(    p, 'MR', ...
    @(x) validateattributes( x, {'MRecon'}, {'scalar'}, mfilename) ); 
add_param_fn(   p, 'patchversion', default.patchVersion, ...
    @(x) validateattributes( x, {'char'}, {}, mfilename) );
add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, MR, varargin{:} );

patchVersion    = p.Results.patchversion;
isVerbose       = p.Results.verbose;


%% Dimensions

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
frameDuration = mrecon_calc_frame_duration( MR ); % NOTE: same for Sweep

% Get Slice Order
numLoc   = size( MR.Data, dim.loca );
locOrder = 1:numLoc; % NOTE: linear for Sweep, not expected to be segmented or other non-linear order
if ( isVerbose )
    fprintf( '    Sweep slice order: ' );
    fprintf( '%i, ', locOrder );
    fprintf( '\b\b  \n' );
end 

% Get Series Duration
% Note: AC_total_scan_time wrong if kt training data interleaved
if strcmp( MR.Parameter.GetValue('EX_ACQ_kt_training_interleaved') , 'yes' )
    tScanTimeStr = MR.Parameter.GetValue('IF_str_total_scan_time');
    ind = strfind( tScanTimeStr, ':' );
    min = str2double( tScanTimeStr( 1:2 ) );
    sec = str2double( tScanTimeStr((ind+1):end) );
    seriesDuration = 60*min + sec;
else
    seriesDuration = MR.Parameter.GetValue('AC_total_scan_time'); % TODO: possibly replace with IF_str_total_scan_time because no bugs with ktTrn interleave
end
   

% Get Slice Duration
numDynDummy     = MR.Parameter.GetValue( 'MP_nr_dummy_dynamic_scans' );
switch patchVersion
    case {'PIH1', 'PIH2'}
        if ( isVerbose )
            fprintf( '    patch version:    %s  \n', patchVersion )
            fprintf( '    aquisition order: dummy, acq_1, acq_2, ..., acq_n, \n' )
            fprintf( '                      dummy, trn_1, trn_2, ..., trn_n  \n' )
        end
        sliceDuration    = seriesDuration / numLoc;
        sliceStartOffset = numDynDummy * frameDuration;
    otherwise
        error( 'patch version (%s) not recognised.  \n', patchVersion )
end

if ( isVerbose )
    fprintf( '    %-25s %-15s %15.3f\n', 'study start time', '(hhmmss.sss)', str2double( tStudyStr ) );
    fprintf( '    %-25s %-15s %15.3f\n', 'series start time', '(hhmmss.sss)', str2double( tSeriesStr ) );
    fprintf( '    %-25s %-15s %15.3f\n', 'study start time', '(s)', tStudy );
    fprintf( '    %-25s %-15s %15.3f\n', 'series start time', '(s)', tSeries );
    fprintf( '    %-25s %-15s %15.3f\n', 'series time offset', '(s)', tSeries - tStudy );
    fprintf( '    %-25s %-15s %15.3f\n', 'series duration', '(s)', seriesDuration );
    fprintf( '    %-25s %-15s %15.3f\n', 'slice durationn', '(s)', sliceDuration );
    fprintf( '    %-25s %-15s %15.3f\n', 'dummy delay', '(s)', sliceStartOffset );
    fprintf( '    %-25s %-15s %15.3f\n', 'frame duration', '(s)', frameDuration );
end

% Calculate Slice Timing
sliceTimeOffset = nan(1,numLoc);
for iLoc = locOrder(:).'  % NOTE: for loop index values should be row vector
    % Slice Time Offset
    sliceTimeOffset(iLoc) = ( tSeries + double(iLoc-1) * sliceDuration + sliceStartOffset - tStudy );  % seconds
    if ( isVerbose )
        fprintf( '    %-25s %-15s %15.3f\n', sprintf( 'slice%02i time offset', iLoc ), '(s)', sliceTimeOffset(iLoc) );
    end
end


%% Save Timing Parameters to Struct

T.tStudyStr         = tStudyStr;
T.tSeriesStr        = tSeriesStr;
T.tStudy            = tStudy;
T.tSeries           = tSeries;
T.tSeriesOffset     = tSeries - tStudy;
T.seriesDuration    = seriesDuration;
T.numLoc            = numLoc;
T.numDynDummy       = numDynDummy;
T.sliceDuration     = sliceDuration;
T.sliceStartOffset  = sliceStartOffset;
T.frameDuration     = frameDuration;
T.sliceOrder        = locOrder(:)';
T.sliceTime         = sliceTimeOffset(:)';
T.sliceTimeOffset   = sliceTimeOffset(:)' - tSeries;


end  % mrecon_getslicetiming(...)
