function T = mrecon_getslicetiming( MR, varargin )
%MRECON_GETSLICETIMING  get timing of m2d slices from MRecon object
%
%   MRECON_GETSLICETIMING( MR )
%
%   T = MRECON_GETSLICETIMING( MR )
%
%   MRECON_WRITENIFTI2D( ... , 'param', val) 
%
%   Output:
%       T                   - struct of timing parameter values
%   
%   Input:
%       MR                  - ReconFrame MRecon object
%
%   Parameter-Value Options:
%       patchVersion        - version string of patch used to acquire data,
%                             used to determine slice timing in m2d stack
%       verbose             - if true, display timing information

% JFPvA (joshua.vanamerom@kcl.ac.uk)


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
frameDuration = mrecon_calc_frame_duration( MR );

% Get Slice Order
locOrder = mrecon_locorder( MR ) + 1;  % NOTE: change from zero index to index 1 
if ( isVerbose )
    fprintf( '    slice order: ' );
    fprintf( '%i, ', locOrder );
    fprintf( '\b\b  \n' );
end 

% Get Slice Duration
seriesDuration  = MR.Parameter.GetValue('AC_total_scan_time');
numLoc          = size( MR.Data, dim.loca );
numDynDummy     = MR.Parameter.GetValue( 'MP_nr_dummy_dynamic_scans' );
switch patchVersion
    case 'PIH1'
        if ( isVerbose )
            fprintf( '    patch version:    %s  \n', patchVersion )
            fprintf( '    aquisition order: dummy, acq_1, dummy, acq_2, ..., dummy, acq_n, \n' )
            fprintf( '                      dummy, trn_1, dummy, trn_2, ..., dummy, trn_n  \n' )
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
    % Slice Timea Offset
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
