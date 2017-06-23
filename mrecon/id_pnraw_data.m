function [ rawDataFilePath, coilSurveyFilePath, senseRefFilePath ] = id_pnraw_data( rawDataDir, seriesNo )
%ID_PNRAW_DATA  identify raw data and reference files
%
%   [ rawDataFilePath, coilSurveyFilePath, senseRefFilePath ] = ID_PNRAW_DATA( rawDataDir, seriesNo );
% 
%   rawDataDir - path to directory containing raw data for exam on pnraw01
%   seriesNo   - series number of acquisition of interest
% 
%   Example: 
%       rawDataDir  = '/home/jva13/mnt/pnraw01-archive/archive-ingenia/2015_10_27/BI_97800/';
%       seriesNo    = 28; 
%       [ r, c, s ] = ID_PNRAW_DATA( rawDataDir, seriesNo );

% jfpva (joshua.vanamerom@kcl.ac.uk)


%% Parse Input


p = inputParser;

addRequired(  p, 'rawDataDir', ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );
addRequired(  p, 'seriesNo', ...
    @(x) validateattributes( x, {'numeric'}, {'scalar'}, mfilename) );

parse( p, rawDataDir, seriesNo );


%% Identify All Files in Raw Data Directory


if ~exist( rawDataDir, 'dir' ),
    error( 'Cannot access raw data: %s', rawDataDir ),
end

D = dir( fullfile( rawDataDir, '*.lab' ) );

for iD = 1:numel(D),
    
    C = strsplit( D(iD).name, '_' );
    D(iD).seriesNo   = str2double( C{4} );
    D(iD).seriesDesc = C{6}(1:(end-4));
    
    if D(iD).seriesNo ~= 1000,
        D(iD).isCoilSurvey = false;
        D(iD).isSenseRef   = false;
    else
        if ~isempty( strfind( D(iD).seriesDesc, 'coilsurvey' ) )
            D(iD).isCoilSurvey = true;
            D(iD).isSenseRef   = false;
        elseif ~isempty( strfind( D(iD).seriesDesc, 'senseref' ) )
            D(iD).isCoilSurvey = false;
            D(iD).isSenseRef   = true;
        else
            D(iD).isCoilSurvey = false;
            D(iD).isSenseRef   = false;
        end
    end
    
end


%% Identify Data File


indRawData = find( [ D.seriesNo ] == seriesNo );
rawDataFilePath = fullfile( rawDataDir, D(indRawData).name );


%% Identify Coil Survey File


% NOTE: use most recent coil survey preceeding seriesNo

indCoilSurvey = find( [ D.isCoilSurvey ] & (1:numel(D)) < indRawData, 1, 'last' );
coilSurveyFilePath = fullfile( rawDataDir, D(indCoilSurvey).name );



%% Identify Sense Ref File

% NOTE: use most recent sense ref preceeding seriesNo

indSenseRef = find( [ D.isSenseRef ] & (1:numel(D)) < indRawData, 1, 'last' );
senseRefFilePath = fullfile( rawDataDir, D(indSenseRef).name );


%% Check Validity of Files

if ~exist( rawDataFilePath, 'file' )
    warning( 'raw k-t data file does not exist: %s', rawDataFilePath )
end

if ~exist( coilSurveyFilePath, 'file' )
    warning( 'coil survey file does not exist: %s', coilSurveyFilePath )
end

if ~exist( senseRefFilePath, 'file' )
    warning( 's file does not exist: %s', senseRefFilePath )
end


%% Verbose Output

fprintf( 'raw k-t data:    %s\n', rawDataFilePath );
fprintf( 'coil suvey data: %s\n', coilSurveyFilePath );
fprintf( 'sense ref data:  %s\n\n', senseRefFilePath );


end  % id_pnraw_data(...)
