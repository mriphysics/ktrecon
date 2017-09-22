function mrecon_setreconparam( MR, varargin )
%MRECON_SETRECONPARAM  set MRecon object reconstruction parameters
%
%   MRECON_SETRECONPARAM( MR )
%
%   MRECON_SETRECONPARAM( ... , 'param', val) 
%
%   Input:
%       MR              - ReconFrame MRecon object
%
%   Parameter-Value Options:
%       optionpairs     - user-specified option names and values as 
%                         n x 2 cell array of character vectors
%
%   Example
%       opts = { 'GeometryCorrection', 'Yes' };
%       MRECON_SETRECONPARAM( MR, 'optionpairs', opts )

% JFPvA (joshua.vanamerom@kcl.ac.uk)


%% Default Values

default.optionpairs  = {};


%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(    p, 'MR', ...
    @(x) validateattributes( x, {'MRecon'}, {'scalar'}, mfilename) );
add_param_fn(   p, 'optionpairs', default.optionpairs, ...
    @(x) validateattributes( x, {'cell'}, {}, mfilename) );

parse( p, MR, varargin{:} );
optionpairs = p.Results.optionpairs;


%% Set Reconstruction Parameters

% Set Reconstructed Voxel Size Equal to Acquired
MR.Parameter.Scan.RecVoxelSize = MR.Parameter.Scan.AcqVoxelSize;

% Set Recon Options
MR.Parameter.Recon.kSpaceZeroFill        = 'No';
MR.Parameter.Recon.ImageSpaceZeroFill    = 'No';
MR.Parameter.Recon.RingingFilter         = 'No';
MR.Parameter.Recon.Gridding              = 'No';
MR.Parameter.Recon.CoilCombination       = 'No';
MR.Parameter.Recon.SENSE                 = 'No';
MR.Parameter.Recon.GeometryCorrection    = 'No';
MR.Parameter.Recon.RotateImage           = 'No';

% Set User-Specified Recon Options
validOptname = fieldnames( MR.Parameter.Recon );
for iOpt = 1:size(optionpairs,1)
    optname = optionpairs{iOpt,1};
    optval  = optionpairs{iOpt,2};
    if any( ismember( optname, validOptname ) )
        if ( ( ischar( MR.Parameter.Recon.(optname) ) == ischar( optval ) ) && ...
             ( isnumeric( MR.Parameter.Recon.(optname) ) == isnumeric( optval ) ) && ...
             ( iscell( MR.Parameter.Recon.(optname) ) == iscell( optval ) ) )
            MR.Parameter.Recon.(optname) = optval;
        else
            warning( 'Reconstruction option ''%s'' value data type mismatch; skipping.', optname )
        end
    else
        warning( 'Reconstruction option ''%s'' does not exist; skipping', optname )
    end
end


end  % mrecon_setreconparam(...)