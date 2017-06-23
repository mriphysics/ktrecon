function niiFilePath = mrecon_writenifti( MR, niiFilePath, varargin )
%MRECON_WRITENIFTI  write MRecon object image-space data to NIfTI 
%
%   MRECON_WRITENIFTI( MR, niiFilePath )
%
%   niiFilePath = MRECON_WRITENIFTI( MR, niiFilePath )
%
%   MRECON_WRITENIFTI( ... , 'param', val) 
%
%   Input:
%       MR              - ReconFrame MRecon object
%       niiFilePath     - path to nifti file to write
%
%   Parameter-Value Options:
%       frameduration   - duration, in seconds, of one dynamic frame
%       complex         - save complex-valued images
%       displayrange    - NIfTI header [min,max] values
%       datascaling     - [interept,slope], i.e., im = slope * MR.Data + intercept

% JFPvA (joshua.vanamerom@kcl.ac.uk)

%% Default Values

default.frameDuration  = 1;
default.saveComplex    = false;
default.displayRange   = [];
default.dataScaling    = [0 1];


%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(    p, 'MR', ...
    @(x) validateattributes( x, {'MRecon'}, {'scalar'}, mfilename) );
addRequired(   p, 'niiFilePath', ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );
add_param_fn(   p, 'frameduration', default.frameDuration, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','real','scalar'}, mfilename) );    
add_param_fn(   p, 'complex', default.saveComplex, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );
add_param_fn(   p, 'displayrange', default.displayRange, ...
    @(x) validateattributes( x, {'numeric'}, {'vector','numel',2}, mfilename) );
add_param_fn(   p, 'datascaling', default.dataScaling, ...
    @(x) validateattributes( x, {'numeric'}, {'vector','numel',2}, mfilename) );

parse( p, MR, niiFilePath, varargin{:} );
frameDuration  = p.Results.frameduration;
saveComplex    = p.Results.complex;
displayRange   = p.Results.displayrange;
dataScaling    = p.Results.datascaling;
if saveComplex && isreal( MR.Data )
    saveComplex = false;
end


%% Anonymous Function

reshape_mrecon_to_nii = @( d ) reshape( permute( d, [ 1 2 8 5 3 4 6 7 9 10 11 12 ] ), size(d,1), size(d,2), size(d,8), size(d,5), [] ); 
    % mrecon: x-y-z-chan-dyn-card-echo-loca-mix-extr1-extr2-aver
    %         1-2-3-4----5---6----7----8----9---10----11----12--
    %
    % nii:    x-y-loca-dyn-other


%% Append .nii to FileName

if ~strcmp( niiFilePath((end-4):end), '.nii' )
    niiFilePath = strcat( niiFilePath, '.nii' );
end


%% Compute Image to World Transformation

A  = mrecon_image_to_nii_world_transform( MR );

% Calculate Scaling (Voxel Size)
S  = sqrt( sum( A(1:3,1:3).^2, 1 ) );


%% Format Image Data

im = single( reshape_mrecon_to_nii( MR.Data ) );
if ~saveComplex
    im = abs( im );
end
im = dataScaling(2) * im + dataScaling(1);


%% Write to File

N = make_nii( im, S );

N.hdr.dime.pixdim(5) = frameDuration;

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

N.hdr.hist.srow_x = A(1,:);
N.hdr.hist.srow_y = A(2,:);
N.hdr.hist.srow_z = A(3,:);
N.hdr.hist.sform_code = 1;

save_nii( N, niiFilePath );

system( sprintf( 'gzip -f %s', niiFilePath ) );
niiFilePath = strcat( niiFilePath, '.gz' );


end  % mrecon_writenifti(...)