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
%       datatype        - 'magnitude' (default), 'phase', 'real', 'imaginary', or 'complex'
%       displayrange    - NIfTI header [min,max] values
%       datascaling     - [interept,slope], i.e., im = slope * MR.Data + intercept

% JFPvA (joshua.vanamerom@kcl.ac.uk)

%% Default Values

default.frameDuration  = 1;
default.dataType       = 'magnitude';
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
add_param_fn(   p, 'datatype', default.dataType, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );
add_param_fn(   p, 'displayrange', default.displayRange, ...
    @(x) validateattributes( x, {'numeric'}, {'vector','numel',2}, mfilename) );
add_param_fn(   p, 'datascaling', default.dataScaling, ...
    @(x) validateattributes( x, {'numeric'}, {'vector','numel',2}, mfilename) );

parse( p, MR, niiFilePath, varargin{:} );
frameDuration  = p.Results.frameduration;
dataType       = p.Results.datatype;
displayRange   = p.Results.displayrange;
dataScaling    = p.Results.datascaling;


%% Anonymous Functions

reshape_mrecon_to_nii = @( d ) reshape( permute( d, [ 1 2 8 5 3 4 6 7 9 10 11 12 ] ), size(d,1), size(d,2), size(d,8), size(d,5), [] ); 
    % mrecon: x-y-z-chan-dyn-card-echo-loca-mix-extr1-extr2-aver
    %         1-2-3-4----5---6----7----8----9---10----11----12--
    %
    % nii:    x-y-loca-dyn-other

switch dataType
    case 'magnitude'
        data2im = @(data) abs(data);
    case 'phase'
        data2im = @(data) angle(data);
    case 'real'
        data2im = @(data) real(data);
    case 'imaginary'
        data2im = @(data) imag(data);
    case 'complex'
        data2im = @(data) complex(real(data),imag(data));
    otherwise
        error( 'dataType ''%s'' not recognised', dataType )
end


%% Append .nii to FileName

if ~strcmp( niiFilePath((end-4):end), '.nii' )
    niiFilePath = strcat( niiFilePath, '.nii' );
end


%% Compute Image to World Transformation

A  = mrecon_image_to_nii_world_transform( MR );

% Sweep Adjustment
if isSweepAcq

    numLoca = size( MR.Data, 8 );
    FOV = MR.Parameter.Scan.FOV;
%     Gap = unique( MR.Parameter.Scan.SliceGap ); %NB: no gap in Sweep acquisitions, just faster Sweep rate if gap > 0
    
    A(1,3) = (FOV(3) + 1) / numLoca; % +1 because counts from zero

end

% Calculate Scaling (Voxel Size)
S  = sqrt( sum( A(1:3,1:3).^2, 1 ) );


%% Format Image Data

im = dataScaling(2) * data2im( reshape_mrecon_to_nii( MR.Data ) ) + dataScaling(1);


%% TAR - Format data to match ISDPACS nifti output:
% 
% % Swap x/y (M/P) in A so it matches .nii format output by ISDPACS
% A(:,1:2) = flip(A(:,1:2),2);   
% % end TAR.
% 
% % permute x/y in image so it matches .nii output by ISDPACS
% im = permute(im,[2,1,3]);
% 
% % warning
% warning('Using adjustments to Josh`s mrecon toolbox: TAR nifti affine correction and image permute.');
%
% end TAR


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