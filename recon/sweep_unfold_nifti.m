function [nii3d, swpParams] = sweep_unfold_nifti( niiFilePath, swpParams, niiOutputFilePath )

%SWEEP_UNFOLD_NIFTI  convert Sweep data from 4D to 3D nifti and save new file
%   
%   save nifti file where the 3rd dimension is space and time combined
% 

%   TAR (t.roberts@kcl.ac.uk)


%% Load 4D nifti

if strfind( niiFilePath, 'nii.gz' )
    nii4d = load_untouch_nii( niiFilePath );
elseif isempty( strfind( niiFilePath, 'nii.gz' ) )
    nii4d = load_untouch_nii( strcat( niiFilePath, 'nii.gz') );
end

nX   = size( nii4d.img, 1 );
nY   = size( nii4d.img, 2 );
nZ   = size( nii4d.img, 3 );
nT   = size( nii4d.img, 4 );
nZnT = nZ*nT;

warning('apodLength hard-coded = 0');
apodLength = 0;


%% Init 3D nifti
nii3d = nii4d;
[nii3d.img, swpParams ] = sweep_window_filter( nii4d.img, swpParams, apodLength );
% nii3d.img = reshape( permute(nii4d.img,[1,2,4,3]), nX, nY, nZnT ); % if identical size to M2D


%% Adjust 4d parameters to 3d parameters

% dim
nii3d.hdr.dime.dim(1) = 3;
nii3d.hdr.dime.dim(4) = size(nii3d.img,3);
nii3d.hdr.dime.dim(5) = 1;

% pixdim
nii3d.hdr.dime.pixdim(1) = 1; % Not sure necessary? 0 or 1?
nii3d.hdr.dime.pixdim(4) = nii4d.hdr.dime.pixdim(4) / swpParams.Encoding.NrDyn(1); % use explicit MRecon NrDyn parameter to avoid conflict with swpWinFullWidth
nii3d.hdr.dime.pixdim(5) = 1;

% affine
nii3d.hdr.hist.srow_x(3) = nii4d.hdr.hist.srow_x(3) / swpParams.Encoding.NrDyn(1);
nii3d.hdr.hist.srow_y(3) = nii4d.hdr.hist.srow_y(3) / swpParams.Encoding.NrDyn(1);
nii3d.hdr.hist.srow_z(3) = nii4d.hdr.hist.srow_z(3) / swpParams.Encoding.NrDyn(1);

% adjust z-location of first slice due to filtering/apodization
% TODO: make the filtering code and associated record keeping better
swpLoc       = unique(nonzeros(swpParams.swpWindows .* swpParams.swpWindowsFilter));
swpNewOrigin = swpLoc(1) - 1;

nii3d.hdr.hist.srow_z(4) = nii3d.hdr.hist.srow_z(4) + ( swpNewOrigin * nii3d.hdr.dime.pixdim(4) );

swpParams.swpNewOrigin = swpNewOrigin;


%% Save
if nargin < 3
    save_untouch_nii( nii3d, [niiFilePath(1:end-7) '_swp3d.nii.gz' ] );
else
    save_untouch_nii( nii3d, niiOutputFilePath );
end


% sweep_unfold_nifti(...)
end