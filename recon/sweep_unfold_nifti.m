function nii3d = sweep_unfold_nifti( niiFilePath, niiOutputFilePath )

%SWEEP_UNFOLD_NIFTI  convert Sweep data from 4D to 3D nifti and save new file
%   
%   save nifti file where the 3rd dimension is space and time combined
% 

%   TAR (t.roberts@kcl.ac.uk)

%% Load 4D nifti

nii4d = load_untouch_nii( niiFilePath );

nX   = size( nii4d.img, 1 );
nY   = size( nii4d.img, 2 );
nZ   = size( nii4d.img, 3 );
nT   = size( nii4d.img, 4 );
nZnT = nZ*nT;

% Create 3D nifti
nii3d     = nii4d;
nii3d.img = reshape( permute(nii4d.img,[1,2,4,3]), nX, nY, nZnT );

% Adjust 4d parameters to 3d parameters
nii3d.hdr.dime.dim(1) = 3;
nii3d.hdr.dime.dim(4) = nZnT;
nii3d.hdr.dime.dim(5) = 1;

% nii3d.hdr.dime.pixdim(1) = 0; % Not sure necessary? 0 or 1?
nii3d.hdr.dime.pixdim(4) = nii4d.hdr.dime.pixdim(4) / nT;

% Save
if nargin == 1
    save_untouch_nii( nii3d, [niiFilePath(1:end-7) '_swp3d.nii.gz' ] );
else
    save_untouch_nii( nii3d, niiOutputFilePath );
end


end