function mrecon_setreconparam( MR )

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

end  % mrecon_setreconparam(...)
