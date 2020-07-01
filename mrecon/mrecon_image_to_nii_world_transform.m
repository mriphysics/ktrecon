function A = mrecon_image_to_nii_world_transform( MR, stackNo )

if ~exist( 'stackIndex', 'var' )
    stackNo = 1;
end

if all( MR.Parameter.ReconFlags.isimspace )
    
    % Field of View
    FOV = [ 0 0 0 ];
    FOV(1) = MR.Parameter.Scan.RecVoxelSize(1) * size(MR.Data,1);
    FOV(2) = MR.Parameter.Scan.RecVoxelSize(2) * size(MR.Data,2);
    FOV(3) = ( MR.Parameter.Scan.RecVoxelSize(3) + MR.Parameter.Scan.SliceGap(1) ) ...
        * MR.Parameter.Scan.ImagesPerStack - MR.Parameter.Scan.SliceGap(1);               
     
    % Set Current Field of View
    curFOV = MR.Parameter.Scan.curFOV;
    MR.Parameter.Scan.curFOV = FOV;
    
    % Get Image to World Transformation
    Aim2worldRAF = MR.Transform( 'ijk', 'RAF', stackNo );  % RL-AP-FH / RAI
    
    % Reset Current Field of View
    MR.Parameter.Scan.curFOV = curFOV;
        
    % Flip to Match NIfTI Orientation
    AflipRL = [ -1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1 ];
    AflipAP = [ 1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1 ];
    A       = AflipRL * AflipAP * Aim2worldRAF;  % NIfTI expects LR-PA-FH
    
%     % TAR:
%     % LAS (aka: LR-AP-HF, aka: Radiological orientation)
%     AflipRL = [ -1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1 ];
%     AflipFH = [ 1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 1 ];
%     A       = AflipRL * AflipFH * Aim2worldRAF;
%     warning('Non-default edit: Using LAS affine in NIFTI.');
%     end TAR
    
    % Shift Origin
    origin = A * [1 1 1 1]';
    A(1:3,4) = origin(1:3);
    
else
    
    A = eye( 4 );
    warning( 'Data is not in image space, returning A = I' );
    
end

end  % mrecon_image_to_nii_world_transform(...)