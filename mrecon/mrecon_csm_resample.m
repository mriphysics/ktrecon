function [csmR, imBodyR, imCoilR, psiSensR] = mrecon_csm_resample( SENS, nZcsm, interpMethod )

%% csm_resample


%% Interpolation method
if nargin < 3
    interpMethod = 'nearest'; % default in imresize3
end


%% Anonymous function

% Swap dimensions to/from ReconFrame
%   ReconFrame: 
%       1-2-3-4----5---6----7----8----9---10----11----12
%       x?y-z?chan?dyn-card?echo?loca?mix?extr1?extr2?aver
%   Otherwise:
%       x-y-dyn-chan-loca-z-card-echo-mix-extr1-extr12-aver
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
ind_reconframe_to_xydcl = [ dim.x dim.y dim.dyn dim.chan dim.loca dim.z dim.card dim.echo dim.mix dim.extr1 dim.extr2 dim.aver ];
[~,ind_xydcl_to_reconframe] = sort( ind_reconframe_to_xydcl );
swap_dim_reconframe_to_xydcl = @( data ) permute( data, ind_reconframe_to_xydcl ); 




%% 

SENS.OutputSizeSensitivity(3) = nZcsm;
SENS.OutputSizeReformated(3)  = nZcsm;
SENS.Perform;

csm     = swap_dim_reconframe_to_xydcl( SENS.Sensitivity );             % coil sensitivity maps
imBody  = swap_dim_reconframe_to_xydcl( SENS.ReformatedBodycoilData );  % body coil image
imCoil  = swap_dim_reconframe_to_xydcl( SENS.ReformatedCoilData );      % array coil images
psiSens = SENS.Psi;                                                     % array coil noise covariance


%% Resample
interpDim = [size(csm,1) size(csm,2) size(csm,3) size(csm,4), nZcsm ];

csmR    = resamp(    csm, interpDim, interpMethod );
imBodyR = resamp( imBody, interpDim, interpMethod );
imCoilR = resamp( imCoil, interpDim, interpMethod );


%% View
% montage_RR(squeeze(abs(csm(:,:,1,1,1:10:end))),'',[0 1]);
% montage_RR(squeeze(abs(csmR(:,:,1,1:100:end))),'',[0 1]);
% montage_RR(squeeze(abs(csm(:,:,1,1,1:10:end))) - squeeze(abs(csmR(:,:,1,1:100:end))) );


end % csm_resample(...)



function imr = resamp( im, dim, interpMethod )

if nargin < 3
    interpMethod = 'nearest'; % default in imresize3
end

[y, x, z, c, loc] = ndgrid( linspace( 1, size(im,1), dim(1) ),...
                            linspace( 1, size(im,2), dim(2) ),...
                            linspace( 1, size(im,3), dim(3) ),...
                            linspace( 1, size(im,4), dim(4) ),...
                            linspace( 1, size(im,5), dim(5) ) );

% nb: dim(3) size usually = 1, hence squeeze:
imr = interpn(squeeze(im), y, x, c, loc, interpMethod);
                
end