
%% Load RCN / DC without performing k-t Reconstruction

% e.g.:
% - intercept mrecon_kt.m in debug mode before k-t reconstruction performed
% - run this code, if ktrecon already performed
% --- requires 's*_slice*_recon.mat' files
%
%

%% Series

ktreconDir = '/fastscratch/tr17/data/sweep_4dflow/c_fcmr333/ktrecon_swp_w96_s32_o0_needsGeoCorrFix';
cd(ktreconDir);

seriesNo = 24;

reconMatFiles = dir(['s' num2str(seriesNo) '_slice*_recon.mat']);
numSlices = numel(reconMatFiles);


%% Load ktRcn/ktDc, assign to RCN / DC

for iSlice = 1:numSlices
    
    load( reconMatFiles(iSlice).name, 'ktRcn', 'ktDc' );
    
    % Put Data in Back in MRecon Objects
    RCN.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:) = single( swap_dim_xydcl_to_reconframe( ktRcn ) );
    DC.Data(:,:,:,:,:,:,:,iSlice,:,:,:,:)  = single( swap_dim_xydcl_to_reconframe( ktDc ) );

    clear ktRcn ktDc
    
end

disp( 'All slices loaded.' );