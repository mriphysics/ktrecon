function [  ] = recon_kt_sc_slw( outputDirPath, kSpaceFilePrefix )
%RECON_KT_SC_SLW   Self-calibrated sliding window k-t reconstruction
%
%   [ ] = RECON_KT_SC_SLW(  )
%
%   output:
%       matFile containing MRecon processed k-space
%
%   input:
%       outputDirPath               path to output directory
%       kSpaceFilePrefix            k-space file name prefix
%
%   see also: mrecon_sc_slw, solverKT

% tar (t.roberts@kcl.ac.uk)
% wrapper for code by Lucilio Cordero-Grande


%% Setup

kSpaceFileSuffix = '_kspace_sc_slw.mat';
kspaceMatFilePath = fullfile( outputDirPath, strcat( kSpaceFilePrefix, kSpaceFileSuffix ) );

cd(outputDirPath);

% GPU
gpu=gpuDeviceCount;
if gpu;gpuF=2;else gpuF=0;end
if gpu;dev=gpuDevice;end


%% Reconstruction Parameters

useDynamicCoils = 0;       % Use dynamic coils flag
sigma           = sqrt(2); %2;%1;%0.2;%sqrt(2);%"Noise level" %2 high SNR/sqrt(2) intermediate SNR/1 low SNR
smoothFactor    = 1;       % Smoothing for spectral estimation
comp            = 1;       % Data compression flag

C   = [8 4 3 2 1 1]; % Sliding Window Hierarchy
dyn = []; % Select series of dynamics
sl  = []; % Select series of slices


%%  Load K-Space

fprintf( 'Loading k-space: %s ... \n', kspaceMatFilePath );
if gpu;wait(dev);end;tsta=tic;
load( kspaceMatFilePath );
if gpu;wait(dev);end;tend=toc(tsta);
fprintf('Time reading MAT: %.2fs\n',tend);

% uncompress kspace data
if comp
    NY=size(y);NY(end+1:5)=1;
    NA=size(A);NA(end+1:5)=1;
    NR=NY./NA;NR(2)=1;
    A=repmat(A,NR);
    z=A;z(:)=0;
    z(A==1)=y;
    y=z;z=[];A=[];
end


%% SOLVER KT

% Set up solver
if gpu;wait(dev);end;tsta=tic;
if ~isempty(sl);y=dynInd(y,sl,3);S=dynInd(S,sl,3);end

NY=size(y);NY(end+1:5)=1;
x=zeros([NY(1:3) length(C) NY(5)],'like',S);
if gpu;x=gpuArray(x);end

fprintf( 'Running SolverKT on %s ... \n', kSpaceFilePrefix );

% Run solver
for z=1:NY(3)
    fprintf('Slice %d of %d\n',z,NY(3));
    yz=dynInd(y,z,3);Sz=dynInd(S,z,3);
    if gpu;yz=gpuArray(yz);Sz=gpuArray(Sz);end
    
    % ESPIRIT
    B=sqrt(normm(Sz,[],4));%Body coil "reference"
    yr=scaleKT(yz,R);
    if ~useDynamicCoils
        yr=sum(yr,5)/NY(5);
        xr=solverKT(yr,Sz);
        Sz=sensitEspiritKT(ifftGPU(yr,2,gpuF),xr*sqrt(normm(yr)./normm(xr)),B,voxsiz);
    else
        xr=solverKT(yr,Sz);
        Sz=repmat(Sz,[1 1 1 1 NY(5)]);
        nor=sqrt(normm(yr,[],1:4)./normm(xr,[],1:4));
        for s=1:NY(5);Sz=dynInd(Sz,s,5,sensitEspiritKT(ifftGPU(dynInd(yr,s,5),2,gpuF),dynInd(xr,s,5)*nor(s),B,voxsiz));end
    end
    
    %HIERARCHICAL SOLVER
    % - this basically does a KT SENSE recon starting with R=1, then using
    %that reconstruction as training data for a R=2 reconstruction, which
    % is used to do an R=3 reconstruction, and so on.
    for c=C
        [yr,A]=scaleKT(yz,c);
        if c==C(1);M=[];else M=calibrationKT(xz,sqrt(R/c)*sigma,smoothFactor);end
        if c==1;[xz,n,~,Xz]=solverKTRegular(yr,Sz,A,[],M);else [xz,n,~,Xz]=solverKT(yr,Sz,A,[],M);end
        
        % c = 1 is basically SENSE. This gives starting/training for
        % following iterations
        
        % yr = data (k space only along PE direction)
        % Sz = sensitivities
        % A  = mask of sampled data
        % M  = calibration kernel (for later?)
        % xz = reconstructed data
        % n  = number of iterations
        % Xz = residuals of reconstruction (useful for later? maybe shows fat regions)
        
        if c==C(1);xTz=xz;else xTz=cat(4,xTz,xz);end
        fprintf('KT %.2f: iterations %d\n',R/c,n);
    end
    x=dynInd(x,z,3,gather(xTz));
end

if gpu;wait(dev);end;tend=toc(tsta);
fprintf('Time reconstructing: %.2fs\n',tend);


%% Save Reconstruction to Matfile

% Final Iteration of Sliding Windows
% TODO: add option to save intermediate levels of sliding window k-t reconstruction
xtRcn = squeeze(x(:,:,:,end,:));

% Reformat data to match JFPVA ktrecon code
xtRcn = gather( permute( xtRcn, [1 2 5 6 4 7 8 3]) );

% Preprocessed flag: signifies whether xtRcn has been geometry corrected,
% etc. to match jfpva reconstruction
isPreprocessed = false;

if gpu;wait(dev);end;tsta=tic;
scSlwMatFilePath = fullfile( outputDirPath, strcat( kSpaceFilePrefix, '_sc_slw_recon.mat' ) );
save( scSlwMatFilePath, 'xtRcn', 'isPreprocessed', '-v7.3' );
if gpu;wait(dev);end;tend=toc(tsta);
fprintf('Time to save sc_slw_recon MAT: %.2fs\n',tend);


%% Finished Message
fprintf( 'SolverKT complete. \n' );


end  % recon_kt_sc_slw(...)