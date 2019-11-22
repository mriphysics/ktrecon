function S=sensitEspiritKT(y,x,B,voxsiz,lr)

%SENSITESPIRITKT   Estimate sensitivities by ESPIRIT
%   S=SENSITESPIRITKT(Y,{X},{B},{VOXSIZ},{LR})  
%   * Y is fully sampled information in image space
%   * {X} are reconstructions in image space
%   * {B} is a normalization from a separate scan
%   * {VOXSIZ} provides the voxel size
%   * {LR} indicates it is a low res spearate scan
%   ** S are the sensitivities
%

if nargin<4 || isempty(voxsiz);voxsiz=ones(1,3);end
if nargin<5 || isempty(lr);lr=0;end

gpu=isa(y,'gpuArray');
if gpu;gpuF=2;else gpuF=0;end

rec.Dyn.Debug=2;
rec.Par.Mine.pedsUn=1;
rec.Enc.AcqVoxelSize=voxsiz;
rec.Par.Labels.SliceGaps=0;
rec.y=y;
if nargin>=2 && ~isempty(x);rec.x=x;end
rec.Alg.parE.NCV=1;%Maximum number of eigenmaps
if lr
    rec.Alg.parE.NC=[10 10 10];%Resolution (mm) of calibration area to compute compression
    rec.Alg.parE.K=[100 100 100];
    rec.Alg.parE.subSp=[1 1 1];%Subsampling in image space to accelerate
    rec.Alg.parE.mirr=[0 0 0];%Whether to mirror along a given dimension
    rec.Alg.parE.Ksph=round((100)^(2/3));%Number of points for spherical calibration area, unless 0, it overrides other K's
else
    rec.Alg.parE.NC=[4 4];%Resolution (mm) of calibration area to compute compression
    rec.Alg.parE.K=[100 100];
    rec.Alg.parE.subSp=[1 1];%Subsampling in image space to accelerate
    rec.Alg.parE.mirr=[8 8];%Whether to mirror along a given dimension
    rec.Alg.parE.Ksph=round(1200^(2/3));%Number of points for spherical calibration area, unless 0, it overrides other K's
end
rec.Alg.parE.eigTh=0;%0.02%Threshold for picking singular vectors of the calibration matrix (relative to the largest singular value)
rec.Alg.parE.absPh=0;%1%Flag to compute the absolute phase
rec.Alg.parE.virCo=1;%Flag to use the virtual coil to normalize the maps
rec.Alg.parE.eigSc=[0.85 0.3];%0.25;
rec.Alg.parE.dimLoc=[];%Dimensions along which to localize to compute virtual coils
rec.Alg.parE.Kmin=6;%Minimum K-value before mirroring
rec.Dyn.Typ2Rec=[];

rec=solveESPIRIT(rec);

S=rec.S;
if nargin>=3 && ~isempty(B);S=bsxfun(@times,S,B./(sqrt(normm(S,[],4)+1e-9)));end

xr=solverKT(fftGPU(y,2,gpuF),S);
if lr
    if abs(angle(mean(xr(:).*conj(x(:)))))>pi/2;S=-S;end
else
    for s=1:size(x,3)
        xxr=dynInd(xr,s,3);xx=dynInd(x,s,3);
        if abs(angle(mean(xxr(:).*conj(xx(:)))))>pi/2;S=dynInd(S,s,3,-dynInd(S,s,3));end
    end
end
%S=bsxfun(@times,S,single(W>rec.Alg.parE.eigSc(1)));
