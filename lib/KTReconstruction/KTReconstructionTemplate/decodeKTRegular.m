function x=decodeKTRegular(x,SH,FH,FHC,M,P)

%DECODEKTREGULAR   Forward model of the reconstruction in the presence of
%motion for regular shots
%   X=DECODEKTREGULAR(X,{SH},{FH},{FHC},{M},{P})  
%   * X is the image
%   * {SH} are the conjugated sensitivities
%   * {F} contains the inverse Fourier operator of the first shot
%   * {FC} contains the conjugate phase correction for the other shots
%   * {M} applies a mask
%   * {P} applies a preconditioner
%   ** X is the encoded data
%

if nargin<2;SH=[];end
if nargin<3;FH=[];end
if nargin<4;FHC=[];end
if nargin<5;M=[];end
if nargin<6;P=[];end

gpu=isa(x,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

if ~isempty(FH);x=matfun(@mtimes,x,FH);end
if ~isempty(SH);x=sum(bsxfun(@times,x,SH),4);end%S^H
if ~isempty(FHC);x=bsxfun(@times,x,FHC);end
if ~isempty(M);x=bsxfun(@times,x,M);end%M
if ~isempty(P);x=ifftGPU(bsxfun(@times,P,fftGPU(x,5,gpuF)),5,gpuF);end
