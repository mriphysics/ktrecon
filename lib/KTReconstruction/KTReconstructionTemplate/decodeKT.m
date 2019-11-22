function x=decodeKT(x,SH,A,M,P)

%DECODEKT   Forward model of the reconstruction in the presence of
%motion
%   X=DECODEKT(X,{SH},{A},{M},{P})  
%   * X is the image
%   * {SH} are the conjugated sensitivities
%   * {A} contains the sampling scheme or the size of the Fourier space
%   * {M} applies a mask
%   * {P} applies a preconditioner
%   ** X is the encoded data
%

if nargin<2;SH=[];end
if nargin<3;A=[];end
if nargin<4;M=[];end
if nargin<5;P=[];end

gpu=isa(x,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

NY=size(x);
if ~isempty(SH);NX=size(SH);else NX=NY;end
if ~isempty(A);x=sum(bsxfun(@times,x,conj(A)),9);end%A^H
for m=2:2;x=ifftGPU(x,m,gpuF);end%F^H
for m=2:2;x=ifold(x,m,NX(m),NY(m));end%U^H
if ~isempty(SH);x=sum(bsxfun(@times,x,SH),4);end%S^H
if ~isempty(M);x=bsxfun(@times,x,M);end%M
if ~isempty(P);x=ifftGPU(bsxfun(@times,P,fftGPU(x,5,gpuF)),5,gpuF);end
