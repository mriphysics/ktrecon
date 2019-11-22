function x=encodeKT(x,S,A,y)

%ENCODEKT   Forward model of the reconstruction in the presence of motion
%   X=ENCODEKT(X,{S},{A},{Y})  
%   * X is the image
%   * {S} is the coil array sensitivity map
%   * {A} contains the sampling scheme or the size of the Fourier space
%   * {Y} contains the samples, for computing the residuals
%   ** X is the encoded data
%

if nargin<2;S=[];end
if nargin<3;A=[];end
if nargin<4;y=[];end
gpu=isa(x,'gpuArray');
if gpu;gpuF=2;else gpuF=0;end

NX=size(x);NX(end+1:3)=1;
if ~isempty(A);NY=size(A);else NY=NX;end

MB=size(A,3);
if MB>1
   NE=size(x,3)/MB;
   if ~isempty(y);y=resPop(y,3,[],7);end
   if ~isempty(S);S=resPop(S,3,[NE MB],7:8);end
   x=resPop(x,3,[NE MB],7:8);
   A=resPop(A,3,[],8);
end 

if ~isempty(S);x=bsxfun(@times,x,S);end%S
for m=2:2;x=fold(x,m,NX(m),NY(m));end%U
for m=2:2;x=fftGPU(x,m,gpuF);end%F
if ~isempty(A);x=sum(bsxfun(@times,x,A),8);end%A
if ~isempty(y);x=bsxfun(@minus,x,y);end%To compute the residuals

if MB>1
    x=resPop(x,7,[],3);
end
