function [x,n,E,Ex]=solverKT(y,S,A,W,M,xast,nX,toler)

%SOLVERKT   Template for reconstructing using KT
%   X=SOLVERKT(Y,S,{A},{W},{M},{XAST},{NX},{TOLER}) 
%   * Y is the measured data (X-KY-Z-C-T)
%   * S is the coil-array sensitivity map (X-Y-Z-C-T)
%   * {A} is a sampling mask (1-KY-1-1-T, 1/0)
%   * {W} is a spatial mask to constrain the solution (X-Y-Z, 1/0)
%   * {M} is a regularization mask (1-KY-1-1-F, 1/0), this is the M^{-2} in
%   the paper
%   * {XAST} is a "hard" prior regularization mask (X-Y-Z-T)
%   * {NX} is the number of iterations of the CG algorithm
%   * {TOLER} is the maximum update for convergence. It defaults to 0
%   ** X the reconstructed image
%   ** N is the number of iterations till convergence
%   ** E is the energy of the residuals
%   ** Ex are the chi2 residuals
%

if nargin<3;A=[];end
if nargin<4;W=[];end
if nargin<5 || isempty(M);M=1e-6;end
if nargin<6;xast=[];end
if nargin<7 || isempty(nX);nX=5000;end
if nargin<8 || isempty(toler);toler=1e-1;end

NY=size(y);NY(end+1:5)=1;
gpu=isa(y,'gpuArray');if gpu;gpuF=2;else gpuF=0;end
%R=buildStandardDFTM(NY(5),1,gpu);R=R{1};RH=R';
n=0;%Effective iterations

MB=size(A,3);
if MB>1
   NSL=size(S,3);
   NE=NSL/MB;
   y=resPop(y,3,[],7);
   S=resPop(S,3,[NE MB],7:8);
   if ~isempty(xast);xast=resPop(xast,3,[NE MB],7:8);end
   if size(M,3)==NSL;M=resPop(M,3,[NE MB],7:8);end
   if ~isempty(W);W=resPop(W,3,[NE MB],7:8);end
   A=resPop(A,3,[],8);
end 
ND=12;%Maximum dimensionality

%PRECONDITIONER
[P,SH]=precondKT(S,M);

if ~isempty(xast);y=-encodeKT(xast,S,A,y);end

%DECODE
r=decodeKT(y,SH,A,W);n=n+1;
NX=size(r);
x=zeros(NX,'like',r);
%z=aplGPU(RH,bsxfun(@times,P,aplGPU(R,r,5)),5);
z=ifftGPU(bsxfun(@times,P,fftGPU(r,5,gpuF)),5,gpuF);
%z=bsxfun(@times,P,r);

p=z;
rsold=multDimSum(conj(z).*r,1:ND);
l=true;
if sqrt(min(abs(rsold(:))))<1e-9;l=false;end

%ITERATIONS
while l
    %SYSTEM MATRIX
    Ap=systemKT(p);n=n+2;
    
    %UPDATES
    al=conj(rsold)./multDimSum(conj(p).*Ap,1:ND);
    xup=bsxfun(@times,al,p);
    x=x+xup; 
    xup=max(abs(xup(:)).^2);
    if xup<toler || n>=nX;break;end
    r=bsxfun(@minus,r,bsxfun(@times,al,Ap));
    %z=aplGPU(RH,bsxfun(@times,P,aplGPU(R,r,5)),5);
    z=ifftGPU(bsxfun(@times,P,fftGPU(r,5,gpuF)),5,gpuF);
    %z=bsxfun(@times,P,r);

    rsnew=multDimSum(conj(z).*r,1:ND);
    if sqrt(min(abs(rsnew(:))))<1e-9;break;end
    be=bsxfun(@times,rsnew,1./rsold);
    p=z+bsxfun(@times,be,p);
    rsold=rsnew;
end

%LOSS
if nargout>=3
    Ex=encodeKT(x,S,A,y);
    E=normm(Ex,[]);
    Ex=decodeKT(Ex,[],A);
    Ex=normm(Ex,[],4);
    if MB>1;Ex=resPop(Ex,[7 8],[],3);end
    n=n+1;
end

if ~isempty(xast);x=x+xast;end
if MB>1;x=resPop(x,[7 8],[],3);end

function x=systemKT(x)
    xS=encodeKT(x,S,A);
    xS=decodeKT(xS,SH,A);  
    if ~isempty(M)
        %xR=aplGPU(RH,bsxfun(@times,M,aplGPU(R,x,5)),5);
        xR=ifftGPU(bsxfun(@times,M,fftGPU(x,5,gpuF)),5,gpuF);
    else
        xR=0;
    end
    %xR=0;
    x=xS+xR;
    if ~isempty(W);x=bsxfun(@times,x,W);end
end

end