function [x,n,E,Ex]=solverKTRegular(y,S,A,W,M,xast,nX,toler)

%SOLVERKTREGULAR   Template for reconstructing using KT for regular shots
%   X=SOLVERKTREGULAR(Y,S,{A},{W},{M},{XAST},{NX},{TOLER}) 
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

gpu=isa(y,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

%TERMS FOR REGULAR SHOTS
NY=size(y);NY(end+1:5)=1;
%R=buildStandardDFTM(NY(5),1,gpu);R=R{1};RH=R';

ND=numDims(y);
[F,FH]=build1DFTM(NY(2),0,gpu);
if (gpu && isaUnderlying(y,'double')) || isa(y,'double');F=double(F);end
if isempty(A);A=ones(1,'like',real(y));end
NA=size(A);NA(end+1:ND)=1;
AA=repmat(A,NY./NA);
y=y(AA==1);
y=reshape(y,[NY(1) NY(2)/(prod(NY)/numel(y)) NY(3:ND)]);
F=repmat(F.',[1 1 1 1 NA(5)]);
FH=repmat(FH,[1 1 1 1 NA(5)]);
NF=size(F);NF(end+1:ND)=1;
AA=repmat(A,NF./NA);
F=F(AA==1);
FH=FH(AA==1);
F=reshape(F,[NF(1) NF(2)/(prod(NF)/numel(F)) NF(3:ND)]);
FH=reshape(FH,[NF(1) NF(2)/(prod(NF)/numel(F)) NF(3:ND)]);
FC=bsxfun(@times,F,1./dynInd(F,1,5));
FC=dynInd(FC,1,2);
FC=permute(FC,[2 1 3:ND]);
FHC=conj(FC);
F=dynInd(F,1,5);
FH=dynInd(FH,1,5).';

n=0;%Effective iterations

%PRECONDITIONER
[P,SH]=precondKT(S,M);

if ~isempty(xast);y=-encodeKTRegular(xast,S,F,FC,y);end

%DECODE
r=decodeKTRegular(y,SH,FH,FHC,W);n=n+1;
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
    Ex=encodeKTRegular(x,S,F,FC,y);
    E=normm(Ex,[]);
    Ex=decodeKTRegular(Ex,[],FH,FHC);
    Ex=normm(Ex,[],4);
    n=n+1;
end

if ~isempty(xast);x=x+xast;end

function x=systemKT(x)
    xS=encodeKTRegular(x,S,F,FC);
    xS=decodeKTRegular(xS,SH,FH,FHC,W);  
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