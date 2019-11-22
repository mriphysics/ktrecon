function [y,T,xou]=groupwiseSliceRegistration(y,W,T,lambda,MB,meanT,fracOrd)

%GROUPWISESLICEREGISTRATION   Per-excitation groupwise registration of dynamic studies
%   [Y,T]=GROUPWISESLICEREGISTRATION(Y,W,{T},{LAMBDA},{MB},{MEANT},{FRACORD})
%   * Y is the data to be registered
%   * W is a mask for the ROI to compute the metric
%   * {T} are the initial motion parameters (defaults to 0)
%   * {LAMBDA} is the regularizer strength (defaults to 1)
%   * {MB} is the multiband factor of the study (defaults to 1)
%   * {MEANT} is a flag to project the solution to the space where the sum of the transform parameters through the volumes is 0, defaults to 1
%   * {FRACORD} is the fractional finite difference metric order
%   * Y is the registered data
%   * T are the estimated motion parameters
%   * XOU is the average and std in time
%

NY=size(y);NY(end+1:8)=1;ND=max(numDims(y),3);
NYV=NY(1:3);

if ~exist('T','var') || isempty(T);T=single(zeros([ones(1,ND-1) NY(ND) NY(3)/MB 6]));end
if ~exist('lambda','var') || isempty(lambda);lambda=1;end%With 0.01/0.1 no perceptible degradation. With 1 a tiny bit. With 10 clear degradation. With 0 is the only case on which a zero stripe appears in the data! We set it to the biggest without degradation
if ~exist('MB','var') || isempty(MB);MB=1;end
if ~exist('meanT','var') || isempty(meanT);meanT=1;end
if ~exist('fracOrd','var') || isempty(fracOrd);fracOrd=0;end

gpu=isa(y,'gpuArray');if gpu;gpuF=2;else gpuF=0;end
NT=size(T);ndT=ndims(T);

[~,kGrid,rkGrid]=generateTransformGrids(NYV,gpu,[],[],1);
[FT,FTH]=buildStandardDFTM(NYV,0,gpu);

if fracOrd~=0;GRes=buildFilter(NYV(1:2),'FractionalFiniteDiscreteIso',ones(1,2),gpu,fracOrd);end%For fractional finite difference motion estimation

Mk=single(eye(NY(3)/MB));

x=single(zeros(NY(1:ND-1)));
convT=single(false([1 prod(NT(ndT-2:ndT-1))]));
if gpu;Mk=gpuArray(Mk);x=gpuArray(x);convT=gpuArray(convT);end
permM=1:ND+1;permM(1:2)=3:4;permM(3)=1;permM(4:ND)=5:ND+1;permM(ND+1)=2;
Mk=permute(Mk,permM);
Mk=repmat(Mk,[ones(1,ND-1) MB 1]);
Mk=reshape(Mk,[ones(1,2) NY(3) ones(1,ND-3) NY(3)/MB]);
Hsm=single(zeros([NY(3) 2]));
Hsm(1,:)=1;
Hsm(2,1)=-1;Hsm(end,2)=-1;
if gpu;Hsm=gpuArray(Hsm);end
Hsm=fftGPU(Hsm,1,gpuF,FT{3}/sqrt(NY(3)));
%Hsm=buildFilter(NY(3),'2ndFiniteDiscrete',1,gpu);
Hsm=sum(Hsm.*conj(Hsm),2);
Hsm=permute(Hsm,[2 3 1]);

%PARAMETERS
lambda=lambda/NT(ndT-2);
nX=3;
toler=1;%1;
a=[1 2 3 1 1 2 1 2 3 1 2 3 1 2 3 4 4 5 4 5 6;
   1 2 3 2 3 3 4 4 4 5 5 5 6 6 6 5 6 6 4 5 6];

NHe=size(a,2);
dHe=single(zeros([NHe prod(NT(1:ndT-1))]));
dH=single(zeros([NT(ndT) prod(NT(1:ndT-1))]));dHEff=dH;

multA=1.2;multB=2;%Factors to divide/multiply the weight that regularizes the Hessian matrix when E(end)<E(end-1)
fina=0;
winic=1;
w=winic*ones(NT(1:ndT-1));
flagw=zeros(NT(1:ndT-1));
w=reshape(w,[ones(1,ND) prod(NT(1:ndT-1))]);
flagw=reshape(flagw,[ones(1,ND) prod(NT(1:ndT-1))]);
Eprev=single(zeros([1 prod(NT(1:ndT-1))]));E=single(zeros([1 prod(NT(1:ndT-1))]));
if gpu;[E,Eprev]=parUnaFun({E,Eprev},@gpuArray);end

%Iterations
BlSz=NT(ndT-1);%BlSz=NT(ndT-2:ndT-1);
NTst=prod(NT(ndT-2:ndT-1));
NTR=[NT(1:ndT-3) NTst NT(ndT)];
cont=0;
NcontTest=1;
while fina~=2
    dHe(:)=0;dH(:)=0;dHEff(:)=0;Eprev=E;    
    if mod(cont,NcontTest)==0
        x=solveXSliceRegistration(x,y,T,Mk,Hsm,lambda,nX,toler);%%%THIS NEEDS ATTENDANCE!
        convT(:)=0;
        cont=0;NcontTest=NcontTest+1;        
    end    
    cont=cont+1;        
    ind2Est=find(~convT);NEst=length(ind2Est);    
    T=reshape(T,NTR); 
    Mk=resPop(Mk,ndT-1,[],ndT-2);
    for s=1:BlSz:NEst;vS=s:min(s+BlSz-1,NEst);vS=ind2Est(vS);
        [ind{1},ind{2}]=ind2sub(NT(ndT-2:ndT-1),vS);             
        [et,etg]=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(T,vS,ndT-2),1,1,[],1);                                
        [xT,xB]=sincRigidTransform(x,et,1,FT,FTH);
        WT=abs(sincRigidTransform(W,et,1,FT,FTH));      
        xT=bsxfun(@minus,xT,dynInd(y,ind{1},ndT-2));%Using the slice thickness may involve different "masks" for x and y
        if fracOrd~=0;xT=filtering(xT,GRes);end
        xT=bsxfun(@times,xT,WT);
        xT=bsxfun(@times,xT,dynInd(Mk,ind{2},ndT-2));      
        Eprev(1,vS)=reshape(multDimSum(real(xT.*conj(xT)),1:ndT-3),[1 length(vS)]);
        G=sincRigidTransformGradient(xB,et,etg,FT,FTH);
        GH=cell(1,NT(ndT));
        for m=1:NT(ndT)
            if fracOrd~=0;G{m}=filtering(G{m},GRes);end
            G{m}=bsxfun(@times,G{m},WT);
            G{m}=bsxfun(@times,G{m},dynInd(Mk,ind{2},ndT-2));
            GH{m}=conj(G{m});
        end
        for m=1:NHe
            GGE=real(G{a(1,m)}.*GH{a(2,m)});           
            dHe(m,vS)=gather(reshape(multDimSum(GGE,1:ndT-3),[1 length(vS)]));            
        end   
        for m=1:NT(ndT)       
            G{m}=real(GH{m}.*xT);         
            dH(m,vS)=gather(reshape(multDimSum(G{m},1:ndT-3),[1 length(vS)]));
        end
    end
            
    MHe=single(eye(NT(ndT)));
    flagw(:)=0;    
    fina=0;
    while fina==0
        for s=ind2Est
            for k=1:NHe
                if a(1,k)==a(2,k)
                    MHe(a(1,k),a(2,k))=(1+w(s))*dHe(k,s);
                else
                    MHe(a(1,k),a(2,k))=dHe(k,s);MHe(a(2,k),a(1,k))=dHe(k,s);
                end                  
                dHEff(:,s)=-winic*single(double(MHe)\double(dH(:,s)))/w(s);
            end  
        end   
        dHEff(:,w>1e10)=0;
        Tupr=reshape(permute(dHEff,[2 1]),NTR);
        Tup=T+Tupr; 
        Tup=restrictTransform(Tup);
        
        ind2EstI=find(~convT & (flagw(:)')~=2);NEst=length(ind2EstI);%Quicker
        %ind2EstI=find(~convT);NEst=length(ind2EstI);%Quicker
        for s=1:BlSz:NEst;vS=s:min(s+BlSz-1,NEst);vS=ind2EstI(vS);
             [ind{1},ind{2}]=ind2sub(NT(ndT-2:ndT-1),vS);          
             et=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(Tup,vS,ndT-2),1,0,[],1);                
             xT=sincRigidTransform(x,et,1,FT,FTH);       
             WT=abs(sincRigidTransform(W,et,1,FT,FTH));
             xT=bsxfun(@minus,xT,dynInd(y,ind{1},ndT-2));%Using the slice thickness may involve different "masks" for x and y
             if fracOrd~=0;xT=filtering(xT,GRes);end
             xT=bsxfun(@times,xT,WT);
             xT=bsxfun(@times,xT,dynInd(Mk,ind{2},ndT-2));                            
             E(1,vS)=reshape(multDimSum(real(xT.*conj(xT)),1:ndT-3),[1 length(vS)]);       
        end

        E(w>1e10)=Eprev(w>1e10);
        flagw(E<=Eprev)=2;  
        fprintf('Energy before: %0.6g / Energy after: %0.6g\n',sum(Eprev),sum(E));
        if any(flagw==1 | flagw==0)  
            w(E>Eprev & ~convT)=w(E>Eprev & ~convT)*multB;
        else
            w(~convT)=w(~convT)/multA;
            w(w<1e-8)=multA*w(w<1e-8);%To avoid numeric instabilities 
            T=Tup;
            fina=2;                     
            traDiff=dynInd(Tupr,1:3,ndT-1);%aux(:,:,1:3);
            rotDiff=convertRotation(dynInd(Tupr,4:6,ndT-1),'rad','deg');  
            [traDiffMaxC,rotDiffMaxC]=parUnaFun({traDiff,rotDiff},@max,[],ndT-2);
            [traDiffMaxS,rotDiffMaxS]=parUnaFun({traDiff,rotDiff},@max,[],ndT-1);
            fprintf('Maximum change in translation (vox): ');fprintf('%0.3f ',traDiffMaxC(:));
            fprintf('/ Maximum change in rotation (deg): ');fprintf('%0.3f ',rotDiffMaxC(:));fprintf('\n');                       
            traLim=0.16;rotLim=0.08;
            acce=1;
            traLim=traLim*acce;rotLim=rotLim*acce;
            convT(traDiffMaxS(:)<traLim & rotDiffMaxS(:)<rotLim)=1;
            fprintf('Not converged motion states: %d of %d\n',NTst-sum(single(convT)),NTst);  
            if any(~convT(:));fina=1;end            
        end
    end
    T=reshape(T,NT);
    Mk=resPop(Mk,ndT-2,[],ndT-1);
    if meanT
        aux=reshape(T,[prod(NT(1:ndT-1)) NT(ndT)]);
        aux=mean(aux,1);
        aux=reshape(aux,[ones(1,ndT-1) NT(ndT)]);
        T=bsxfun(@minus,T,aux);
    end
end

%xou=x;
%RECONSTRUCTION
nX=50;x(:)=0;
lambda=lambda*NT(ndT-2);
for m=1:NY(ND);y=dynInd(y,m,ND,solveXSliceRegistration(x,dynInd(y,m,ND),dynInd(T,m,ndT-2),Mk,Hsm,lambda,nX,toler));end
xou=mean(y,ND);%MEANING OF AVERAGE SLIGHTLY CHANGED, SEE COMMENTED LINE BEFORE
xou=cat(ND,xou,std(y,0,ND));

