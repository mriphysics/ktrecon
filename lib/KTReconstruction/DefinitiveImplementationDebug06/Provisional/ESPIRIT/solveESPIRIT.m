function rec=solveESPIRIT(rec)

%SOLVEESPIRIT   Estimates the sensitivities using the ESPIRIT method based
%on [1] M Uecker, P Lai, MJ Murphy, P Virtue, M Elad, JM Pauly, SS 
%Vasanawala, M Lustig, "ESPIRiT—An eigenvalue approach to autocalibrating
%parallel MRI: where SENSE meets GRAPPA," Magn Reson Med, 71:990-1001, 
%2014, [2] M Uecker, M Lustig, "Estimating absolute-phase maps using 
%ESPIRiT and virtual conjugate coils," Magn Reson Med, 77:1201-1207 (2017),
%[3] M Uecker, P Virtue, SS Vasanawala, M Lustig, "ESPIRiT reconstruction 
%using soft SENSE," ISMRM (2013), [4] M. Buehrer, P. Boesiger, S. Kozerke,
%"Virtual body coil calibration for phased-array imaging," ISMRM (2009).
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   ** REC is a reconstruction structure with estimated sensitivities rec.S
%   and eigenmaps rec.W
%

useFull=0;

debugOr=rec.Dyn.Debug;
%rec.Dyn.Debug=2;

parE=rec.Alg.parE;
gpu=isa(rec.y,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

tsta=tic;
%WE USE THE MEDIAN OVER THE REPEATS TO ESTIMATE
ND=length(parE.NC);%Number of dimensions
N=size(rec.y);N(end+1:12)=1;
NPE=min(length(rec.Par.Mine.pedsUn),prod(N(5:12)));
if isfield(rec,'x');xf=permute(rec.x,[1 2 3 5 4]);end
if NPE==1
    %for s=1:N(4);yf=dynInd(yf,s,4,multDimMed(dynInd(rec.y,s,4),5:rec.Plan.NDims));end
    %if isfield(rec,'x');xf=multDimMed(xf,5:rec.Plan.NDims);end
    yf=rec.y(:,:,:,:,1);
    if isfield(rec,'x');xf=xf(:,:,:,:,1);end
else
    yf=rec.y(:,:,:,:,:);
end
if isfield(rec,'x')
    %for m=1:ND;xf=fold(xf,m,size(x,m),N(m));end
    xf=resampling(xf,N(1:ND),2);
end

for n=1:length(rec.Alg.parE.mirr)
    if rec.Alg.parE.mirr(n)>0
        mirr=rec.Alg.parE.mirr(n);
        yext1=dynInd(yf,1:mirr,n);
        yext2=dynInd(yf,N(n)-mirr+1:N(n),n);
        yf=cat(n,flip(yext1,n),yf,flip(yext2,n));
        xext1=dynInd(xf,1:mirr,n);
        xext2=dynInd(xf,N(n)-mirr+1:N(n),n);
        xf=cat(n,flip(xext1,n),xf,flip(xext2,n));
    end
end
N=size(yf);N(end+1:12)=1;
yf=gather(yf);
if isfield(rec,'x');xf=gather(xf);end


voxsiz=rec.Enc.AcqVoxelSize;
voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);
DeltaK=(1./(N(1:ND).*voxsiz(1:ND)));

DeltaK=DeltaK(1:ND);
parE.NC=(1./parE.NC)./DeltaK;%Resolution of k-space points for calibration
parE.NC=min(parE.NC,N(1:ND)-1);%Never bigger than the image resolution
parE.NC=parE.NC-mod(parE.NC,2)+1;%Nearest odd integer
%parE.NC=parE.NC.*(1+rec.Alg.parE.mirr(1:ND));%Mirroring

if parE.Ksph>0
    parE.K=ones(1,ND)*parE.Ksph^(1/ND);
else
    parE.K=(1./parE.K)./DeltaK;%Resolution of target coil profiles
    %parE.K=parE.K-mod(parE.K,2)+1;
    parE.K=max(parE.K,parE.Kmin);
    parE.K=min(parE.K,N(1:ND));%Never bigger than the image resolution
    %parE.K=parE.K-mod(parE.K,2)+1;
    %parE.K=parE.K.*(1+rec.Alg.parE.mirr(1:ND))-rec.Alg.parE.mirr(1:ND);
end
NK=round(prod(parE.K));

if rec.Dyn.Debug>=2
    fprintf('Size of calibration area:%s\n',sprintf(' %d',parE.NC));
    fprintf('Approximate size of kernels:%s\n',sprintf(' %d',round(parE.K)));
end

tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time preprocessing: %.3f s\n',tend);end

rec.S=[];rec.W=[];
for pe=1:NPE
    y=dynInd(yf,pe,5);
    if gpu;y=gpuArray(y);end
    if isfield(rec,'x')
        x=dynInd(xf,pe,5);
        if gpu;x=gpuArray(x);end
    end
    N=size(y);N(end+1:4)=1;
    
    %VIRTUAL BODY COIL
    if parE.virCo>0
        tsta=tic;
        if parE.virCo>=2 || ~isfield(rec,'x')
            %COMPRESSION
            dimNoLoc=setdiff(1:3,parE.dimLoc);
            perm=[1 2 3 5 4];    
            yH=permute(conj(y),perm);
            Ainv=bsxfun(@times,y,(sum(abs(y).^2,4)+1e-9).^(-1));
            perm=[4 5 1 2 3];           
            P=permute(multDimSum(bsxfun(@times,Ainv,yH),dimNoLoc),perm);
            yH=[];Ainv=[];
        
            NP=size(P);
            P=resSub(P,3:5);       
            [~,P]=svdm(P);
            P=reshape(P,NP);
        
            if parE.virCo==3%WE ESTIMATE IN ORTHOGONAL SPACE
                y=permute(y,perm);
                y=matfun(@mtimes,matfun(@ctranspose,P),y);
                y=ipermute(y,perm);
                pa=y;
            else%WE ADD A BODY COIL CHANNEL
                pa=permute(y,perm);   
                pa=matfun(@mtimes,matfun(@ctranspose,P),pa);
                pa=ipermute(pa,perm);
            end
            %PHASE CORRECTION
            pa=dynInd(pa,1,4);pb=pa;pb(:)=1;        
            for n=1:length(parE.dimLoc)           
                phDif=sign(multDimSum(dynInd(pa,1:N(parE.dimLoc(n)),parE.dimLoc(n)).*conj(dynInd(pa,[1 1:N(parE.dimLoc(n))-1],parE.dimLoc(n))),setdiff(1:3,parE.dimLoc(n))));                       
                pb=bsxfun(@times,pb,conj(cumprod(phDif,parE.dimLoc(n))));
                pa=bsxfun(@times,pa,pb);
            end
            if parE.virCo==3
                y=bsxfun(@times,pb,y);
            else 
                pa=bsxfun(@times,pb,pa);
                y=cat(4,pa,y);
            end
        else
            y=cat(4,x,y);
        end
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time compressing channels: %.3f s\n',tend);end
    end
    NPm=N;
    
    
    %N(1:ND)=N(1:ND).*(1+rec.Alg.parE.mirr);
    
        
    tsta=tic;
    for n=1:ND
        %mirr=zeros(1,ND);mirr(n)=rec.Alg.parE.mirr(n);
        %y=mirroring(y,mirr,1);
        y=fftGPU(y,n,gpuF)/sqrt(N(n));
        y=fftshift(y,n);   
    end
    if parE.absPh;y=cat(4,y,fftflip(y,1:ND));end
    y=resampling(y,parE.NC,2);
    if parE.absPh;yor=y(:,:,:,1:end/2);end

    nc=size(y,4);
    if ND==2;NZ=size(y,3);else NZ=1;end
    if isempty(parE.NCV);parE.NCV=nc;end
    assert(parE.NCV<=nc,'Not prepared to deal with number of eigenmaps bigger than the number of channels');
    NMaps=max(round(N(1:ND)./parE.subSp),parE.NC);
    NMapsPm=NMaps;
    %if ~useFull;NMapsPm=NMaps./(1+rec.Alg.parE.mirr);end
    S=zeros([nc parE.NCV prod(NMapsPm) NZ],'like',rec.y);
    W=zeros([1 parE.NCV prod(NMapsPm) NZ],'like',rec.y);   
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time extracting calibration and booking memory: %.3f s\n',tend);end
        
    %HANKEL INDEXES
    if parE.Ksph>0
        NGR=ceil(NK.^(1/ND))*ones(1,ND);
        NGD=2*NGR+1;
        %WEIGHTS OF THE PATCH STRUCTURE
        xx=generateGrid(NGD,gpu,NGD,ceil((NGD+1)/2));
        r=xx{1}(1);r(1)=double(0);
        for n=1:ND;r=bsxfun(@plus,r,abs(xx{n}).^2);end
        [~,ir]=sort(r(:));
        irM=ir(1:NK);
        irs=ind2subV(NGD,irM);
        xLim=zeros(ND,2,'like',xx{1});
        for n=1:ND;xx{n}=xx{n}(irs(:,n));xx{n}=xx{n}(:);xLim(n,:)=[min(xx{n}) max(xx{n})];end 
        sw=cat(2,xx{:});
        iX=cell(1,ND);   
        for s=1:ND
            iX{s}=1-xLim(s,1):parE.NC(s)-xLim(s,2);
            iX{s}=bsxfun(@plus,sw(:,s),iX{s})';
        end
    else
        iw=1:NK;
        sw=ind2subV(parE.K,iw);
        iX=cell(1,ND);
        for s=1:ND
            iX{s}=0:parE.NC(s)-parE.K(s);
            iX{s}=bsxfun(@plus,sw(:,s),iX{s})';
        end
    end
    iXi=cell(1,ND);

    NCh=size(y,4);%Number of coils
    for z=1:NZ%Slices
        %COMPUTE CALIBRATION MATRIX
        tsta=tic;      
        if ND==2;yz=dynInd(y,z,3);else yz=y;end
        if z==NZ;y=[];end
        A=zeros([NK*NCh NK*NCh],'like',yz);
        NZZ=1;
        if ND==3;NZZ=size(iX{3},1);end
        NSa=size(iX{1},1)*size(iX{2},1);
        BlSz=1;        
        for zz=1:BlSz:NZZ;vZ=zz:min(zz+BlSz-1,NZZ);         
            Ab=zeros([NSa*length(vZ) NK NCh],'like',A);             
            for w=1:NK
                for s=1:2;iXi{s}=iX{s}(:,w);end      
                if ND==3;iXi{3}=iX{3}(vZ,w);end                
                Ab(:,w,:)=reshape(dynInd(yz,iXi,1:ND),[NSa*length(vZ) 1 NCh]);
            end
            Ab=Ab(:,:);
            A=A+Ab'*Ab;
        end
        A=(A+A')/2;      
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time building Hankel (size %d) slice %d/%d: %.3f s\n',size(A,1),z,NZ,tend);end                
        
        tsta=tic;
        A=gather(A);
        [V,D]=schur(A);
        if gpu;[V,D]=parUnaFun({V,D},@gpuArray);end
        D=diagm(sqrt(abs(D)));
        D=flip(D,2);%We use decreasingly sorted singular values
        V=flip(V,2);
        [D,iS]=sort(D,2,'descend');  
        V=indDim(V,iS,2);
        D(D<1e-9)=1;
        NV=size(V);          
        V=reshape(V,[NK NCh NV(2)]);%Kernel  
        if parE.eigTh<=0;[~,nv]=screePoint(D.^2);else nv=find(D>=D(1)*parE.eigTh,1,'last');end%Automatic/tuned detection
        if rec.Dyn.Debug==2;fprintf('Threshold for SV (relative to the largest SV): %.3f\n',D(nv)/D(1));end
        %CROP KERNELS AND COMPUTE EIGEN-VALUE DECOMPOSITION IN IMAGE SPACE TO GET MAPS
        NDV=numDims(V);
        V=dynInd(V,1:nv,NDV);        
        %ROTATE KERNEL TO ORDER BY MAXIMUM VARIANCE
        perm=[1 3 2];
        V=permute(V,perm);
        V=reshape(V,NK*nv,[]);
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing k-space kernels slice %d/%d: %.3f s\n',z,NZ,tend);end                
        if rec.Dyn.Debug>=2;fprintf('Size of k-space kernels:%s\n',sprintf(' %d',size(V)));end        
        
        tsta=tic;
        if size(V,1)<size(V,2);[~,~,v]=svd(V);else [~,~,v]=svd(V,'econ');end     
        V=V*v;
        if parE.Ksph>0
            Vaux=reshape(V,[NK nv nc]);
            K=xLim(:,2)-xLim(:,1)+1;
            V=zeros([K' nv nc],'like',Vaux);
            swi=bsxfun(@minus,sw',xLim(:,1))+1;
            swa=cell(1,ND);
            for w=1:NK
                for s=1:ND;swa{s}=swi(s,w);end
                V=dynInd(V,swa,1:ND,reshape(Vaux(w,:,:),[ones(1,ND) nv nc]));
            end
        else
            V=reshape(V,[parE.K nv nc]);
        end        
        NDV=numDims(V);
        perm=1:NDV;perm([NDV-1 NDV])=[NDV NDV-1];
        V=conj(permute(V,perm));        
        for n=1:ND;V=ifftshift(V,n);end         
        V=V*(prod(NMaps)/sqrt(NK));        
        perm(1:2)=ND+(1:2);perm(3:ND+2)=1:ND;
        V=permute(V,perm);
        if ND==3
            for n=ND
                NMapsOr=size(V);
                NMapsOr(n+2)=NMaps(n);
                [~,FH]=build1DFTM(NMapsOr(n+2),0,gpu);
                FH=resampling(FH,[NMapsOr(n+2) size(V,n+2)],1);
                %if ~useFull && rec.Alg.parE.mirr(n);FH=FH(1:end/2,:);end
                V=aplGPU(FH,V,n+2);
            end
        end
        NSp=size(V,5);
        for a=1:NSp
            Va=dynInd(V,a,5);
            for n=1:2
                NMapsOr=size(Va);
                NMapsOr(n+2)=NMaps(n);
                [~,FH]=build1DFTM(NMapsOr(n+2),0,gpu);
                FH=resampling(FH,[NMapsOr(n+2) size(Va,n+2)],1);
                %if ~useFull && rec.Alg.parE.mirr(n);FH=FH(1:end/2,:);end
                Va=aplGPU(FH,Va,n+2);
            end
            Va=resSub(Va,3:ND+2);
            NVa=size(Va,3);
            [D,U]=svdm(Va);
            D=D(:,1:parE.NCV,:);U=U(:,1:parE.NCV,:);
            D=real(D);
            if parE.virCo==1 && isfield(rec,'x');Uref=U(1,1,:);else Uref=sign(U(1,1,:));end           
            %if parE.virCo==1 && isfield(rec,'x');Uref=U(1,:,:);else Uref=sign(U(1,:,:));end           
            U=matfun(@mtimes,v,bsxfun(@rdivide,U,Uref));
            D=reshape(D,[1 parE.NCV NVa]);U=reshape(U,[nc parE.NCV NVa]);
            vA=(1:NVa)+(a-1)*NVa;
            S(:,:,vA,z)=U;W(:,:,vA,z)=D;
        end;Va=[];V=[];
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing spatial maps slice %d/%d: %.3f s\n',z,NZ,tend);end     
    end
    
    tsta=tic;
    S=ipermute(S,perm);
    W=ipermute(W,perm);


    if ND==3    
        S=reshape(S,[NMapsPm nc parE.NCV NZ]);
        W=reshape(W,[NMapsPm 1 parE.NCV NZ]);    
    else
        S=reshape(S,[NMapsPm NZ nc parE.NCV]);
        W=reshape(W,[NMapsPm NZ 1 parE.NCV]); 
        
        %perm=[1:ND ND+[3 1:2]];
        %S=permute(S,perm);
        %W=permute(W,perm);
    end

    if parE.absPh
        NS=size(S);NS(end+1:5)=1;
        S=resPop(S,4,[NS(4)/2 2]);
        NDS=numDims(S);
        ph=exp(1i*angle(sum(dynInd(S,1,NDS).*dynInd(S,2,NDS),NDS-1))/2);
        S=dynInd(S,1,NDS);
        perm=1:NDS;perm([4 NDS-1])=[NDS-1 4]; 
        S=permute(S,perm);
        S=bsxfun(@times,S,conj(ph));

        yor=resampling(yor,NMaps,2);
        for n=1:ND
            yor=ifftshift(yor,n);
            yor=ifftGPU(yor,n,gpuF);
        end    
        yor=abs(angle(bsxfun(@times,S,conj(yor))))>pi/2;
        S(yor)=-S(yor);yor=[];
    end

    if ~useFull
        S=resampling(S,NPm(1:ND));
        W=abs(resampling(W,NPm(1:ND)));
    else    
        S=resampling(S,N(1:ND));
        W=abs(resampling(W,N(1:ND)));
    
        %S=mirroring(S,rec.Alg.parE.mirr,0);
        %W=mirroring(W,rec.Alg.parE.mirr,0);
    end    
     
    if parE.virCo
        if parE.virCo==3
            S=bsxfun(@times,S,conj(pb));
            perm=[4 5 1 2 3];
            S=permute(S,perm);   
            S=matfun(@mtimes,P,S);
            S=ipermute(S,perm);
        else
            NS=size(S);
            S=dynInd(S,2:NS(4),4);
        end
    end
    
    for n=1:length(rec.Alg.parE.mirr)
        if rec.Alg.parE.mirr(n)>0
            mirr=rec.Alg.parE.mirr(n);
            S=dynInd(S,mirr+1:size(S,n)-mirr,n);
            W=dynInd(W,mirr+1:size(W,n)-mirr,n);
        end
    end
    
    if isfield(rec,'x')
        N=size(rec.x);N(end+1:3)=1;N=N(1:3);
        S=resampling(S,N,2);
        W=resampling(W,N,2);
    end
    perm=[1:4 6 5];
    S=permute(S,perm);W=permute(W,perm);
    
    %Fat-shitt
    %if rec.Par.Labels.WFS==0;SH{2}=-round(29.5507);else SH{2}=-round(rec.Par.Labels.WFS);end%Water-Fat shift  
    %Saux=shifting(dynInd(S,1,6),SH,[],1);
    %Saux=shifting(dynInd(S,1,6),SH);
    %for s=1:size(Saux,4);Saux(:,:,:,s)=extrapolating(Saux(:,:,:,s),abs(Saux(:,:,:,s))~=0,'PG-CG');end
    %S=cat(6,S,Saux);
    %W=dynInd(W,size(W,6)+1,6,0);
    
    if isempty(rec.S);rec.S=S;rec.W=W;else rec.S=cat(5,rec.S,S);rec.W=cat(5,rec.W,W);end
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing absolute phase and normalized maps: %.3f s\n',tend);end
end

indM=[7 27];
for m=indM
    if ~any(ismember(rec.Dyn.Typ2Rec,m));rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,m);end
    rec.Dyn.Typ2Wri(m)=1;
end

rec.Dyn.Debug=debugOr;