function rec=solveESPIRITOld(rec)

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

debugOr=rec.Dyn.Debug;
rec.Dyn.Debug=2;

parE=rec.Alg.parE;
gpu=isa(rec.y,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

tsta=tic;
N=size(rec.y);

voxsiz=rec.Enc.AcqVoxelSize;
voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);
DeltaK=(1./(N(1:3).*voxsiz));

ND=length(parE.NC);%Number of dimensions
DeltaK=DeltaK(1:ND);
parE.NC=(1./parE.NC)./DeltaK;%Resolution of k-space points for calibration
parE.NC=min(parE.NC,N(1:ND)-1);%Never bigger than the image resolution
parE.NC=parE.NC-mod(parE.NC,2)+1;%Nearest odd integer
parE.K=(1./parE.K)./DeltaK;%Resolution of target coil profiles
parE.K=parE.K-mod(parE.K,2)+1;
parE.K=max(parE.K,3);
parE.K=min(parE.K,N(1:ND));%Never bigger than the image resolution
parE.K=parE.K-mod(parE.K,2)+1;

if rec.Dyn.Debug>=2
    fprintf('Size of calibration area:%s\n',sprintf(' %d',parE.NC));
    fprintf('Size of kernels:%s\n',sprintf(' %d',parE.K));
end

%WE USE THE MEDIAN OVER THE REPEATS TO ESTIMATE
NPE=min(length(rec.Par.Mine.pedsUn),N(5));
yf=zeros([N(1:4) NPE],'like',rec.y);
if isfield(rec,'x');xf=permute(rec.x,[1 2 3 5 4]);end
if NPE==1
    for s=1:N(4);yf=dynInd(yf,s,4,multDimMed(dynInd(rec.y,s,4),5:rec.Plan.NDims));end
    if isfield(rec,'x');xf=multDimMed(xf,5:rec.Plan.NDims);end
end
yf=gather(yf);
if isfield(rec,'x');xf=gather(xf);end
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
            for m=1:ND;x=fold(x,m,size(x,m),N(m));end
            y=cat(4,x,y);
        end
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time compressing channels: %.3f s\n',tend);end
    end
        
    tsta=tic;
    for n=1:ND
        y=fftGPU(y,n,gpuF)/sqrt(N(n));
        y=fftshift(y,n);   
    end
    if parE.absPh;y=cat(4,y,fftflip(y,1:ND));end
    y=resampling(y,parE.NC,2);
    if parE.absPh;yor=y(:,:,:,1:end/2);end
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time extracting calibration: %.3f s\n',tend);end

    tsta=tic;
    %COMPUTE CALIBRATION MATRIX
    %y=gather(y);
    %K=parE.K;
    %size(y)    
    %save('/home/lcg13/Work/DataDefinitiveImplementationDebug06/Hankel.mat','y','K');
    %1
    %pause    
    y=hankel(gather(y),parE.K);        
    if ND==2;y=permute(y,[1 2 4 3]);end%Coils to the third dimension
    NY=size(y);NY(end+1:4)=1;
    y=reshape(y,[NY(1) prod(NY(2:3)) NY(4)]);
    nc=NY(3);

    if isempty(parE.NCV);parE.NCV=nc;end
    assert(parE.NCV<=nc,'Not prepared to deal with number of eigenmaps bigger than the number of channels');
    NMaps=max(round(N(1:ND)./parE.subSp),parE.NC);
    S=zeros([nc parE.NCV prod(NMaps) NY(4)],'like',rec.y);
    W=zeros([1 parE.NCV prod(NMaps) NY(4)],'like',rec.y);   
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing Hankel matrix: %.3f s\n',tend);end
    if rec.Dyn.Debug>=2;fprintf('Size of Hankel matrix:%s\n',sprintf(' %d',size(y)));end
    
    for z=1:NY(4)%Slices
        tsta=tic;
        %PERFORM 1ST SVD AND CONVERT SINGULAR VECTORS INTO K-SPACE KERNELS
        %[~,D,V]=svd(gather(y(:,:,z)),'econ');
        %[~,D,V]=svd(y(:,:,z),'econ');
        
        %size(y(:,:,z))       
        %visReconstruction(log(y(1:2451,:,z)))        
        %%%%HEREHEREHERE---THIS IS NOT THE SAME AS IN SOLVEESPIRIT
        
        
        if size(y,1)<size(y,2);[D,~,V]=svdm(y(:,:,z));else [D,V]=svdm(y(:,:,z)');D=diagm(D);end
        %figure
        %plot(D)
        %visReconstruction(V)
        if gpu;[D,V]=parUnaFun({D,V},@gpuArray);end
        %[~,D,V]=svd(y(:,:,z),'econ');
        %if gpu;[D,V]=parUnaFun({D,V},@gpuArray);end
        if z==NY(4);y=[];end
        NV=size(V);
        V=reshape(V,[parE.K NY(3) NV(2)]);%Kernel
        D=diag(D);D=D(:);
        D=gather(D);    
        if parE.eigTh<=0;[~,nv]=screePoint(D.^2);else nv=find(D>=D(1)*parE.eigTh,1,'last');end%Automatic/tuned detection
        if rec.Dyn.Debug==2;fprintf('Threshold for SV (relative to the largest SV): %.3f\n',D(nv)/D(1));end

        %CROP KERNELS AND COMPUTE EIGEN-VALUE DECOMPOSITION IN IMAGE SPACE TO GET MAPS
        NDV=numDims(V);
        V=dynInd(V,1:nv,NDV);
        
        %ROTATE KERNEL TO ORDER BY MAXIMUM VARIANCE
        perm=[1:ND ND+2 ND+1];
        V=permute(V,perm);
        V=reshape(V,prod(parE.K)*nv,[]);
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing k-space kernels slice %d/%d: %.3f s\n',z,NY(4),tend);end        
        if rec.Dyn.Debug>=2;fprintf('Size of k-space kernels:%s\n',sprintf(' %d',size(V)));end
        
        tsta=tic;
        if size(V,1)<size(V,2);[~,~,v]=svd(V);else [~,~,v]=svd(V,'econ');end        
        V=V*v;    
        V=reshape(V,[parE.K nv nc]);
        NDV=numDims(V);
        perm=1:NDV;perm([NDV-1 NDV])=[NDV NDV-1];
        V=conj(permute(V,perm));
        
        if ND==3;V=gather(V);end%It generally won't fit the gpu
        for n=1:ND
            NMapsOr=size(V);
            V=ifftshift(V,n);
            NMapsOr(n)=NMaps(n);
            V=resampling(V,NMapsOr,1);
            if ND==3;V=ifftGPU(V,n,0);else V=ifftGPU(V,n,gpuF);end
        end    
        V=V*(prod(NMaps)/sqrt(prod(parE.K)));

        perm(1:2)=ND+(1:2);perm(3:ND+2)=1:ND;
        V=permute(V,perm);
        V=resSub(V,3:ND+2);
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing order by variance slice %d/%d: %.3f s\n',z,NY(4),tend);end        
        if rec.Dyn.Debug>=2;fprintf('Size of ordered kernels:%s\n',sprintf(' %d',size(V)));end

        tsta=tic;
        NSp=size(V,3);BlSp=10000;
        for a=1:BlSp:NSp;vA=a:min(a+BlSp-1,NSp);
            NVa=length(vA);
            Va=V(:,:,vA);
            if gpu;Va=gpuArray(Va);end
            [D,U]=svdm(Va);      
            D=D(:,1:parE.NCV,:);U=U(:,1:parE.NCV,:);
            if gpu;D=gpuArray(D);U=gpuArray(U);end
            D=real(D);
            if parE.virCo==1 && isfield(rec,'x');Uref=U(1,:,:);else Uref=sign(U(1,:,:));end           
            U=matfun(@mtimes,v,bsxfun(@rdivide,U,Uref));
            D=reshape(D,[1 parE.NCV NVa]);U=reshape(U,[nc parE.NCV NVa]);
            S(:,:,vA,z)=U;W(:,:,vA,z)=D;
        end;V=[];
        S=ipermute(S,perm);
        W=ipermute(W,perm);        
        S=reshape(S,[NMaps nc parE.NCV NY(4)]);
        W=reshape(W,[NMaps 1 parE.NCV NY(4)]);

        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing spatial maps slice %d/%d: %.3f s\n',z,NY(4),tend);end     
    end
    %S=bsxfun(@rdivide,S,dynInd(S,1,4));
    
    tsta=tic;
    if ND==2
        perm=[1:ND ND+[3 1:2]];
        S=permute(S,perm);
        W=permute(W,perm);
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

    if parE.eigSc<1 && parE.eigSc>0
        y=2*(W-parE.eigSc)./(1-parE.eigSc)-1;
        y=0.5+0.5*tanh(3*y);%"Standard" sigmoid
        S=bsxfun(@times,S,y);
        W=bsxfun(@times,W,y);
    end
    S=resampling(S,N(1:ND));
    W=abs(resampling(W,N(1:ND)));
     
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

    if isfield(rec,'x')
        N=size(rec.x);N(end+1:3)=1;N=N(1:3);
        S=resampling(S,N,2);
        W=resampling(W,N,2);
    end
    perm=[1:4 6 5];
    S=permute(S,perm);W=permute(W,perm);
    
    if isempty(rec.S);rec.S=S;rec.W=W;else rec.S=cat(5,rec.S,S);rec.W=cat(5,rec.W,W);end
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing absolute phase and normalized maps: %.3f s\n',tend);end
end

indM=[7 27];
for m=indM
    if ~any(ismember(rec.Dyn.Typ2Rec,m));rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,m);end
    rec.Dyn.Typ2Wri(m)=1;
end

rec.Dyn.Debug=debugOr;