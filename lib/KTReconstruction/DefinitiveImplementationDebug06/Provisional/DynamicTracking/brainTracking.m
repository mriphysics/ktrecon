function rec=brainTracking(rec)

% BRAINTRACKING performs temporal tracking of volumes
%   REC=BRAINTRACKING(REC)
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), and the data 
%   information (rec.(rec.Plan.Types))
%   * REC is a reconstruction structure with information in tracked coordinates
%

if rec.Fail;return;end

gpu=isa(rec.w,'gpuArray');
NDims=numDims(rec.w);NDims=min(NDims,3);
voxsiz=rec.Par.Scan.AcqVoxelSize(1,1:NDims);
voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);
distan=5*ones(1,3);%For soft masking
N=size(rec.w);N(end+1:4)=1;

if ~isempty(regexp(rec.Names.Name,'dhcp','ONCE')) || ~isempty(regexp(rec.Names.Name,'wip','ONCE')) || ~isempty(regexp(rec.Names.Name,'dev','ONCE'))%NOT A FETAL CASE-SAFER TO GO THIS ROUTE
    mirr=zeros(1,3);
    rec.b=rec.w;
    rec.w=gather(rec.w);
    rec.M=refineMask(multDimSum(abs(rec.b),4:5),rec.Alg.parS,voxsiz);
    %SOFT MASK, SIXTH COMPONENT
    Msoft=morphFourier(rec.M,distan,voxsiz,mirr,1);%SOFT MASK USED TO COMPUTE THE ROI
    Msoft(Msoft>1-1e-6)=1;Msoft(Msoft<1e-6)=0;    
    rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,[8;26]);rec.Dyn.Typ2Wri([8 26])=1;
    
    %ROI COMPUTATION
    rec.Enc.ROI=computeROI(Msoft);
    rec.Enc.ROI(3,:)=[1 N(3) N(3) N(3) 0 0];%To avoid problems later on with MB replication
    fprintf('ROI processing:\n%s',sprintf(' %d %d %d %d %d %d\n',rec.Enc.ROI'));
    
    %ROI EXTRACTION
    typ2Rec=rec.Dyn.Typ2Rec;
    for n=typ2Rec';datTyp=rec.Plan.Types{n};
        if ~ismember(n,24);rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1,1:3);end
    end
    
    %TRANSFORM
    ND=numDims(rec.b);ND=max(ND,4);
    rec.T=single(zeros([ones(1,ND-1) N(ND) 6]));    
else    
    AdHocArray=rec.Par.Mine.AdHocArray;
    assert(AdHocArray(1)~=102,'SIR reconstruction is not supported any longer'); 
    if AdHocArray(1)==101;MB=AdHocArray(4);else MB=1;end

    NE=length(rec.Par.Labels.TE);
    if rec.Par.Labels.TE(end)<rec.Par.Labels.TE(1);NE=1;end%Sometimes it appears there is a second TE while there isn't
    x=abs(rec.w);
    rec.w=gather(rec.w);
    N=size(x);N(end+1:4)=1;
    if NE>1 && numDims(rec.w)<=4;x=reshape(x,[N(1:3) NE N(4)/NE]);end
    if NE==1 && numDims(rec.w)>=5;x=reshape(x,[N(1:3) prod(N(4:end))]);end

    %if numDims(rec.u)>4;fprintf('Motion correction not defined for %d dimensions',numDims(rec.u));return;end
    ND=numDims(x);ND=max(ND,4);
    if ND==5;x=permute(x,[1 2 3 5 4]);end
    N=size(x);N(end+1:4)=1;
    mirr=[0 0 0 1];mirr(end+1:ND)=0;

    %APPROXIMATION TO THE CENTER
    x=gather(x);
    feat=multDimMed(x,4:ND);
    if gpu;feat=gpuArray(feat);end
    %[ME,sphcen0,sphrad]=ellipsoidalHoughTransform(feat,voxsiz,[10 25]);%This was before using radious in mmconversion to voxels
    [ME,sphcen0,sphrad]=ellipsoidalHoughTransform(feat,voxsiz,[20 60]);
    rec.Par.Mine.sphcen0=sphcen0;
    rec.Par.Mine.sphrad=sphrad;
    if N(4)>10
        writeNII('/home/lcg13/Work/DataDefinitiveImplementationDebug06/x',{'ME','feat'},{ME,feat});
        %1
        %pause
    end

    %TORSO MASK
    %These are the parameters in use, probably this mask could be made a bit tighter, but haven't explored that yet
    %rec.Alg.parS.nErode=6;rec.Alg.parS.nDilate=16;
    Ml=refineMask(feat,rec.Alg.parS,voxsiz);

    %TRACKING
    if gpu;x=gpuArray(x);end
    xl=x;
    for s=1:N(4)
        %xl=dynInd(xl,s,4,ellipsoidalHoughTransform(dynInd(x,s,4),voxsiz,[10 25],3,[],Ml));%This returns the features
        xl=dynInd(xl,s,4,ellipsoidalHoughTransform(dynInd(x,s,4),voxsiz,[20 60],3,[],Ml));%This returns the features
    end
    x=gather(x);
    NH=N.*(1+mirr);
    %NF=ones(1,length(N))/8;
    NF=voxsiz/32;%RECENT CHANGE
    %NF(4)=1/4;
    NF(4)=1/8;%RECENT CHANGE
    H=buildFilter(NH,'tukeyIso',NF,0,1,mirr);%It does not fit in the GPU, otherwise it should be done there
    if gpu;H=gpuArray(H);end
    xl=abs(filtering(xl,H,mirr));H=[];

    %if N(4)>10
    %    writeNII('/home/lcg13/Work/DataDefinitiveImplementationDebug06/x',{'xl'},{xl});
    %end

    xl=reshape(xl,[prod(N(1:3)) N(4)]);
    [~,indMaccum]=max(xl,[],1);xl=[];
    sphcen=ind2subV(N(1:3),gather(indMaccum));%Centers after tracking
    sphcen=bsxfun(@minus,sphcen,sphcen0);%Shifts to be applied

    %SHIFTING THE DATA
    T=cell(1,3);
    perm=1:4;perm([1 4])=[4 1];
    sphcenp=permute(sphcen,perm);
    for c=1:3
        T{c}=-dynInd(sphcenp,c,2);
        if gpu;T{c}=gpuArray(T{c});end
    end
    if gpu;x=gpuArray(x);end
    x=abs(shifting(x,T));

    %BRAIN MASK
    [x,feat]=parUnaFun({x,feat},@gather);
    feat=multDimMed(x,4:ND);
    if gpu;feat=gpuArray(feat);end
    K=32;
    S=128;
    %rRange=[5 30];
    rRange=[5 35];%RECENT CHANGE
    like=[0 -1 -1 0];
    %pri=[0 0 0 0 10 2];
    %pri=[0 0 0 0 10 1];%RECENT CHANGE
    pri=[0 0 0 0 5 1];%RECENT CHANGE
    [rec.M,par]=brainSegmentation(feat,sphcen0,rRange,K,S,sphrad,like,pri);
    rec.Par.Mine.EllipsoidParameters(1,:)=par;
    rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,8);rec.Dyn.Typ2Wri(8)=1;
    if gpu;x=gpuArray(x);end

    % %SHIFTING THE DATA TO THE CENTER OF THE ELLIPSOID
    % sphcenaux=sphcen;sphcenaux(:)=0;
    % sphcenaux=bsxfun(@plus,sphcenaux,sphcen0-par(1:3));%Shifts to be applied to x
    % T=cell(1,3);
    % perm=1:4;perm([1 4])=[4 1];
    % sphcenp=permute(sphcenaux,perm);
    % for c=1:3
    %     T{c}=-dynInd(sphcenp,c,2);
    %     if gpu;T{c}=gpuArray(T{c});end
    % end
    % if gpu;x=gpuArray(x);end
    % x=abs(shifting(x,T));
    % sphcen=bsxfun(@plus,sphcen,sphcen0-par(1:3));%Shifts to be applied in general

    %if N(4)>10
    %    writeNII('/home/lcg13/Work/DataDefinitiveImplementationDebug06/x',{'M'},{rec.M});
    %    1
    %    pause
    %end

    %MASK USED TO COMPUTE THE ROI-FIFTH COMPONENT OF THE MASK IMAGE
    dilate=10;dilate=dilate*[1 1.2 1];
    rec.M=cat(4,rec.M,morphFourier(dynInd(rec.M,3,4),dilate,voxsiz,mirr(1:3)));

    %SOFT MASK-SIXTH COMPONENT OF THE MASK IMAGE
    Msoft=morphFourier(dynInd(rec.M,5,4),distan,voxsiz,mirr(1:3),1);%SOFT MASK USED TO COMPUTE THE ROI
    Msoft(Msoft>1-1e-6)=1;Msoft(Msoft<1e-6)=0;
    rec.M=cat(4,rec.M,Msoft);

    %ROI COMPUTATION
    rec.Enc.ROI=computeROI(dynInd(rec.M,6,4));
    rec.Enc.ROI(3,:)=[1 N(3) N(3) N(3) 0 0];%To avoid problems later on with MB replication
    fprintf('ROI processing:\n%s',sprintf(' %d %d %d %d %d %d\n',rec.Enc.ROI'));
    rec.Par.Mine.sphcen0=rec.Par.Mine.sphcen0-rec.Enc.ROI(1:3,1)';
    rec.Par.Mine.EllipsoidParameters(1,1:3)=rec.Par.Mine.EllipsoidParameters(1,1:3)-(rec.Enc.ROI(1:3,1)').*voxsiz;

    %ROI EXTRACTION FOR REGISTRATION
    x=extractROI(x,rec.Enc.ROI,1,1:3);
    M=extractROI(dynInd(rec.M,3,4),rec.Enc.ROI,1:3);

    %REGISTRATION
    ND=numDims(x);ND=max(ND,4);
    rec.T=single(zeros([ones(1,ND-1) N(ND) 6]));
    NDT=numDims(rec.T);
    %lev=2;
    %parComp=[1 1 1 0 0 0];
    %[~,rec.T,~]=groupwiseVolumeRegistration(x,M,rec.T,[],lev,rec.Alg.parU.fractionOrder,parComp);x=[];
    kmax=[3 3];
    dk=[2 1];
    lev=[2 2];
    [~,rec.T,~]=integerShiftRegistration(x,M,rec.T,kmax,dk,lev,rec.Alg.parU.fractionOrder);x=[];

    %ADDING PREVIOUS LOCATION OF THE CENTERS
    perm=1:NDT;perm([1 2 ND NDT])=[ND NDT 1 2];
    sphcen=permute(sphcen,perm);
    rec.T=dynInd(rec.T,1:3,NDT,dynInd(rec.T,1:3,NDT)+sphcen);
    rec.T=dynInd(rec.T,4:6,NDT,0);

    %SHIFTING
    for c=1:3
        T{c}=-round(dynInd(rec.T,c,NDT));
        if gpu;T{c}=gpuArray(T{c});end
    end
    if gpu;rec.w=gpuArray(rec.w);end
    rec.b=shifting(rec.w,T);
    rec.w=gather(rec.w);
    rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,26);rec.Dyn.Typ2Wri(26)=1;

    %ROI EXTRACTION
    typ2Rec=rec.Dyn.Typ2Rec;
    for n=typ2Rec';datTyp=rec.Plan.Types{n};
        if ~ismember(n,24);rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1,1:3);end
    end

    %FOR TRACKING VISUALIZATION
    xGrid=generateGrid(N(1:3),gpu,N(1:3),[0 0 0],[0 0 0]);
    for n=1:3;xGrid{n}=bsxfun(@minus,xGrid{n},bsxfun(@plus,dynInd(rec.T,n,NDT),sphcen0(n)));end
    rGrid=sqrt(bsxfun(@plus,bsxfun(@plus,xGrid{1}.^2,xGrid{2}.^2),xGrid{3}.^2));
    rec.t=gather(single(rGrid<sphrad));rGrid=[];
    rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,25);rec.Dyn.Typ2Wri(25)=1;
end