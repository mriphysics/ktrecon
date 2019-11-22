pnrawPath='/home/lcg13/Data/pnrawOrRelease03/pnraw/raw-ingenia/';
pathData='../KTReconstructionTemplateData/';
datFile=strcat(pnrawPath,'2018_02_27/SA_489768/sa_27022018_1617221_20_2_pih1kt8i0bffem2dtrtraV4.raw');
refFile=strcat(pnrawPath,'2018_02_27/SA_489768/sa_27022018_1549237_1000_17_pih1senserefscanclearV4.raw');
%refFile=strcat(pnrawPath,'2018_02_27/SA_489768/sa_27022018_1617058_1000_38_pih1senserefscanclearV4.raw');
surFile=strcat(pnrawPath,'2018_02_27/SA_489768/sa_27022018_1551205_1000_20_pih1coilsurveyscanV4.raw');
%surFile=strcat(pnrawPath,'2018_02_27/SA_489768/sa_27022018_1616419_1000_35_pih1coilsurveyscanV4.raw');
addpath(genpath('../DefinitiveImplementationDebug06'));
addpath(genpath('.'));
%addpath(genpath('/home/lcg13/Work/MRecon-3.0.557'));
addpath(genpath('/home/lcg13/Work/MRecon-3.0.515'));
gpu=gpuDeviceCount;
if gpu;gpuF=2;else gpuF=0;end
if gpu;dev=gpuDevice;end
useDynamicCoils=0;%Whether to use dynamic coils
sigma=sqrt(2);%2;%1;%0.2;%sqrt(2);%"Noise level" %2 high SNR/sqrt(2) intermediate SNR/1 low SNR
smoothFactor=1;%Smoothing for spectral estimation

comp=1;%To compress the data when storing
C=[8 4 3 2 1 1];
dyn=[];%To select a series of dynamics, data needs to be reinstantiated by setting inst=1
sl=8;%To select a series of slices
%sl=[];

%ARRANGE DATA
inst=0;
if inst
    if gpu;wait(dev);end;tsta=tic;
    [y,S,A,~,R,fps,voxsiz]=arrangeKT(datFile,refFile,surFile,dyn,comp);%Calibration data seems acquired in a different time, so we don't use it although there's some support inside, last parameter serves to set a series of dynamics       
    if gpu;wait(dev);end;tend=toc(tsta);
    fprintf('Time reading RF: %.2fs\n',tend);
    if gpu;wait(dev);end;tsta=tic;
    S=gather(S);y=gather(y);A=gather(A);save(strcat(pathData,'data.mat'),'S','y','A','R','fps','voxsiz','-v7.3');
    if gpu;wait(dev);end;tend=toc(tsta);
    fprintf('Time writing MAT: %.2fs\n',tend);
else
    if gpu;wait(dev);end;tsta=tic;        
    load(strcat(pathData,'data.mat'));
    %load(strcat(pathData,'dataSlice8.mat'));
    if gpu;wait(dev);end;tend=toc(tsta);
    fprintf('Time reading MAT: %.2fs\n',tend);
end

%UNCOMPRESS DATA
if comp
    NY=size(y);NY(end+1:5)=1;
    NA=size(A);NA(end+1:5)=1;
    NR=NY./NA;NR(2)=1;
    A=repmat(A,NR);
    z=A;z(:)=0;
    z(A==1)=y;
    y=z;z=[];A=[];
end

if gpu;wait(dev);end;tsta=tic;
if ~isempty(sl);y=dynInd(y,sl,3);S=dynInd(S,sl,3);end
%y=gather(y);S=gather(S);save(strcat(pathData,'dataSlice8.mat'),'S','y','R','fps','voxsiz','-v7.3');pause

NY=size(y);NY(end+1:5)=1;
x=zeros([NY(1:3) length(C) NY(5)],'like',S);
if gpu;x=gpuArray(x);end
for z=1:NY(3)
    fprintf('Slice %d of %d\n',z,NY(3));
    yz=dynInd(y,z,3);Sz=dynInd(S,z,3);      
    if gpu;yz=gpuArray(yz);Sz=gpuArray(Sz);end

    %ESPIRIT
    B=sqrt(normm(Sz,[],4));%Body coil "reference"
    yr=scaleKT(yz,R);
    if ~useDynamicCoils
        yr=sum(yr,5)/NY(5);
        xr=solverKT(yr,Sz);
        Sz=sensitEspiritKT(ifftGPU(yr,2,gpuF),xr*sqrt(normm(yr)./normm(xr)),B,voxsiz);
    else
        xr=solverKT(yr,Sz);
        Sz=repmat(Sz,[1 1 1 1 NY(5)]);
        nor=sqrt(normm(yr,[],1:4)./normm(xr,[],1:4));    
        for s=1:NY(5);Sz=dynInd(Sz,s,5,sensitEspiritKT(ifftGPU(dynInd(yr,s,5),2,gpuF),dynInd(xr,s,5)*nor(s),B,voxsiz));end
    end

    %HIERARCHICAL SOLVER
    for c=C
        [yr,A]=scaleKT(yz,c);
        if c==C(1);M=[];else M=calibrationKT(xz,sqrt(R/c)*sigma,smoothFactor);end
        if c==1;[xz,n,~,Xz]=solverKTRegular(yr,Sz,A,[],M);else [xz,n,~,Xz]=solverKT(yr,Sz,A,[],M);end
        if c==C(1);xTz=xz;else xTz=cat(4,xTz,xz);end
        fprintf('KT %.2f: iterations %d\n',R/c,n);  
    end
    x=dynInd(x,z,3,gather(xTz));
end
if gpu;wait(dev);end;tend=toc(tsta);
fprintf('Time reconstructing: %.2fs\n',tend);

if gpu;wait(dev);end;tsta=tic;
NX=size(x);NX(1:2)=round(NX(1:2)/2);
%x=resampling(x,NX,2);
x=gather(x);
save(sprintf('%s/x.mat',pathData),'x');
writeGifKT(strcat(pathData,'video.gif'),x(:,:,1:2:end,:,:),fps);
if gpu;wait(dev);end;tend=toc(tsta);
fprintf('Time writing GIF: %.2fs\n',tend);
