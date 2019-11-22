addpath(genpath('.'))
addpath(genpath('DefinitiveImplementationDebug06'))

pathData='/home/lcg13/Work/GiulioKTData';
fileData='data.mat';

tic
load(strcat(pathData,filesep,fileData));
toc

gpu=gpuDeviceCount;
if gpu;gpuF=2;else gpuF=0;end
if gpu;[y,M,Akt1,Akt2,AMB2,AMB3,S]=parUnaFun({y,M,Akt1,Akt2,AMB2,AMB3,S},@gpuArray);end

%TO USE SAME SENSITIVITIES FOR ALL CASES 1, TO ESTIMATE THEM FROM THE ACCELERATED DATA, 0
useSameSensit=0;

%REVERSE ENGINEER TO ORIGIN OF SPECTRUM
if useSameSensit~=2;y=shifting(y,{17},1);end

%ESTIMATE SENSITIVITIES
if useSameSensit==1;S=sensitFullKT(fftGPU(y,2,gpuF));elseif useSameSensit==2;S=sensitFullKT(fftGPU(S,2,gpuF));end%Mine vs Giulio

%SLICE EXTRACTION
slV=1:6;
[S,y,M]=parUnaFun({S,y,M},@dynInd,slV,3);
NY=size(y);NY(end+1:6)=1;

%ROI COMPUTATION
ROI=computeROI(M);
M=extractROI(M,ROI,1,1);%To accelerate for testing

%BUILD ENCODINGS AND RECONSTRUCTION SETTINGS
%Akt1(:,2:end,:,:,:)=Akt1(:,end:-1:2,:,:,:);%SLIGHTLY PROBLEMATIC AS IT LEAVES ASYMMETRIC SPECTRUM
A{1}=Akt1;A{1}(:)=1;B{1}=zeros(1,NY(2),'like',A{1});
A{2}=AMB2;B{2}=zeros(1,NY(2),'like',A{2});
A{3}=AMB3;B{3}=zeros(1,NY(2),'like',A{3});
A{4}=Akt1;B{4}=single(sum(Akt1,5)==NY(5));
A{5}=Akt2;B{5}=single(sum(Akt2,5)==NY(5));
A{6}=AMB2.*Akt1;B{6}=single(sum(Akt1,5)==NY(5));
A{7}=AMB3.*Akt2;B{7}=single(sum(Akt2,5)==NY(5));
RV=[1 2 3 5 5 10 15];
RI=[1 2 3 5 5 5 5];
RR=[1 4 6 7];
%sigma=2e-4*ones(1,7);
if useSameSensit;sigma=2e-5*ones(1,7);else sigma=5e-6*ones(1,7);end
sigma=5e-5*ones(1,7);%Regularization weight %%%PARAMETER (try 1e-4 for more regularized solutions)
smoothFactor=16;%To smooth in space to estimate the Power Spectral Density of x-f%%%PARAMETER
gF=ones(1,7);
sigma=sigma.*gF;
leg={'GT','MB2','MB3','KT5a','KT5b','MB2KT5a','MB3KT5b'};
tic
for r=RR
    %SYNTHESIS
    Ar=A{r};Br=B{r};
    yr=encodeKT(y,[],Ar);
    
    %PARAMETERES
    R=RV(r);%Acceleration rate
    C=[RI(r):-1:1 1];%Hierarchical plan

    %ESTIMATE SENSITIVITIES
    if useSameSensit;Sr=S;else Sr=sensitFullKT(yr,Ar);end
    
    %ROI EXTRACTION
    [Sr,yr]=parUnaFun({Sr,yr},@extractROI,ROI,1,1);%To accelerate for testing, does not give exact same results due to boundary conditions when filtering to compute calibration

    %HIERARCHICAL KT
    for c=C
        [yrc,Arc]=scaleKT(yr,c,Ar,Br);
        if c==C(1);P=[];else P=calibrationKT(xr{r},sqrt(R/c)*sigma(r),smoothFactor);end%Second output should be for ktPCA
        [xr{r},n]=solverKT(yrc,Sr,Arc,[],P,[],[],1e-7);
        fprintf('KT %.2f: iterations %d\n',R/c,n);          
    end    
end
toc
plotVelocities(xr,M,leg(RR));
x=cat(1,xr{:});
x=extractROI(x,ROI,1,2);%To visualize LV ROI
visReconstruction(resPop(x,2:3,[],2),0)
