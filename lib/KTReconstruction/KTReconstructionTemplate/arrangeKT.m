function [y,S,A,no,R,fps,voxsiz,ca]=arrangeKT(datFile,refFile,surFile,dyn,comp,gpu)

%ARRANGEKT   ReconFrame methods to arrange data
%   [Y,S,A,NO,CA,R,FPS,VOXSIZ,CA,K]=ARRANGEKT(DATFILE,REFFILE,SURFILE,DYN)  
%   * DATFILE is the name of the data file
%   * REFFILE is the name of the sensitivity file
%   * {SURFILE} is the name of the coil survey file
%   * {DYN} serves to read a series of dynamics
%   * {COMP} serves to compress the data when writing
%   * {GPU} maps to gpuArray
%   ** Y is the acquired data in image space
%   ** S are the sensitivities
%   ** A is a mask with the sampled data
%   ** NO are the noise samples
%   ** R is the acceleration factor
%   ** FPS are the frames per second
%   ** VOXSIZ is the effective voxel size
%   ** CA is the calibration data in image space
%

if nargin<3;surFile=[];end
if nargin<4;dyn=[];end
if nargin<5 || isempty(comp);comp=0;end
if nargin<6 || isempty(gpu);gpu=gpuDeviceCount;end

if gpu;gpuF=2;else gpuF=0;end
ca=[];
perm=[1 2 8 4 5 6 7 3 9 10 11 12];

MR{1}=MRecon(surFile);%Body coil low res
MR{2}=MRecon(surFile);%Surface coil low res
MR{3}=MRecon(refFile);%Reference coil low res
MR{4}=MRecon(datFile);%Data
if length(MR{4}.Parameter.Parameter2Read.mix)==2 && nargout>6;MR{5}=MRecon(datFile);end

for s=[1 2 3 4]%1:length(MR)
    if s==1
        MR{s}.Parameter.Parameter2Read.loca=0;
        MR{s}.Parameter.Parameter2Read.typ=1;
    elseif s==2
        MR{s}.Parameter.Parameter2Read.loca=1;
        MR{s}.Parameter.Parameter2Read.typ=1;
    elseif s==3
        MR{s}.Parameter.Parameter2Read.typ=1;
    else
         if ~isempty(dyn);MR{s}.Parameter.Parameter2Read.dyn=dyn(:);end       
         if s==4
            MR{s}.Parameter.Parameter2Read.mix=0;
         else
            MR{s}.Parameter.Parameter2Read.mix=1;
            MR{s}.Parameter.Parameter2Read.typ=1;
         end
    end            
    MR{s}.ReadData;
    MR{s}.RandomPhaseCorrection;
    %MR{s}.RemoveOversampling; %TAR for fcmr recons
    MR{s}.PDACorrection;
    MR{s}.DcOffsetCorrection;
    MR{s}.MeasPhaseCorrection;
    MR{s}.SortData;
    MR{s}.GridData;
    %MR{s}.RingingFilter;
    %MR{s}.ZeroFill;
    MR{s}.K2IM;
    MR{s}.EPIPhaseCorrection;
    MR{s}.K2IP;
    MR{s}.GridderNormalization;
    MT{s}=MR{s}.Transform('ijk','RAF',1);
    if s==1
        x=MR{s}.Data;
        if gpu;x=gpuArray(x);end
        %x=ifftshift(x,2);
    elseif s==2
        y=MR{s}.Data;
        if gpu;y=gpuArray(y);end
        %y=ifftshift(y,2);
        S=sensitEspiritKT(y,x*sqrt(normm(y)/normm(x)),[],MR{s}.Parameter.Scan.AcqVoxelSize,1);
    elseif s==3
        y=MR{s}.Data;
        if gpu;y=gpuArray(y);end
        %y=ifftshift(y,2);
        S=mapVolume(S,y,MT{2},MT{3},1);
        B=sqrt(normm(S,[],4)+1e-9);
        x=sum(conj(S).*y,4)./(B.^2);
        S=sensitEspiritKT(y,x*sqrt(normm(y)/normm(x)),B,MR{s}.Parameter.Scan.AcqVoxelSize,1);
    elseif s==4
        no=MR{s}.Data{5};
        if gpu;no=gpuArray(no);end
        y=MR{s}.Data{1};
        y=permute(y,perm);       
        for z=1:size(y,3)
            yz=dynInd(y,z,3);
            if gpu;yz=gpuArray(yz);end
            %yz=ifftshift(yz,2);
            yz=fftGPU(yz,2,gpuF);%Data in Fourier domain along PE
            yz=standardizeCoils(yz,no);
            y=dynInd(y,z,3,gather(yz));
        end
        S=dynInd(S,MR{s}.Parameter.Labels.CoilNrsPerStack{1},4);
        S=mapVolume(S,y,MT{3},MT{4},1);
        S=standardizeCoils(S,no);
                
        A=normm(y,[],[1 3:4]);
        A=single(A>1e-3);%NOTE THIS IS AD-HOC, IT MAY BE NUMERICALLY UNSTABLE IF SUMMING UP OVER A DIFFERENT NUMBER OF POINTS    
        if comp%Compress the data
            NY=size(y);NY(end+1:5)=1;
            NA=size(A);NA(end+1:5)=1;
            NR=NY./NA;
            y=y(repmat(A,NR)==1);
            NY(2)=sum(dynInd(A,1,5));
            y=reshape(y,NY);
        end

        voxsiz=MR{s}.Parameter.Scan.AcqVoxelSize;
        voxsiz(3)=voxsiz(3)+MR{s}.Parameter.Labels.SliceGaps(1);
        fps=1000/(MR{s}.Parameter.Labels.RepetitionTime(1)*sum((diff(MR{s}.Parameter.Encoding.KyRange,1,2)+1)./[MR{s}.Parameter.Scan.KtFactor;1]));
        R=MR{s}.Parameter.Scan.KtFactor;  
    elseif s==5
        ca=MR{s}.Data;
        if gpu;ca=gpuArray(ca);end
        ca=permute(ca,perm);
        %ca=ifftshift(ca,2);
        NC=size(ca);NC(end+1:3)=1;
        ca=fftGPU(ca,2,gpuF);
        ca=resampling(ca,[NC(1) diff(MR{s}.Parameter.Encoding.KyRange(2,:),1,2)+1 NC(3:end)],1);
        ca=standardizeCoils(ca,no);
    end
end
