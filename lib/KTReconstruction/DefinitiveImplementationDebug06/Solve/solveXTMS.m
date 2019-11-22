function rec=solveXTMS(rec)

%SOLVEXTMS   Performs a variable projection based aligned
%SENSE reconstruction following L Cordero-Grande, EJ Hughes, J Hutter, AN Price, 
%JV Hajnal, Three-Dimensional Motion Corrected Sensitivity Encoding 
%Reconstruction for Multi-Shot Multi-Slice MRI: Application to Neonatal Brain 
%Imaging. Magn Reson Med 79:1365-1376, 2018.
%   REC=SOLVEXTMS(REC)
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the encoding information (rec.Enc)
%   * REC is a reconstruction structure with reconstructed data and
%   surrogate information (rec.(rec.Plan.Types))
%

%PERMUTE PHASE ENCODES TO THE FIRST DIMENSIONS
typ2Rec=rec.Dyn.Typ2Rec;
perm=1:rec.Plan.NDims;perm(1:2)=[2 1];
for n=typ2Rec'
    if n~=5;datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end
end
gpu=rec.Dyn.GPU;gpuIn=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);if gpuIn;gpuFIn=2;else gpuFIn=0;end
debug=2;

%USE ONLY THE FIRST AVERAGE
y=dynInd(rec.y,1,12);

%BUILD SHOTS
NEchos=max(rec.Par.Labels.TFEfactor,rec.Par.Labels.ZReconLength);
ktraj=unique(rec.Assign.z{2}(rec.Assign.z{8}==0),'stable');
%figure
%plot(ktraj)
%pause
NProfs=numel(ktraj);
if NEchos==1;NEchos=NProfs;end
NShots=NProfs/NEchos;
if rec.Dyn.Debug>=2
    fprintf('Number of echoes: %d\n',NEchos);
    fprintf('Number of shots: %d\n',NShots);
end
ktraj=reshape(ktraj,[NEchos NShots]);
kRange=rec.Enc.kRange{2};
A=extractShots(ktraj,kRange);%subTree before

%COIL ARRAY COMPRESSION AND RECONSTRUCTED DATA GENERATION
parXT=rec.Alg.parXT;
%NX=size(rec.M);
%if prod(NX(1:3))>rec.Dyn.MaxMem(2);y=gather(y);rec.M=gather(rec.M);rec.S=gather(rec.S);end
[S,y,eivaS]=compressCoils(rec.S,parXT.perc,y);
S=dynInd(S,1:eivaS(1),4);y=dynInd(y,1:eivaS(1),4);
NS=size(S);
if gpuIn;y=gpuArray(y);end
if rec.Dyn.Debug>=2 && ~isempty(parXT.perc)
    if parXT.perc(1)<1;fprintf('Number of compressed coil elements at%s%%: %d\n',sprintf(' %0.2f',parXT.perc*100),NS(4));else fprintf('Number of compressed coil elements: %d\n',NS(4));end
end

%SLICE PROFILE PARAMETERS
parXT.SlTh=rec.Par.Scan.AcqVoxelSize(3);%Slice thickiness
parXT.SlOv=rec.Par.Scan.AcqVoxelSize(3)+rec.Par.Scan.SliceGap(1);%Slice overlap

%CREATE RECONSTRUCTION DATA ARRAY
NX=size(rec.M);NX(end+1:3)=1;
rec.d=zeros(NX,'like',y);
if ~any(rec.Dyn.Typ2Rec==16);rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,16);rec.Dyn.Typ2Wri(16)=1;end
typ2Rec=rec.Dyn.Typ2Rec;

%INITIALIZATION
NA=size(A);NX=size(rec.M);
rec.T=single(zeros([1 1 1 1 NA(5) NX(3) 6]));
outlD=single(ones([1 1 NX(3) 1 NA(5)]));
if gpu;outlD=gpuArray(outlD);end
NS=size(rec.S);NY=size(y);

%ALIGNED RECONSTRUCTION
if parXT.correct==0;estT=0;else estT=[1 0];end
L=length(estT);
for l=1:L
    if debug==2;fprintf('Resolution level: %d\n',l);tstart=tic;end
    res=2^(L-l);
    NXres(3)=NX(3);NYres(3)=NY(3);
    NXres(1:2)=floor(NX(1:2)/res);NYres(1:2)=floor(NY(1:2)/res);        
    [xRes,SRes,MRes]=parUnaFun({rec.d,S,rec.M},@resampling,NXres);    
    MRes=single(abs(MRes)>0.5);
    yRes=resampling(y,NYres);    
    ARes=resampling(A,NYres(1),1);
    cent=[1 1 1];
    %xGrid=generateGrid(NXres,gpu,NX,cent);
    xGrid=generateGridMS(NX,0,NXres,cent,0);     
    if debug==2;telapsed=toc(tstart);fprintf('Time setting multires: %.4fs\n',telapsed);end

    [xRes,rec.T,outlD]=optimizeLevelMS2D(xRes,yRes,rec.T,SRes,MRes,ARes,xGrid,parXT,debug,res,estT(l),outlD,parXT.alpha(L+1-l));
    rec.d=rec.M.*resampling(xRes,NX); 
end

%PERMUTE BACK
for n=typ2Rec'
    if n~=5;datTyp=rec.Plan.Types{n};
        rec.(datTyp)=permute(rec.(datTyp),perm);
    end
end

%POSTPROCESSING
if rec.Alg.MargosianFilter%PARTIAL FOURIER
    Enc=rec.Enc;
    rec.d=margosianFilter(rec.d,Enc);
end
rec.d=removeOverencoding(rec.d,rec.Alg.OverDec);%REMOVE OVERDECODING  
rec.Dyn.Typ2Wri(17)=1;%To write the motion estimates to file

