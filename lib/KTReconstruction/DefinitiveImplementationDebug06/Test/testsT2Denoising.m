clearvars
%addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug06/Provisional/HeterogeneousSVDFilter'));
addpath(genpath('/home/lcg13/git/complexSVDShrinkageDWI/complexSVDShrinkageDWI'));

pathFile='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SVR/2018_04_09/MC_2731/An-S2/mc_09042018_0942279_9_2_sf7zt2mbzax0sensemi8s1V4';

load('/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SVR/2018_04_09/MC_2731/ZZ-HE/mc_09042018_0942279_9_2_sf7zt2mbzax0sensemi8s1V4.mat');

kyRange=Par.Encoding.KyRange;

rec.Names.Name='mc_09042018_0942279_9_2_sf7zt2mbzax0sensemi8s1V4';
rec.Par=Par;
rec.Par.Mine.Proce=1;          
rec.Fail=0;rec=reconPlanner(rec);rec=reconAlgorithm(rec);rec=reconSpecific(rec);rec.Dyn.Typ2Wri(:)=0; 
recB{1}=rec;
rec=invertData(recB,1);

suff={'Aq','Re','No'};
gpu=1;
if gpu;gpuF=2;else gpuF=0;end
[x,MS,MT]=readNII(pathFile,suff,gpu);

no=dynInd(x{3},2,4);
y=x{1};
xref=x{2};
[M,N,O]=size(y);

pFF=rec.Par.Scan.HalfScanFactors(1);
H=buildHalfScanFilter(N,pFF,0);H=H(:)';
C=1./(H.^2+eps);%Noise covariance in the spectrum
Neff=sum(H~=0);

y(no==0)=0;
%visReconstruction(y,0)
%vis3DArray(y,0,6)


x=y;x(:)=0;
yi=y./(no+eps);%Noise decorrelation in space

%yi=plugNoise(yi);

yi=fftGPU(yi,2,gpuF)/sqrt(Neff);
yi=yi(:,H~=0,:);
C=C(H~=0);

yi=double(gather(yi));
C=double(gather(C));
xi=yi;xi(:)=0;amse=yi;amse(:)=0;

test=1;
if test==1%TEST 1--GLOBAL DENOISING--this has made certain sense but it introduces blurring, at least when applied globally, perhaps good to define a temporal window in z   
    yi=yi(:,:);
    C=repmat(C,[1 O]);
    xi=heterogeneousSVShrinkage(yi,[],C);
end
%TEST 3---SLICE BY SLICE
if test==2
    for o=1:size(yi,3)
        yii=yi(:,:,o);
        [xii,amsei]=heterogeneousSVShrinkage(yii.',C.');
        xii=xii.';amsei=amsei.';
        xi(:,:,o)=xii;
        amse(:,:,o)=amsei;
    end
end

xi=single(xi);
if gpu;xi=gpuArray(xi);end
xi=reshape(xi,[M Neff O]);
x(:,H~=0,:)=xi;
x=ifftGPU(x,2,gpuF)*sqrt(Neff);
x=x.*no;
visReconstruction(x,0)
%vis3DArray(x,0,6)


J=1;
bhat=1;%0.5;%0.5;%[];%0.05;
wi=[8 32];
NY=size(y);NY(end+1:3)=1;
sH=buildShearlet(NY(1:2),J,gpu,ones(1,J),'kos');
xref=y;
if gpu;[xref,no]=parUnaFun({xref,no},@gpuArray);end      
xref=shearletFilter(xref,sH,no,bhat,wi);
%visReconstruction(xref,0,[],[])
%vis3DArray(xref,0,6)

pathFile='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SVR/2018_04_09/MC_2731/An-S2/mc_09042018_0942279_9_2_sf7zt2mbzax0sensemi8s1V4Test';


xx{1}=cat(4,y,xref,x);
suffix{1}='Test';
writeNII(pathFile,suffix,xx,{MS{1}},{MT{1}})


return


return



visReconstruction(x)
visReconstruction(y)

return

for m=1:2
    y=fftGPU(y,m,gpuF);
    y=fftshift(y,m);  
end
y=log(abs(y));
visReconstruction(y)




return
%MATRIX DIMENSIONS
N=1024;
gamma=0.2;
M=floor(N*gamma);
MN=[M N];

%SIGNAL
s=[100];%0.5:0.5:10;%[100 100];%Column vector of size R, rows representing groups

NC=4;
r=100;
s=bsxfun(@times,s,exp(linspace(0,log(r^2),NC)'));
[NC,R]=size(s);
dimHetero=1;
X=[];

U=randU([M R],1);
V=randU([R N],1);
if dimHetero==2
    V=reshape(V,[R N/NC NC]);
    for n=1:NC;X=cat(2,X,matfun(@mtimes,bsxfun(@times,U,s(n,:)),V(:,:,n)));end
else
    U=reshape(U,[M/NC NC R]);
    for n=1:NC;X=cat(1,X,matfun(@mtimes,bsxfun(@times,squeeze(U(:,n,:)),s(n,:)),V));end
end

%THIS WOULD BUILD AN OVERLAPPED WINDOW, IT HAS BEEN USELESS
%E=16;%2*N/NC;
%In=tukeywin(E,0);
%In=circconvmtx(In,N);
%subR=4;
%In=In(:,1:subR:end)';
%In=bsxfun(@times,In,1./sum(In,1));

%NOISE
k=10;
C=1*linspace(1,k,M)';
W=bsxfun(@times,sqrt(C),plugNoise(X,1));

%MATRIX AND PRECOMPUTATIONS
Y=X+W;
C=C*N;
Y=Y*sqrt(N);
Y=double(Y);
X=X*sqrt(N);
Xmax=max(X(:))^2;

%TEST 1-WE DECORRELATE, ESTIMATE AND CORRELATE
Ys=bsxfun(@times,Y,1./sqrt(C));
[amses,~,Xhats]=SVShrinkage(Ys);
Xhats=bsxfun(@times,Xhats,sqrt(C));
amses=amses.*sqrt(C);
fprintf('MSE (dB): %.2f\n',10*log10(Xmax/mean(abs(Xhats(:)-X(:)).^2)));
%fprintf('AMSE (dB): %.2f\n',10*log10(Xmax/mean(amses(:))));

%TEST 2-WE APPLY OPTSHRINK
%esd=ESDMixAndMix(C,N,[],2);
esd=ESDSimulated(C,N);
[amse,~,Xhat]=SVShrinkage(Y,[],[],esd);
fprintf('MSE (dB): %.2f\n',10*log10(Xmax/mean(abs(Xhat(:)-X(:)).^2)));
fprintf('AMSE (dB): %.2f\n',10*log10(Xmax/mean(amse(:))));

%TEST 3-WE APPLY STANDARD WEIGHTED ESTIMATION
[Xhath,amseh]=heterogeneousSVShrinkage(Y,C);
fprintf('MSE (dB): %.2f\n',10*log10(Xmax/mean(abs(Xhath(:)-X(:)).^2)));
fprintf('AMSE (dB): %.2f\n',10*log10(Xmax/mean(amseh(:))));

% %TEST 4-WE APPLY ESTIMATION OF THE NUMBER OF REGIONS
% if dimHetero==2
%     [Xhatl,amsel]=heterogeneousSVShrinkage(Y,C,[],[],0);
% else
%     [Xhatl,amsel]=heterogeneousSVShrinkage(Y,C,[],0);
% end
% fprintf('MSE (dB): %.2f\n',10*log10(Xmax/mean(abs(Xhatl(:)-X(:)).^2)));
% fprintf('AMSE (dB): %.2f\n',10*log10(Xmax/mean(amsel(:))));
% 
% %TEST 4-WE APPLY ESTIMATION OF THE NUMBER OF REGIONS WITH RECOMPUTATION OF
% %SVD
% if dimHetero==2
%     [Xhatlr,amselr]=heterogeneousSVShrinkage(Y,C,[],[],0,1);
% else
%     [Xhatlr,amselr]=heterogeneousSVShrinkage(Y,C,[],0,[],1);
% end
% fprintf('MSE (dB): %.2f\n',10*log10(Xmax/mean(abs(Xhatlr(:)-X(:)).^2)));
% fprintf('AMSE (dB): %.2f\n',10*log10(Xmax/mean(amselr(:))));

%TEST 5-WE APPLY LOCALIZED ESTIMATION WITH PERFECT LOCALIZATION OF SUBDIVISIONS
NS=NC;
InBase=ones(1,MN(dimHetero)/NS);
In=[];
for n=1:NS;In=blkdiag(In,InBase);end
if dimHetero==1;[Xhatll,amsell]=heterogeneousSVShrinkage(Y,C,[],In');
else [Xhatll,amsell]=heterogeneousSVShrinkage(Y,C,[],[],In);
end
fprintf('MSE (dB): %.2f\n',10*log10(Xmax/mean(abs(Xhatll(:)-X(:)).^2)));
fprintf('AMSE (dB): %.2f\n',10*log10(Xmax/mean(amsell(:))));

% %TEST 5-WE APPLY LOCALIZED ESTIMATION WITH MISSED LOCALIZATION OF SUBDIVISIONS
% InBase=ones(1,N/NC);
% In=[];
% for n=1:NC;In=blkdiag(In,InBase);end;
% In=circshift(In,[0 N/(2*NC)]);
% [Xhatllw,amsellw]=heterogeneousSVShrinkage(Y,C,[],[],In);
% fprintf('MSE (dB): %.2f\n',10*log10(Xmax/mean(abs(Xhatllw(:)-X(:)).^2)));
% fprintf('AMSE (dB): %.2f\n',10*log10(Xmax/mean(amsellw(:))));


% 
% figure
% imshow(Y,[])
% figure
% imshow(X,[])
% figure
% imshow(Xhatl,[])
% figure
% imshow(Xhatll,[])