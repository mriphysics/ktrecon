addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug06'));

%TIKH=[0 0.2];
%TIKH=[0.5 1];
%TIKH=0;
%TIKH=10;
%TIKH=[0.2 0.5 1 2];
%TIKH=[1 2];
TIKH=80;



%TIKH=2;
%[10 20 5 50 2];%[2 5];
RESO=1.1;%[1.0 1.1 1.2 1.3 1.4];%1.2;

%WORST SNR
path='2018_04_09/MC_2731';
fil='mc_09042018_0942279_9_2_sf7zt2mbzax0sensemi8s1V4';

%GOOD SNR THIN CORTEX
%path='2018_04_06/SC_2630';
%fil='sc_06042018_1338179_8_2_sf7zt2mbzax0sensemi8s1V4';

%GOOD SNR
%path='2018_03_02/MI_43210';
%fil='mi_02032018_0943009_8_2_sf7zt2mbzax0sensemi8s1V4';

%RELATIVELY GOOD SNR MORE MOTION ARTIFACTED
%path='2018_02_23/MI_41210';
%fil='mi_23022018_1032032_7_2_sf7zt2mbzax0sensemi8s1V4';


for o=1:length(TIKH)
    for p=1:length(RESO)
        fprintf('Case %s - TIKH-%.2f-RESO-%.2f\n',path,TIKH(o),RESO(p));
        %load('/home/lcg13/Data/rawDestin/ReconstructionsDebug06/2018_04_09/MC_2731/An-S2/mc_09042018_0942279_9_2_sf7zt2mbzax0sensemi8s1V4_SVR.mat');
        %load('/home/lcg13/Data/rawDestin/ReconstructionsDebug06/2018_04_09/MC_2731/An-S2/mc_09042018_0942279_9_2_sf7zt2mbzax0sensemi8s1V4_SVRPrev.mat');
        load(sprintf('/home/lcg13/Data/rawDestin/ReconstructionsDebug06/%s/An-S2/%s_SVR.mat',path,fil));
        load(sprintf('/home/lcg13/Data/rawDestin/ReconstructionsDebug06/%s/An-S2/%s_SVRPrev.mat',path,fil));

        svrPrev.sphcen0=svr.sphcen0;
        svrPrev.TrSp=svr.TrSp;
        svrPrev.ShSp=svr.ShSp;
        svrPrev.sphrad=svr.sphrad;
        svrPrev.TP=svr.TP;
        svrPrev.TV=svr.TV;
        svrPrev.TE=svr.TE;
        if isfield(svr,'TPHist');svrPrev.TPHist=svr.TPHist;end
        if isfield(svr,'TVHist');svrPrev.TVHist=svr.TVHist;end
        if isfield(svr,'TEHist');svrPrev.TEHist=svr.TEHist;end
        W=svr.W;
        svr=svrPrev;
        svr.ParSVR.FOVSize=144;%128;%IT WAS 160 BEFORE
        svr.ParSVR.MS=1;
        svr.ParSVR.targetResFact=RESO(p);%TARGET RESOLUTION AS A FACTOR OF THE GRID RESOLUTION, THE BIGGER THE LOWER RESOLUTION
        %assert(svr.ParSVR.targetResFact>1,'The target resolution factor should be strictly higher than 1 and it is %.3f',svr.ParSVR.targetResFact);
        %ordFracDer=cos(pi/(2*svr.ParSVR.targetResFact)).^(-2);

        svr.ParSVR.FWHM=4;
        svr.ParSVR.regFracOrd=[0 16];%[2 1 0];%REGULARIZATION FRACTIONAL ORDER
        svr.ParSVR.ti=[1 TIKH(o)];%[1 1 1];%TIKHONOV REGULARIZER
        svr.ParSVR.spti=[1 1/svr.ParSVR.targetResFact];%[1 1 1];%SPACING TIKHONOV REGULARIZER
        svr.ParSVR.GibbsRingi=0;

        svr.ParSVR.SlBefore=0;%SLICE PROFILE AT ACQUIRED RESOLUTION (0). SLICE PROFILE AT RECONTRUCTED RESOLUTION (1)

        fprintf('Setting up SVR\n');tsta=tic;     
        svr=svrSetUp(svr);
        tend=toc(tsta);fprintf('Time setting up SVR: %.3f s\n\n',tend);

        svr.W=W;

        svr=svrCG(svr,[],[],1e-4);

        svrWriteData(svr,sprintf('TIKH-%.2f-RESO-%.2f',svr.ParSVR.ti(2),svr.ParSVR.targetResFact),svr.x,[],[],'16VR_');
        svr=[];
        svrPrev=[];
    end
end


% 
% 
% fileRECON='/home/lcg13/Data/rawDestin/ReconstructionsDebug06/2018_04_09/MC_2731/An-S2/mc_09042018_0942279_9_2_sf7zt2mbzax0sensemi8s1V4_AVR';
% suff={'Test1.00'};
% 
% [y,MS,MT]=readNII(fileRECON,suff,1);
% 
% x=y{1};
% res=0.5;
% 
% N=size(x);
% Nr=round(N/res);
% 
% x=real(resampling(x,Nr));
% MS{1}=MS{1}*res;
% for n=1:3
%     MT{1}(n,n)=MT{1}(n,n)*res;
% end
% 
% fileName='/home/lcg13/Data/rawDestin/ReconstructionsDebug06/2018_04_09/MC_2731/An-S2/mc_09042018_0942279_9_2_sf7zt2mbzax0sensemi8s1V4_AVR';
% 
% suff={'Test0.80Res'};
% 
% y{1}=x;
% 
% writeNII(fileName,suff,y,MS,MT); 
% 
% return
