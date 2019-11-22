function svr=svrRunSegmentation(svr,typ)

%SVRRUNSEGMENTATION   Runs a MRF segmentation
%   SVR=SVRRUNSEGMENTATION(SVR)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%

if typ==0%Quick segmentation---Fitting ellipsoids averaging
    svr.Mx=single(svr.Mx>0.5);
else%MRF segmentation---Removed for simplicity   
    K=32;
    S=128;
    if ~svr.PreTr;rRange=[10 55]/mean(svr.MSS);else rRange=[10 60]/mean(svr.MSS);end
    like=[0 -1 -1 0];
    kappa=[5 5];
    %if ~svr.PreTr;pri=[0 0 0 0 5 2];else pri=[0 0 0 0 5 5];end
    pri=[0 0 0 0 10 20];
    extK=[1 1];
    eroDilaFactor=[8 12];
    mirr=[1 1 1];
    %if svr.rec{1}.Par.Mine.Modal==5;factRad=1;else factRad=1.5;end
    factRad=1;
    N=size(svr.x);
    x=resampling(svr.x,svr.NNEncode,2);
    NE=size(x);
    M=brainSegmentation(x,ceil((NE+1)/2),rRange,K,S,multDimMea(svr.Elp(:,4:6,:),1:3)/(factRad*mean(svr.MSS)),like,pri,[],kappa,extK,[],eroDilaFactor);
    M=resampling(M,svr.NN,2);
    svr.Mx=dynInd(M,3,4);
    svr.Mx=single(abs(resampling(svr.Mx,N,[],mirr))>0.5);
end
svrWriteData(svr,sprintf('Step%02d-SegmentationMask',svr.ParSVR.Step),svr.Mx,[],[],[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;
if svr.ParSVR.Visualize;svrVisualization(svr.x,svr.Mx,'Masking',round(mean(svr.MSS(:))));end
