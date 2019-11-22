function svr=svrCNNBrainDetectionOld(svr)

%SVRCNNBRAINDETECTION   Detects the brain based on V-Net CNN training
%   SVR=SVRCNNBRAINDETECTION(SVR)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), 
%   orientations (svr.MT) and reconstruction parameters (svr.rec)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), 
%   orientations (svr.MT) and reconstruction parameters (svr.rec)
%

gpu=svr.rec{1}.Dyn.GPU;
svr.NV=length(svr.y);

load(sprintf('%s/fetalbrainregressnettcR2.mat',svr.ParSVR.pathNt));

NN=[16 10 6]*16;%Standard volume size
Res=2;%Resolution to obtain the segmentations
NNr=NN/Res;%Standard volume size at the resolution to obtain the segmentations
NNr2=[NNr(1) prod(NNr(2:3)) 1];%Size of the input to the network

for v=1:svr.NV    
    x=svr.y{v};
    if gpu;x=gpuArray(x);end
    NX=size(x);  
    
    %SAMPLE TO CNN SIZES
    x=resampling(x,NN,2);%Cropping/Enlarging
    
    xM=multDimMea(x,1:2);%This should use the abs once the bug in training is corrected
    xM=reshape(xM,[Res NN(3)/Res]);
    [~,iMax]=max(xM,[],1);
    nMax=(1:Res:NN(3))+iMax-1;     
    x=x(:,:,nMax);           
    x=resampling(x,NNr,0,ones(1,3));
    x=abs(x);    
    x=x/maxX;%Normalize to network dynamic range    
    x=min(max(x,0),1);   
    x=reshape(x,NNr2);%Reshape to input to the network            

    y=predict(fetalbrainregressnettc,im2uint16(gather(x))); 
    y(1:3)=y(1:3)+ceil((NX-NN)/2);
    z=buildGeometry(NX,'Ellipsoid',[y zeros(1,3)],gpu);
    svr.El{v}=gather(z);
    svr.Elp{v}=y;
    svr.y{v}=svr.y{v};
end
if svr.ParSVR.Visualize;svrVisualization(svr.y,svr.El,'CNN brain detection');end
svrWriteData(svr,'El',svr.El,1,[],'',0);
