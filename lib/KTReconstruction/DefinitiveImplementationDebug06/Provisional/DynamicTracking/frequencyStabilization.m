function rec=frequencyStabilization(rec)

% FREQUENCYSTABILIZATION performs temporal tracking of volumes
%   REC=FREQUENCYSTABILIZATION(REC)
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), and the data 
%   information (rec.(rec.Plan.Types))
%   * REC is a reconstruction structure with information in tracked coordinates
%

if rec.Fail || rec.Par.Mine.Modal~=9;return;end%Only run if modality is fMRI

gpu=isa(rec.x,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

AdHocArray=rec.Par.Mine.AdHocArray;
assert(AdHocArray(1)~=102,'SIR reconstruction is not supported any longer'); 

NE=length(rec.Par.Labels.TE);
if rec.Par.Labels.TE(end)<rec.Par.Labels.TE(1);NE=1;end%Sometimes it appears there is a second TE while there isn't
x=abs(rec.x);
rec.x=gather(rec.x);
N=size(x);N(end+1:4)=1;
if NE>1 && numDims(rec.x)<=4;x=reshape(x,[N(1:3) NE N(4)/NE]);end
if NE==1 && numDims(rec.x)>=5;x=reshape(x,[N(1:3) prod(N(4:end))]);end

ND=numDims(x);ND=max(ND,4);
if ND==5;x=permute(x,[1 2 3 5 4]);end
N=size(x);N(end+1:4)=1;

subPyr=4;
mirr=[0 0 0 1];mirr(end+1:ND)=0;

T{2}=zeros([ones(1,3) N(4)],'like',x);
NL=size(subPyr);
for l=1:NL
    %GLOBAL SEARCH
    sR=subPyr(l)*ones(1,4);sR(5:ND)=1;
    sR(4)=32;%8
    NR=ceil(N./sR);
    sRR=N./NR;
    %if numel(x)>4e8;x=gather(x);end    
    xl=abs(resampling(x,NR,[],mirr,0));  
    if gpu;xl=gpuArray(xl);x=gpuArray(x);end
    xl=fftGPU(xl,2,gpuF);  
    kGrid=generateGrid([1 NR(2)],gpu,2*pi,ceil(([1 NR(2)]+1)/2));
    kGrid=ifftshift(kGrid{2},2);
    sub=2*(4.^(0:4));
    ran=1./(4.^(1:5));
    TR{2}=real(resampling(T{2},[1 1 1 NR(4)],[],mirr(1:4),0));
    for q=1:length(ran)
        xM=mean(shifting(xl,TR),4);
        rGrid=-ran(q)*NR(2):1/(sub(q)):NR(2)*ran(q);
        rGrid=permute(rGrid,[1 3 4 5 6 2]);
        rGrid=bsxfun(@plus,rGrid,TR{2});
        HH=exp(-1i*bsxfun(@times,rGrid,kGrid));
        xll=bsxfun(@times,xl,HH);  
        xll=bsxfun(@minus,xll,xM);
        %xll=ifftGPU(xll,2,gpuF);%To constrain to a specific ROI
        %xll=xll(:,floor(NR(2)/2)+1:end,:,:,:,:);          
        xll=abs(xll).^2;
        xll=multDimSum(xll,[1:3 5]);     
        [~,shi]=min(xll,[],6);
        TR{2}=indDim(rGrid,shi,6);        
    end
    T{2}=real(resampling(TR{2}*sRR(2),[1 1 1 N(4)],[],mirr(1:4),0));
end
T{2}=bsxfun(@minus,T{2},T{2}(1));
x=[];xl=[];xM=[];xll=[];
    
%TE
if NE==1;TE=rec.Par.Labels.TE(1);elseif NE==2;TE=diff(rec.Par.Labels.TE);else error('Distortion reversal only working for single or double echo data');end          
%ES
NPE=length(rec.Enc.kGrid{2});
if isfield(rec.Par.Mine,'ES')
    ES=rec.Par.Mine.ES;
else    
    WFS=rec.Par.Labels.WFS;
    ES=1000*(WFS/(3*42.576*3.35*(NPE+1)));%Using NPE+1 is really strange but has provided best match (not complete with ES read from GoalC parameters...   
    %echo spacing in msec = 1000 * (water-fat shift (per pixel)/(water-fat shift (in Hz) * echo train length))
    %    echo train length (etl) = EPI factor + 1
    %    water-fat-shift (Hz) = fieldstrength (T) * water-fat difference (ppm) * resonance frequency (MHz/T)
    %    water-fat difference (ppm) = 3.35 [2]
    %    resonance frequency (MHz/T) = 42.576  (for 1H; see Bernstein pg. 960)
end
fprintf('Echo spacing (ms): %.2f\n',ES);
N=size(rec.x);
ESFOV=ES*NPE/N(2);

B=convertB0Field(T{2},TE,ESFOV,'pix','rad',N(2));
if gpu;rec.x=gpuArray(rec.x);end
rec.x=bsxfun(@times,exp(-1i*B),rec.x);%Discard this component from the phase
%if numel(rec.x)>4e8;rec.x=gather(rec.x);T{2}=gather(T{2});end    
rec.w=shifting(rec.x,T);
if gpu;rec.w=gpuArray(rec.w);end
rec.x=gather(rec.x);
rec.Dyn.Typ2Rec(rec.Dyn.Typ2Rec==12)=[];

rec.F=convertB0Field(gather(T{2}),TE,ESFOV,'pix','Hz',N(2));
rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,24);
rec.Dyn.Typ2Wri(23:24)=1;

