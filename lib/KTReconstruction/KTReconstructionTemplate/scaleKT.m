function [yo,Ao]=scaleKT(y,R,A,B)

%SCALEKT   Generates y in Fourier domain, also generates the mask for a 
%given scale R
%   [Y,A]=SCALEKT(Y,{R},{A},{C})  
%   * Y is sampled information in image space
%   * {R} is the temporal grouping kernel size, defaults to 1
%   * {A} is a sampling mask
%   * {B} represents a calibration area
%   ** {Y} is sampled information in ky space
%   ** {A} is a sampling mask at this scale
%

if nargin<2 || isempty(R);R=1;end

gpu=isa(y,'gpuArray');
if gpu;gpuF=2;else gpuF=0;end

NY=size(y);NY(end+1:5)=1;
if nargin>=3%IF THERE'S CALIBRATION CODE SHOULD GO THROUGH THIS LOOP, SO A MUST BE INSERTED
    MB=size(A,3);
    %If there's calibration data we try to go high res in that area
   
    if nargin>=4;RKC=[R min(R,MB)];else RKC=R;end%Non-calibrated/Calibrated
    Rmax=R;
    if MB==1;Rmax=1;end
    NY=size(y);NY(end+1:8)=1;
    yo=zeros([NY(1:8) Rmax],'like',y);
    NA=size(A);NA(end+1:8)=1;
    Ao=zeros([NA(1:8) Rmax],'like',A);
    for c=1:length(RKC)
        if nargin>=4;indC=find(B==(c-1));else indC=1:NY(2);end
        if MB>1
            ya=circshift(dynInd(y,indC,2),[zeros(1,4) floor(RKC(c)/2)]);
            Aa=circshift(dynInd(A,indC,2),[zeros(1,4) floor(RKC(c)/2)]);
            for r=2:RKC(c)
                ya=dynInd(ya,r,9,circshift(dynInd(ya,r-1,9),[zeros(1,4) -1]));
                Aa=dynInd(Aa,r,9,circshift(dynInd(Aa,r-1,9),[zeros(1,4) -1]));
            end
            yo=dynInd(yo,{indC,1:RKC(c)},[2 9],ya);
            Ao=dynInd(Ao,{indC,1:RKC(c)},[2 9],Aa);                        
        else
            w=single(circconvmtx(ones(RKC(c),1),NY(5),1));
            if gpu;w=gpuArray(w);end
            ya=aplGPU(w,dynInd(y,indC,2),5);
            Aa=aplGPU(w,dynInd(A,indC,2),5);
            ya=bsxfun(@times,ya,1./(Aa+1e-6));
            Aa=single(Aa>1e-3);
            yo=dynInd(yo,indC,2,ya);
            Ao=dynInd(Ao,indC,2,Aa);
        end
    end
    if MB>1 && Rmax>MB%This helps to accelerate, it may be disabled if it introduces problems
        yo=bsxfun(@times,fold(yo,9,Rmax,MB),1./(fold(abs(dynInd(Ao,1,3)),9,Rmax,MB)+1e-6));
        Ao=bsxfun(@times,fold(Ao,9,Rmax,MB),1./(fold(abs(dynInd(Ao,1,3)),9,Rmax,MB)+1e-6));
    end
else
    w=single(circconvmtx(ones(R,1),NY(5),1));
    if gpu;w=gpuArray(w);end
    yo=aplGPU(w,y,5);
    if nargout>=2
        Ao=normm(yo,[],[1 3:4]);
        Ao=single(Ao>1e-3);%NOTE THIS IS AD-HOC, IT MAY BE NUMERICALLY UNSTABLE IF SUMMING UP OVER A DIFFERENT NUMBER OF POINTS
    end
end
