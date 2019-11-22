function writeGifKT(name,x,fps,ups,gr)

%WRITEGIFKT   Writes a gif with a given sampling rate and upsampling factor
%   WRITEGIFKT(NAME,X,FPS,UPS,GR)  
%   * NAME is a file name
%   * X is the data
%   ** {FPS} is the sampling rate, defaults to 20
%   ** {UPS} is an upsampling factor, defaults to 2
%   ** {GR} is a gibbs ringing smoothing factor, it defaults to 0
%

if nargin<3 || isempty(fps);fps=20;end
if nargin<4 || isempty(ups);ups=2;end
if nargin<5 || isempty(gr);gr=0;end

gpu=isa(x,'gpuArray');
NX=size(x);NX(end+1:5)=1;

if gr~=0 && gr<=1;x=filtering(x,buildFilter(NX(1:3),'tukeyIso',[],gpu,gr));end
x=resampling(x,NX(1:2)*ups);

x=abs(x);
NX=size(x);NX(end+1:5)=1;
if NX(3)<NX(4);x=permute(x,[1 3 2 4 5]);else x=permute(x,[1 4 2 3 5]);end
NX=size(x);NX(end+1:5)=1;
x=reshape(x,[NX(1)*NX(2) NX(3)*NX(4) NX(5)]);
x=x/max(x(:));

for s=1:size(x,3)
    [SIf,cm]=gray2ind(gather(x(:,:,s)),256);  
    if s==1
        imwrite(SIf,cm,name,'Delay',1/fps,'Loop',inf);
    else
        imwrite(SIf,cm,name,'Delay',1/fps,'WriteMode','append');
    end      
end