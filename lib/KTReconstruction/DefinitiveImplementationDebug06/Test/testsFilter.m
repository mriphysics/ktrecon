addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug06'));


fil{1}='/home/lcg13/Data/ci_06032019_1645012_4_2_mv2mpragebostondisordernomotionsenseV4_Tr.mat';
fil{2}='/home/lcg13/Data/ci_06032019_1649339_5_2_mv2mpragebostondisordermotionsenseV4_Tr.mat';

for n=2:2    
    load(fil{n});
    rec{n}=MotionInfo;    
    NL=length(rec{n}.Residuals);%Number of levels    
    N=size(rec{n}.Residuals{NL});N=N(1:2);%Size of k-space
    itL=zeros(1,NL);
    for l=1:NL
        rec{n}.Residuals{l}=resampling(rec{n}.Residuals{l},N,1);
        rec{n}.States{l}=resampling(rec{n}.States{l},N,1);
        itL(l)=size(rec{n}.Residuals{l},4);
    end
    r{n}=cat(4,rec{n}.Residuals{:});
    r{n}=dynInd(r{n},cumsum(itL),4);
    A{n}=cat(4,rec{n}.States{:});
    S{n}=extractShots(rec{n}.ktraj{NL},rec{n}.kRange,1,1,1);
    for m=1:2
        r{n}=fftshift(r{n},m);
        A{n}=fftshift(A{n},m);        
        S{n}=fftshift(S{n},m);
    end
    r{n}=dynInd(r{n},7,4);
    A{n}=dynInd(A{n},7,4);
    
    figure
    imshow(log(r{n}(:,:)),[])
    set(gcf, 'Position', get(0,'Screensize'))  
    figure
    imshow(A{n}(:,:),[])
    set(gcf, 'Position', get(0,'Screensize'))
    
    rS{n}=bsxfun(@times,r{n},S{n});
    rS{n}=multDimSum(rS{n},1:2);
    rS{n}=repmat(rS{n},[1 1 1 1 1 6]);
    visMotion(rS{n},[],0)
    
    %Tsta=rec{n}.Motions{1};
    %for l=1:size(Tsta,1)
    %    visMotion(dynInd(Tsta,l,1),[],0)
    %end
    
    %Tsta=rec{n}.Motions{3};
    %for l=1:size(Tsta,1)
    %    visMotion(dynInd(Tsta,l,1),[],0)
    %end
end


