function [A,echoLim]=extractShots(ktraj,kRange,NShots,NWithin)

%EXTRACTSHOTS generates the shot sampling structure for a given acquisition
%   [A,ECHOLIM]=EXTRACTSHOTS(KTRAJ,KRANGE,NSHOTS,NSUBDIV)
%   * KTRAJ is the spectral trajectory (NEchosxNShotsx2 matrix)
%   * KRANGE is the spectral coverage (2x2 matrix: Rows: PE dimensions. Columns: min/max range)
%   * NSHOTS is the number of shots (or sweeps)
%   * {NWITHIN} is the number of subdivisions per shot
%   ** A is a binary matrix with the shot structure in the fifth dimension
%   ** ECHOLIM are the echo limits in within shot subdivisions
%

if nargin<4 || isempty(NWithin);NWithin=ones(1,NShots);end

kShift=floor((diff(kRange,1,2)+1)/2);
kSize=diff(kRange,1,2)+1;

A=single(zeros([kSize' 1 1 sum(NWithin)]));

NProf=size(ktraj,1);%Number of profiles
%%%HEREHEREHERE

%GROUPING OF SAMPLES
for ns=1:NShots
    for ws=1:sub
        if ws==1%Right-hand side, decreasing from low to high frequencies
            contwt=subTree;
            for wt=maxSubTree:-1:1
                for ne=echoLim(1,(ws-1)*maxSubTree+wt):echoLim(2,(ws-1)*maxSubTree+wt)  
                    if NDimens==2
                        if ktraj(ne,ns,1)~=inf;A(ktraj(ne,ns,1)+kShift(1)+1,ktraj(ne,ns,2)+kShift(2)+1,contwt,ws,ns)=1;end
                    else
                        if ktraj(ne,ns,1)~=inf;A(ktraj(ne,ns,1)+kShift(1)+1,1,contwt,ws,ns)=1;end
                    end
                end
                if maxSubTree-wt+1<subTree;contwt=contwt-1;end
            end
        else%Left-hand side, increasing from low to high frequencies
            contwt=1;
            for wt=1:maxSubTree
                for ne=echoLim(1,(ws-1)*maxSubTree+wt):echoLim(2,(ws-1)*maxSubTree+wt)         
                    if NDimens==2
                        if ktraj(ne,ns,1)~=inf;A(ktraj(ne,ns,1)+kShift(1)+1,ktraj(ne,ns,2)+kShift(2)+1,contwt,ws,ns)=1;end
                    else
                        if ktraj(ne,ns,1)~=inf;A(ktraj(ne,ns,1)+kShift(1)+1,1,contwt,ws,ns)=1;end
                    end
                end
                if wt<subTree;contwt=contwt+1;end
            end
        end            
    end
end
for n=1:2;A=ifftshift(A,n);end
