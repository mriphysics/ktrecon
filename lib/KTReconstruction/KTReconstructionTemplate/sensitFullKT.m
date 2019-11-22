function S=sensitFullKT(S,A)

%SENSITFULLKT   Computes the sensitivities for KT reconstruction
%   S=SENSITFULLKT(S,{A})  
%   * S is the sampled data in (Fourier domain for PE)
%   * {A} is the sampling pattern
%   ** S are the sensitivities
%

if nargin<2;A=[];end

gpu=isa(S,'gpuArray');
if gpu;gpuF=2;else gpuF=0;end

%S=dynInd(S,1,6);
S=mean(S,6);
NS=size(S);NS(end+1:8)=1;
if ~isempty(A)
    MB=size(A,3);
    if MB>1
        S=resPop(S,3,[],7);
        A=resPop(A,3,[],8);
        S=resPop(S,5,[MB NS(5)/MB],9:10);
        A=resPop(A,5,[MB NS(5)/MB],9:10);
        S=bsxfun(@times,conj(A),S);
        S=bsxfun(@rdivide,sum(S,10),sum(abs(A),10));
        S=sum(S,9)/MB;

        %S=bsxfun(@times,conj(A),S);
        %S=bsxfun(@rdivide,sum(S,5),sum(abs(A),5));
        S=resPop(S,7:8,[],3);
    else
        S=bsxfun(@times,conj(A),S);
        S=bsxfun(@rdivide,sum(S,5),sum(abs(A),5));
    end
else
    S=mean(S,5);
end
S=ifftGPU(S,2,gpuF);
p=1;
B=multDimSum(abs(S).^p,4)/prod(NS(4));
S=bsxfun(@rdivide,S,B.^(1/p));
%S=filtering(S,buildFilter(NS(1:2),'tukeyIso',1,gpu,1));%Spatial smoothing
