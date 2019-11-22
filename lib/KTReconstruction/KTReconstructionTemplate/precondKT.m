function [P,SH]=precondKT(S,M)

%PRECONDKT   Builds the preconditioner
%   P=PRECONDKT(S,{M})  
%   * S are the sensitivities
%   * {M} is a regularization for inversion
%   ** P is the preconditioner
%   ** SH are the conjugated sensitivity maps
%

if nargin<2 || isempty(M);M=0;end
P=1./bsxfun(@plus,sum(normm(S,[],4),5)/size(S,5),M);
if nargout>=2;SH=conj(S);end
