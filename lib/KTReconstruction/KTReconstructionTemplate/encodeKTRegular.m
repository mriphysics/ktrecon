function x=encodeKTRegular(x,S,F,FC,y)

%ENCODEKTREGULAR   Forward model of the reconstruction in the presence of 
%motion in case of regular shots
%   X=ENCODEKTREGULAR(X,{S},{F},{FC},{Y})  
%   * X is the image
%   * {S} is the coil array sensitivity map
%   * {F} contains the Fourier operator of the first shot
%   * {FC} contains the phase correction for the other shots
%   * {Y} contains the samples, for computing the residuals
%   ** X is the encoded data
%

if nargin<2;S=[];end
if nargin<3;F=[];end
if nargin<4;FC=[];end
if nargin<5;y=[];end

if ~isempty(FC);x=bsxfun(@times,x,FC);end
if ~isempty(S);x=bsxfun(@times,x,S);end%S
if ~isempty(F);x=matfun(@mtimes,x,F);end
if ~isempty(y);x=bsxfun(@minus,x,y);end%To compute the residuals
