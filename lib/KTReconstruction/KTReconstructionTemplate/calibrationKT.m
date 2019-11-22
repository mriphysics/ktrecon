function  M=calibrationKT(x,sigma,smoothFactor)

%CALIBRATIONKT   Calculates the squared regularization kernels for kt
%   M=CALIBRATIONKT(X,{SIGMA},{SMOOTHFACTOR})
%   * X is an estimate of the data
%   * {SIGMA} is the noise level, defaults to 1
%   * {SMOOTHFACTOR} is the noise level, defaults to 16
%   ** M is the squared regularization kernel
%

if nargin<2 || isempty(sigma);sigma=1;end
if nargin<3 || isempty(smoothFactor);smoothFactor=16;end

gpu=isa(x,'gpuArray');
if gpu;gpuF=2;else gpuF=0;end

NX=size(x);NX(end+1:6)=1;
M=filtering(x,buildFilter(NX(1:2),'tukeyIso',1/smoothFactor,gpu,1));%Spatial smoothing to estimate the spectrum
%R=buildStandardDFTM(NX(5),1,gpu);R=R{1};
%M=abs(aplGPU(R,M,5));
M=abs(fftGPU(M,5,gpuF));
M=(M/sigma+1e-3).^(-2);


