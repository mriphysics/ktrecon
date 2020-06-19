function [ ktRcn, ktBln, noiseCov, xfTrn, xfMask, xfPri, xtTrn, xtTrnCus, xtRcn ] = recon_ktsense( ktAcq, ktTrn, csm, varargin )
%RECON_KTSENSE   k-t SENSE reconstruction of 2D dynamic MRI data
%
%   ktRcn = RECON_KTSENSE( ktAcq, ktTrn, csm )
%
%   [ ktRcn, ktBln, noiseCov ] = RECON_KTSENSE( ktAcq, ktTrn, csm )
%
%   RECON_KTSENSE( ..., 'param', val)
%
%   output:
%       ktRcn               reconstructed k-space data      (complex,   x-y-time)
%       ktBln               baseline k-space data           (complex,   x-y-time)
%       noiseCov            noise covariance matrix         (complex,   channel-channel)
%
%   input: 
%       ktAcq               k-t undersampled k-space data   (complex,   x-y-time-channel)
%       ktTrn               training k-space data           (complex,   x-y-time-channel)
%       csm                 coil/channel sensitivity maps   (complex,   x-y-1-channel)
%
%   option parameter-value pairs:
%       lambda0             k-t SENSE regularisation        (real/positive, scalar)
%       lambdaroi           regularisation inside ROI       (real/positive, scalar)
%       mask                ROI for reduced lambdaROI       (logical,   x-y)
%       noisevov            noise covariance matrix         (complex,   channel-channel)
%       knoise             noise k-space data              (complex,   sample-1-1-channel)
%       removeoversampling  remove frequency oversampling   (logical,   scalar)
%       fn_rmv_os           remove oversampling function    (anonymous function)
%       verbose             output intermediate steps       (logical,   scalar)
%
%   see also: recon_xtsense, recon_xfsense

% jfpva (joshua.vanamerom@kcl.ac.uk)

%% Notes

% TODO: change to do global regularisation using scalar ktRegStrength or
% spatially-varying regularisation using ktRegStrength map


%% Data Dimensions

dimX = 1;  nX = size( ktAcq, dimX );    % frequency-encode
dimY = 2;  nY = size( ktAcq, dimY );    % phase-encode
dimT = 3;  nT = size( ktAcq, dimT );    % time
dimC = 4;  nC = size( ktAcq, dimC );    % channel


%% Parse Input

default.lambda0             = 0.0014;
default.lambdaROI           = 0.01 * default.lambda0;
default.mask                = false( nX, nY );
default.noiseCov            = [];
default.kNoise              = [];
default.removeoversampling  = false;
default.xtTrnCus            = []; % TAR
default.makeHarmonicFilter  = false; % TAR
default.frameDuration       = []; %TAR
default.fn_rmv_os           = @( x ) x((round(size(x,dimX)/4)+1):(size(x,dimX)-round(size(x,dimX)/4)),:,:,:);
default.isVerbose           = false;

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(  p, 'ktAcq', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,nT,nC]}, mfilename ) );

addRequired(  p, 'ktTrn', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,nT,nC]}, mfilename ) );

addRequired(  p, 'csm', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[NaN,nY,1,nC]}, mfilename ) );

add_param_fn( p, 'lambda0', default.lambda0, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','scalar'}, mfilename ) );

add_param_fn( p, 'lambdaroi', default.lambdaROI, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','scalar'}, mfilename ) );
    
add_param_fn( p, 'mask', default.mask, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'2d'}, mfilename ) );

add_param_fn( p, 'noisecov', default.noiseCov, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nC,nC]}, mfilename ) );
    
add_param_fn( p, 'knoise', default.kNoise, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'size',[NaN,NaN,NaN,nC]}, mfilename ) );
    
add_param_fn( p, 'removeoversampling', default.removeoversampling, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'scalar'}, mfilename ) );

add_param_fn(  p, 'xtTrnCus', default.xtTrnCus, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,nT]}, mfilename ) );
    
add_param_fn( p, 'makeHarmonicFilter', default.makeHarmonicFilter, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'scalar'}, mfilename ) );
    
add_param_fn(  p, 'frameDuration', default.frameDuration, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive'}, mfilename ) );
    
add_param_fn( p, 'fn_rmv_os', default.fn_rmv_os, ...
        @(x) validateattributes( x, {'function_handle'}, ...
        {}, mfilename ) );

add_param_fn( p, 'verbose',     default.isVerbose, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

parse( p, ktAcq, ktTrn, csm, varargin{:} );

lambda0             = p.Results.lambda0;
lambdaROI           = p.Results.lambdaroi;
mask                = p.Results.mask;
noiseCov            = p.Results.noisecov;
kNoise              = p.Results.knoise;
isRmvOversampling   = p.Results.removeoversampling;
xtTrnCus            = p.Results.xtTrnCus;
makeHarmonicFilter  = p.Results.makeHarmonicFilter;
frameDuration       = p.Results.frameDuration;
fn_rmv_os           = p.Results.fn_rmv_os;
isVerbose           = p.Results.verbose;

assert( ~isreal( ktAcq ), 'ktAcq not complex-valued' )
assert( ~isreal( ktTrn ), 'ktTrn not complex-valued' )
assert( ~isreal( csm )  ,   'csm not complex-valued' )

if ~isempty( kNoise )
    assert( ~isreal( kNoise ), 'kNoise not complex-valued'  )
end


%% Setup

if ( isVerbose )
    fprintf( '\n%s()\n\n', mfilename );
end


%% Anon Fns

kt2xt = @( kt ) ifft2( ifftshift( ifftshift( kt, dimX ), dimY ) );

xt2kt = @( kt ) ifftshift( ifftshift( fft2( kt ), dimX ), dimY );

phase_correct = @( k ) abs(k) .* exp( sqrt(-1) * ( angle(k) + bsxfun( @times, pi/2 * ones( size(k) ), repmat( [+1;-1], size(k,1)/2, 1 ) ) ) );
 
inv_phase_correct = @( k ) abs(k) .* exp( sqrt(-1) * ( angle(k) + bsxfun( @times, pi/2 * ones( size(k) ), repmat( [-1;+1], size(k,1)/2, 1 ) ) ) );


%% Phase Correction

ktAcq = phase_correct( ktAcq );

ktTrn = phase_correct( ktTrn );


%% Sampling Pattern

ktSmp = single( sum( sum( ktAcq, dimC ), dimX ) ~= 0 );


%% Separate Baseline and Difference

ktBln = bsxfun( @rdivide, sum( ktAcq, dimT ), sum( ktSmp, dimT ) );

ktDff = ktAcq - bsxfun( @times, ktBln, ktSmp );


%% Transform to x-t Space

xtTrn = kt2xt( ktTrn );
xtAcq = kt2xt( ktAcq );
xtBln = kt2xt( ktBln );
xtDff = kt2xt( ktDff );
xtSmp = kt2xt( ktSmp ); 


%% Estimate Noise Covariance

if ( isempty( kNoise ) )    % estimate from acquired k-t undersampled data
    noiseData = reshape( xtAcq([1,nX],:,:,:), [], nC ) * sqrt( nX * nY * nT );
else                        % estimate from acquired noise data
    noiseData = reshape( kNoise, [], nC );
end

noiseCov = cov( noiseData );


%% Remove Oversampling

if ( isRmvOversampling )
    csm   = fn_rmv_os( csm ); % TAR
    xtTrn = fn_rmv_os( xtTrn );
    xtBln = fn_rmv_os( xtBln );
    xtDff = fn_rmv_os( xtDff );
    if size( mask, dimX ) > size( xtTrn, dimX )
        mask = fn_rmv_os( mask );
    end
end


%% Use Custom Training Data for Reconstruction

if ~isempty( xtTrnCus )

    % Scale xtTrnCus to match original xtTrn
    xtTrnScaleFactor = mean(abs(xtTrn(:))) ./ mean(abs(xtTrnCus(:)));
    
    xtTrnCus = xtTrnCus .* xtTrnScaleFactor;
    if ( isRmvOversampling )
        xtTrnCus = fn_rmv_os( xtTrnCus );
    end
end


%% k-t SENSE Reconstruction

safetyMargin    = 1/sqrt(lambda0);
alpha           = lambda0 / lambdaROI;

% Reconstruction with Different Regularization Methods
if ~any( mask(:) )  % if no mask
    
    % Uniform Regularization    
    [ xtRcn, xtBln, ~, ~, ~, xfRcn, ~, ~, xfTrn, xfMask, ~, ~, ~, ~, ~, xfPri ] = recon_xtsense( xtDff, xtSmp, xtTrn, csm, noiseCov, 'xtBln', xtBln, 'safetyMargin', safetyMargin, 'mask', mask, 'alpha', alpha, 'verbose', isVerbose );    
    
end
    
if any( mask(:) ) % if any part of mask true
    
    % Adaptive / Harmonic Regularization
    makeAdaptiveFilter = true;
        
    % with Standard Training Data
    if isempty( xtTrnCus )
        [ xtRcn, xtBln, ~, ~, ~, xfRcn, ~, ~, xfTrn, xfMask, ~, ~, ~, ~, ~, xfPri ] = recon_xtsense( xtDff, xtSmp, xtTrn, csm, noiseCov, 'xtBln', xtBln, 'safetyMargin', safetyMargin, 'mask', mask, 'alpha', alpha, 'makeAdaptiveFilter', makeAdaptiveFilter, 'makeHarmonicFilter', makeHarmonicFilter, 'frameDuration', frameDuration, 'verbose', true );
    end
    
    % with Custom Training Data
    if ~isempty( xtTrnCus )
        [ xtRcn, xtBln, ~, ~, ~, xfRcn, ~, ~, xfTrn, xfMask, ~, ~, ~, ~, ~, xfPri ] = recon_xtsense( xtDff, xtSmp, xtTrnCus, csm, noiseCov, 'xtBln', xtBln, 'safetyMargin', safetyMargin, 'mask', mask, 'alpha', alpha, 'makeAdaptiveFilter', makeAdaptiveFilter, 'makeHarmonicFilter', makeHarmonicFilter, 'frameDuration', frameDuration, 'verbose', true );
    end
    
end


%% Transform to k-t Space

% Zero-Pad to Match input
if isRmvOversampling,
    padSize = round(  ( [nX,nY] - size(xtRcn(:,:,1)) ) / 2 );
    xtRcn   = padarray( xtRcn, padSize, 0 );
    xtBln   = padarray( xtBln, padSize, 0 );
end

% Transform
ktRcn = xt2kt( xtRcn );
ktBln = xt2kt( xtBln );


%% Inverse Phase Correction

ktRcn = inv_phase_correct( ktRcn );
ktBln = inv_phase_correct( ktBln );


end  % recon_ktsense(...)