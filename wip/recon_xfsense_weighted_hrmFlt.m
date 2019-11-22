function [ xtRlt, xtBln, safetyMargin, xfFlt, unmix, xfRlt, xfBln, xfAcq, xfTrn, xfMask, psi, csm, mask, dt, alpha, xfPri ] = recon_xtsense( xtAcq, xtSmp, xtTrn, csm, psi, varargin )
%RECON_XTSENSE  k-t SENSE dynamic MRI reconstruction from x-t data
% 
%   [ xtRcn, xtBln ] = RECON_XTSENSE( xtAcq, xtSmp, xtTrn, csm, psi )
%
%   input:
%       xtAcq:  acquired undersampled complex image-space data;     x-y-t-c
%       xtSmp:  image-space transform of k-space sampling pattern;  x-y-t-1
%       xtTrn:  training complex image-space data;                  x-y-t-c
%       csm:    coil sensitivity matrices;                          x-y-1-c
%       psi:    noise covariance matrix;                                c-c
%   optional:
%       xtBln:          image-space baseline (DC) of acquired data; x-y-1-c
%                       default zeros
%       safetyMargin:   safety margin;                                    1
%       mask:           fetal heart ROI;                            x-y-1-1 
%                       default full FOV
%       additional mask-related inputs (alpha,beta,dt,loFreq,xfMaskKernelWidth)
%
%       xToRecon:           reconstruct subset of x values; 
%                           default all x
%       makeAdaptiveFilter: create adaptive filter; 
%                           default false
%       isVerbose:          verbose output; 
%                           default false
%
%   ouput:
%       xtRcn:  reconstructed complex images;                       x-y-t-1
%
%   See also recon_ktsense, recon_xfsense

% jfpva (joshua.vanamerom@kcl.ac.uk)  


%% Data Dimensions

dimX = 1;  nX = size( xtAcq, dimX );
dimY = 2;  nY = size( xtAcq, dimY );
dimT = 3;  nT = size( xtAcq, dimT );  dimF = dimT;  nF = nT;
dimC = 4;  nC = size( xtAcq, dimC );


%% Parse Input

default.xtBln        = complex( zeros( nX, nY ) );
default.safetyMargin = [];

default.mask         = false( nX, nY );
default.alpha        = 4;
default.beta         = 1;
default.dt           = 0.077;  % seconds
default.loFreq       = 0.5;    % Hz
default.xfMaskKernelWidth = 1;

default.xToRecon     = 1:nX;
default.makeAdaptiveFilter = false;

% TAR: Harmonics Mask Parameters
default.makeHarmonicFilter = false;
default.ktFactor     = 8;
default.hPadding     = 1;
default.hrRange      = [105, 180];
default.frameDuration = [];

default.isVerbose    = false;

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(  p, 'xtAcq', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,nT,nC]}, mfilename ) );

addRequired(  p, 'xtSmp', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[NaN,nY,nT,NaN]}, mfilename ) );
    
addRequired(  p, 'xtTrn', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,nT,NaN]}, mfilename ) );
    
addRequired(  p, 'csm', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,1,nC]}, mfilename ) );

addRequired(  p, 'psi', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nC,nC]}, mfilename ) );

add_param_fn(  p, 'xtBln', default.xtBln, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,NaN,nC]}, mfilename ) );
    
add_param_fn( p, 'safetyMargin', default.safetyMargin, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','scalar'}, mfilename ) );
   
add_param_fn( p, 'mask',          default.mask, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

add_param_fn(  p, 'alpha', default.alpha, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','nonnegative'}, mfilename ) );

add_param_fn(  p, 'beta', default.beta, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','nonnegative'}, mfilename ) );

add_param_fn(  p, 'dt', default.dt, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','positive'}, mfilename ) );

add_param_fn( p, 'loFreq', default.loFreq, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','scalar'}, mfilename ) );
    
add_param_fn( p, 'xfMaskKernelWidth', default.xfMaskKernelWidth, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive'}, mfilename ) );   

add_param_fn( p, 'xToRecon', default.xToRecon, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','vector','>=',0,'<=',nX}, mfilename ) );    

add_param_fn( p, 'makeAdaptiveFilter', default.makeAdaptiveFilter, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

add_param_fn( p, 'makeHarmonicFilter', default.makeHarmonicFilter, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );
    
add_param_fn(  p, 'ktFactor', default.ktFactor, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','nonnegative'}, mfilename ) );
    
add_param_fn(  p, 'hPadding', default.hPadding, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','nonnegative'}, mfilename ) );
    
add_param_fn(  p, 'hrRange', default.hrRange, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'size',[1 2]}, mfilename) );
    
add_param_fn(  p, 'frameDuration', default.frameDuration, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive'}, mfilename ) );
    
add_param_fn( p, 'verbose',     default.isVerbose, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

parse( p, xtAcq, xtSmp, xtTrn, csm, psi, varargin{:} );

xtBln        = p.Results.xtBln;
safetyMargin = p.Results.safetyMargin;

mask         = p.Results.mask;
alpha        = p.Results.alpha;
beta         = p.Results.beta;
loFreq       = p.Results.loFreq;
dt           = p.Results.dt;
xfMaskKernelWidth = p.Results.xfMaskKernelWidth;

xToRecon     = p.Results.xToRecon;
makeAdaptiveFilter = p.Results.makeAdaptiveFilter;

makeHarmonicFilter = p.Results.makeHarmonicFilter;
ktFactor     = p.Results.ktFactor;
hPadding     = p.Results.hPadding;
hrRange      = p.Results.hrRange;
frameDuration = p.Results.frameDuration;

isVerbose    = p.Results.verbose;


%% Validate

% dt

if ( dt > 0.2 ),
    warning( 'Large temporal resolution (dt) specified: %g seconds', dt ),
end


% xfMaskKernelWidth

if ( numel( xfMaskKernelWidth  ) == 1 ), 
    xfMaskKernelWidth = round( xfMaskKernelWidth ) * ones(1,3);
elseif ( numel( xfMaskKernelWidth ) == 2 ),
    xfMaskKernelWidth = [ reshape( round( xfMaskKernelWidth ), 1, 2 ), 1 ];
else
    warning( 'Using xfMaskKernelWidth(1:3)' )
    xfMaskKernelWidth = reshape( round( xfMaskKernelWidth(1:3) ), 1, 3 );
end


%% Setup

if ( isVerbose ),
    fprintf( '\n%s()\n\n', mfilename );
end

% Transformation Functions

xt2xf = @( xt ) fftshift( fft( xt, [], dimT ), dimT );
xf2xt = @( xf ) ifft( ifftshift( xf, dimT ), [], dimT );


%% Transform to x-f Space 

xfAcq = xt2xf( xtAcq );  

xfSmp = xt2xf( xtSmp );

xfBln = xtBln * nF;  % NOTE: xt2xf not required since xtBln has singular time component, so scaled to match xt2xf of xtBln with size nF in time/freq.

xfTrn = xt2xf( xtTrn );


%% Identify DC Frequency

[ ~, iFdc ] = max( abs( squeeze( sum( sum( sum( xfTrn, dimX ), dimY ), dimC ) ) ) );

fMax     = 1/(2*dt);        % Hz
df       = fMax/(nF/2);     % Hz
nFLoFreq = ceil( min( fMax, loFreq ) / df );
iLoFreq  = mod( iFdc + (-(nFLoFreq-1):(nFLoFreq)) - 1, nF ) + 1;


%% Combine Channels in Training and Baseline Data

if ( size( xfBln, dimC ) > 1 || size( xfTrn, dimC ) > 1 ),
    
    % Calculate SENSE unfolding matrix for R=1 (for coil combination), 
    % using same approach as x-f unaliasing below
    unmix = complex(zeros(nX,nY,1,nC));
    [cX,cY] = find( sum( sum( csm, dimC ), dimT ) ~= 0 );
    for iR = 1:numel(cX);       
        S = reshape( csm( cX(iR), cY(iR), : ), nC, 1 );
        unmix(cX(iR),cY(iR),1,:) = inv( S'*inv(psi)*S ) * S' * inv(psi);
    end
    
    % Create coil combination anonymous function
    combine_coils = @( x ) sum( bsxfun( @times, x, unmix ), dimC );  
    
    % Combine channels in xfBln
    if ( size( xfBln, dimC ) > 1 ),
        xfBln = combine_coils( xfBln );
    end

    % Combine channels in xfTrn
    if ( size( xfTrn, dimC ) > 1 ),
        xfTrn = combine_coils( xfTrn );
    end

end


%% Create x-f Estimate

xfPri = xfTrn;
    
% Remove DC if Baseline Provided

if any( xfBln(:) )  
    xfPri(:,:,iFdc,:) = 0;
end

% Create xfMask

xfMask = ones( nX, nY, nF );

if ~( makeHarmonicFilter )

    % Adaptive Filter Mask
    if any( mask(:) )  % only if any part of mask is true

        xfMask = beta * xfMask;

        % NOTE: not setting low frequecies in xfMask to 1
        % xfMask(:,:,iLoFreq) = 1;

        xfMask( repmat( mask, [1,1,nF] ) ) = alpha;

        % Smooth Transition from mask to background

        if ( any( xfMaskKernelWidth > 1 ) )

            xfKernelX = gausswin( xfMaskKernelWidth(dimX) );  % TODO: compare to Tsao2003 frequency (Hanning) filter, e.g., xfKernel1D = hann( xfMaskKernelWidth )
            xfKernelY = gausswin( xfMaskKernelWidth(dimY) ); 
            xfKernelF = gausswin( xfMaskKernelWidth(dimF) ); 

            xfKernel = bsxfun( @times, xfKernelX * xfKernelY.', reshape( xfKernelF, 1,1,[] ) );
            xfKernelNorm = xfKernel / sum( xfKernel( : ) );
            xfMask = convn( padarray( xfMask, floor(xfMaskKernelWidth/2), 'replicate' ), xfKernelNorm, 'valid' );

        end

    end
    
    % Apply x-f Mask to x-f Estimate

    xfPri = bsxfun( @times, xfPri, xfMask );
    
    
% Low Frequency & Cardiac Harmonics x-f Mask
elseif ( makeHarmonicFilter )
    
    % Initialise
    xfLowFreqMask      = ones( nX, nY, nF );
    xfFundamentalsMask = ones( nX, nY, nF );
    xfHarmonicsMask    = ones( nX, nY, nF );

    warning('off');
    xfMask.xfLowFreqMask      = xfLowFreqMask;
    xfMask.xfFundamentalsMask = xfFundamentalsMask;
    xfMask.xfHarmonicsMask    = xfHarmonicsMask;
    warning('on');
    
%     % Create Low Frequency x-f Mask using Gaussian Kernel
%     xfLowFreqKernelCoords = ( iFdc-nT/ktFactor ):( iFdc+nT/ktFactor );
%     xfLowFreqKernelWidth  = 1/( ktFactor/2 ) * numel( xfLowFreqKernelCoords );
%     gaussWindow = gausswin( nF, xfLowFreqKernelWidth )';
%     xfLowFreqMask = permute( repmat( gaussWindow, nX, 1, nY ), [1 3 2]);

    % Create Low Frequency x-f Mask using Tukey Window
    xfLowFreqKernelCoords = ( iFdc-nT/ktFactor ):( iFdc+nT/ktFactor );
    xfLowFreqKernelWidth  = 1/( ktFactor/2 ) * numel( xfLowFreqKernelCoords );
    tukeyWindow = zeros( 1,nF );
    tukeyWindow( xfLowFreqKernelCoords ) = tukeywin( nF/( ktFactor/2 )+1, 0.5);
    xfLowFreqMask = permute( repmat( tukeyWindow, nX, 1, nY ), [1 3 2]);

    
    % Create Heartbeat Harmonics x-f Mask

    % Dilate Mask
%     mask = imdilate( double(mask), strel('disk',5) ) > 0;
    
    % Find Edge Coordinates of x-t Mask
    [ mRows, mCols ] = find(mask);
    mRows = unique( mRows ); mCols = unique( mCols );
    xMin = min(mCols); xMax = max(mCols);
    yMin = min(mRows); yMax = max(mRows);

    % Estimate Heart Rate for Identifying x-f Harmonic Locations
    [ rrInterval, ~, ~, ~, ~, fPeak ] = estimate_heartrate_xf( xtTrn, frameDuration, 'roi', mask>0, 'hrRange', hrRange, 'verbose', true, 'visible', 'off', 'useUpsampleFreq', false );


    % Find Harmonic Locations

    % Initalize
    hPeaksNeg = 1;
    hPeaksPos = 1;
    hCtr = 1;
    hSeparation = floor(diff(fPeak)/2); % Peak Separation
    % floor in case peaks not integer separated
    % this will mean slightly misaligned mask
    % TODO: improved +/- 1 voxel peak shift
    
    while hPeaksNeg > 0

        % Find Peak Locations
        hPeaksNeg(hCtr) = fPeak(1) - (hCtr * hSeparation);
        hPeaksPos(hCtr) = fPeak(2) + (hCtr * hSeparation);

        if hPeaksNeg(hCtr) < 0
            hPeaksNeg(end) = [];
            hPeaksPos(end) = [];
            break;
        end

        hCtr = hCtr + 1;

    end
    
    numHarmonics = numel(hPeaksNeg);

    % Create Masks of Fundamentals and Harmonics
    fPeakNegPad = fPeak(1)-hPadding:fPeak(1)+hPadding;
    fPeakPosPad = fPeak(2)-hPadding:fPeak(2)+hPadding;
    
    fMask = zeros( nX, nF );
    fMask( yMin+hPadding:yMax+hPadding, [fPeakNegPad, fPeakPosPad] ) = 1;

    for iH = 1:numHarmonics
        hPeakNegPad(iH,:) = hPeaksNeg(iH)-hPadding:hPeaksNeg(iH)+hPadding;
        hPeakPosPad(iH,:) = hPeaksPos(iH)-hPadding:hPeaksPos(iH)+hPadding;
    
        hPeaksNegAndPos = [hPeakNegPad(iH,:), hPeakPosPad(iH,:)];
        hPeaksNegAndPos = hPeaksNegAndPos( hPeaksNegAndPos > 0 & hPeaksNegAndPos <= nF );
        
        hMask(:,:,iH) = zeros( nX, nF );
%         hMask( yMin+hPadding:yMax+hPadding, hPeaksNegAndPos, iH ) = iH + 1;
        hMask( yMin+hPadding:yMax+hPadding, hPeaksNegAndPos, iH ) = 1/(2 * 2^(iH-1));
    end
    
    % Get Mean x-f Signal in Fundamental Peaks
    xfFundamentalsMask = permute( repmat( fMask, 1, 1, nY ), [1 3 2]);
    xfFundamentalsMask( :, [1:xMin-1, xMax+1:nY], : ) = 0; % Apply mask heart region only

    fPeaksMeanSig = mean( nonzeros( abs(xfPri) .* xfFundamentalsMask ) );
    fPeaksMaxSig  =  max( nonzeros( abs(xfPri) .* xfFundamentalsMask ) );
    
    % Scale Harmonics Masks
    harmonicMultiplier = fPeaksMaxSig;
    for iH = 1:numHarmonics
        hMask(:,:,iH) = hMask(:,:,iH) .* harmonicMultiplier;
    end
    
    hPeaksMaxSig = max(hMask(:));
        
    % Collapse Harmonics Masks
    harmonicsMask = sum( cat( 3, hMask ), 3);
    fundamentalsMask = fMask;
    

    % Convolve with Gaussian
    % TODO: 2D Gaussian width currently hard-coded to 4 (= 2*ktFactor)
    gaussKernel2D = gausswin( nF, xfLowFreqKernelWidth * 8 ) * gausswin( nF, xfLowFreqKernelWidth * 8 ).';
    
    harmonicsMaskConv = convn( harmonicsMask, gaussKernel2D, 'same' );
    harmonicsMaskConv = harmonicsMaskConv ./ max(harmonicsMaskConv(:)); % Normalise
    harmonicsMaskConv = harmonicsMaskConv .* hPeaksMaxSig; 

    xfHarmonicsMask = permute( repmat( harmonicsMaskConv, 1, 1, nY ), [1 3 2]);
    xfHarmonicsMask( :, [1:xMin-1, xMax+1:nY], : ) = 0; % Apply mask heart region only
    
    fundamentalsMaskConv = convn( fundamentalsMask, gaussKernel2D, 'same' );
    fundamentalsMaskConv = ( fundamentalsMaskConv ) ./ max(fundamentalsMaskConv(:)); % Normalise

    xfFundamentalsMask = permute( repmat( fundamentalsMaskConv, 1, 1, nY ), [1 3 2]);
    xfFundamentalsMask( :, [1:xMin-1, xMax+1:nY], : ) = 0; % Apply mask heart region only

    % Regularize Fundamentals & Harmonics Masks
    xfFundamentalsMask = xfFundamentalsMask .* alpha;
    xfHarmonicsMask    = xfHarmonicsMask .* alpha;
    
    % Combine Low Freq & Harmonics Masks
%     xfMask = max( xfHarmonicsMask , xfLowFreqMask );
%     xfMask = max( xfHarmonicsMask .* alpha, xfLowFreqMask );    
    
    % Collate x-f Mask into Single Struct
    
    xfMask.xfLowFreqMask      = xfLowFreqMask;
    xfMask.xfFundamentalsMask = xfFundamentalsMask;
    xfMask.xfHarmonicsMask    = xfHarmonicsMask;

    % Apply x-f Mask to x-f Estimate

    xfLowFreqFundMask = bsxfun( @max, xfLowFreqMask, xfFundamentalsMask );
    xfPri = bsxfun( @times, xfPri, xfLowFreqFundMask );   
    xfPri = bsxfun( @max, xfPri, xfHarmonicsMask );

end


%% Zero-Pad Baseline Data in Frequency

if ( size( xfBln, dimF ) == 1 )
    xfBln = cat( 3, zeros(nX,nY,iFdc-1,1), xfBln, zeros(nX,nY,nF-iFdc,1) );
end


%% Recon

[ xfRcn, safetyMargin, xfFlt ] = recon_xfsense( xfAcq, xfPri, xfSmp, csm, psi, 'safetyMargin', safetyMargin, 'xToRecon', xToRecon, 'makeAdaptiveFilter', makeAdaptiveFilter, 'verbose', isVerbose );

if isempty( xfFlt )
    xfFlt = zeros(nX,nY,nF,nC);
end


%% Combine Baseline and Reconstructed 

xfRlt = xfBln + xfRcn;



%% Transform x-f to x-t

xtRlt = xf2xt( xfRlt );

xtBln = xfBln(:,:,iFdc) / nF;  % NOTE: xf2xt not required since xfBln is zero outside DC frequency, so extract DC frame and scale to match xf2xt with size nF in time/freq.


%% Visualise
% if ( isVerbose )
%     H = vis_xtsense( xfRlt, xfBln, xfAcq, xfTrn, xfMask, xfFlt, psi, csm, mask, dt, safetyMargin, alpha, unmix );
%     fprintf('Figures:  \n\n')
%     for iH = 1:length(H),
%         fprintf('![](figs/%s.png)  \n',H{iH}.Name)
%     end
% end


end  % recon_xtsense(...)