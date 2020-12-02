function [ktAcqWin, ktTrnWin, csmWgt, swpParams] = mrecon_ktsweep_window( ACQ, TRN, swpWinFullWidth, varargin  )
%MRECON_KTSWEEP_WINDOW  custom k-space window from k-t SWEEP acquisition
%
%   MRECON_KT( ktAcqSwp, ktTrnSwp, csmSwp, swpWinFullWidth )
%
%   MRECON_KT( ..., 'outputdir', outputDirPath )
%
%   MRECON_KT( ..., 'outputname', outFilePrefix )
%
%   MRECON_KT( ..., 'sweepwindowlocations', swpWinLoca )
%
%   MRECON_KT( ..., 'sweepwindowspacing', swpWinSpacing )
%
%   MRECON_KT( ..., 'sweepwindowindex', iWin )
%
%   MRECON_KT( ..., 'verbose', true )
%
%   Specify window location and sizes within a k-t SWEEP acquisition and return the
%   corresponding k-space and coil sensitivity data for k-t SENSE reconstruction
%

% TAR (t.roberts@kcl.ac.uk)


%% Data Dimensions

dimX = 1;  nX = size( ACQ.Data, dimX );    % frequency-encode
dimY = 2;  nY = size( ACQ.Data, dimY );    % phase-encode
dimT = 5;  nT = size( ACQ.Data, dimT );    % time / dynamics
dimC = 4;  nC = size( ACQ.Data, dimC );    % channel
dimZ = 8;  nZ = size( ACQ.Data, dimZ );    % loca / slices

nZnT     = nT * nZ;                        % no. loca in Sweep volume


%% Optional Input Argument Default Values

default.outputDirPath       = pwd;
default.outFilePrefix       = '';
default.swpWinSpacing       = [];
default.swpWinLoca          = [];
default.iWin                = [];
default.isVerbose           = false;


%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end


addRequired(    p, 'ACQ', ...
    @(x) validateattributes( x, {'MRecon'}, {'scalar'}, mfilename) );

addRequired(    p, 'TRN', ...
    @(x) validateattributes( x, {'MRecon'}, {'scalar'}, mfilename) );

addRequired(    p, 'swpWinFullWidth', ...
    @(x) validateattributes( x, {'numeric'}, {'positive','scalar'}, mfilename) );

add_param_fn(   p, 'outputdir', default.outputDirPath, ...
    @(x) validateattributes( x, {'char'}, {'nonempty','vector'}, mfilename) );

add_param_fn(   p, 'outputname', default.outFilePrefix, ...
    @(x) validateattributes( x, {'char'}, {'nonempty','vector'}, mfilename) );

add_param_fn(   p, 'sweepwindowspacing', default.swpWinSpacing, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','scalar'}, mfilename) ); % TAR

add_param_fn(   p, 'sweepwindowlocations', default.swpWinLoca, ...
    @(x) validateattributes( x, {'numeric'}, {'positive','vector'}, mfilename) ); % TAR

add_param_fn(   p, 'sweepwindowindex', default.iWin, ...
    @(x) validateattributes( x, {'numeric'}, {'positive'}, mfilename) );

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, ACQ, TRN, swpWinFullWidth, varargin{:} );

outputDirPath       = p.Results.outputdir;
outFilePrefix       = p.Results.outputname;
swpWinSpacing       = p.Results.sweepwindowspacing;
swpWinLoca          = p.Results.sweepwindowlocations;
iWin 			    = p.Results.sweepwindowindex;
isVerbose           = p.Results.verbose;


%% Anonymous Functions

% Swap dimensions to/from ReconFrame
%   ReconFrame:
%       1-2-3-4----5---6----7----8----9---10----11----12
%       x?y-z?chan?dyn-card?echo?loca?mix?extr1?extr2?aver
%   Otherwise:
%       x-y-dyn-chan-loca-z-card-echo-mix-extr1-extr12-aver
dim.x       = 1;
dim.y       = 2;
dim.z       = 3;
dim.chan    = 4;
dim.dyn     = 5;
dim.card    = 6;
dim.echo    = 7;
dim.loca    = 8;
dim.mix     = 9;
dim.extr1   = 10;
dim.extr2   = 11;
dim.aver    = 12;
ind_reconframe_to_xydcl = [ dim.x dim.y dim.dyn dim.chan dim.loca dim.z dim.card dim.echo dim.mix dim.extr1 dim.extr2 dim.aver ];
[~,ind_xydcl_to_reconframe] = sort( ind_reconframe_to_xydcl );
swap_dim_reconframe_to_xydcl = @( data ) permute( data, ind_reconframe_to_xydcl );
swap_dim_xydcl_to_reconframe = @( data ) permute( data, ind_xydcl_to_reconframe );


%% Create Matrix of Sweep Window Indices

clear swpWindows ktAcqSwpWin ktTrnSwpWin

% Window Half Width
swpWinHalfWidth = ceil( swpWinFullWidth / 2 );

% Window Locations (centre points)
if swpWinLoca
    swpWinSpacing = unique( diff( swpWinLoca ) );
elseif isempty( swpWinSpacing )
    swpWinSpacing = swpWinFullWidth; % M2D equivalent
end

if isempty( swpWinLoca )
    swpWinLoca = swpWinHalfWidth:swpWinSpacing:nZnT;
end

numSwpWindows = numel( swpWinLoca );

% Sweep Windows Array
for iW = 1:numSwpWindows
    swpWindows(:,iW) = ...
        swpWinLoca(iW)-swpWinHalfWidth+1:swpWinLoca(iW)+swpWinHalfWidth;
end

% Ensure Windows within Acquisition FOV
[~,swpWinOutOfBounds,~] = find(swpWindows > nZnT);
swpWindows( :, unique(swpWinOutOfBounds) ) = [];
numSwpWindows = size( swpWindows, 2 );


%% Visualise Sweep Windows
if isVerbose
    
    % Sweep Windows
    figure; hold on;
    plot(swpWindows','.b');
    
    % M2D Frames
    m2dWinLoca = 0:nT:nZnT;
    for iW = 1:numel(m2dWinLoca)
        plot(1:numSwpWindows, repmat(m2dWinLoca(iW),1,numSwpWindows),'k--');
    end
    
    xlabel('Sweep Window No.'); ylabel('Frame Index');
    legend('Sweep Windows','Location','NorthWest');
    axis([1 numSwpWindows 1 nZnT]);
    
    % Save
    hFig = gcf; hFig.Name = strcat( outFilePrefix, '_swp_windows' );
    saveas( hFig, [outputDirPath '/' hFig.Name, '.fig' ] );
    saveas( hFig, [outputDirPath '/' hFig.Name, '.png' ] ); clear hFig;
    
end


%% Partition Sweep Data using Custom K-Space Windows

% Partition Whole Sweep Volume if Window Not Specified
if isempty( iWin )
    iWin = 1:numSwpWindows;
end

% Window Check
if unique( diff( iWin ) ) > 1
    warning('Sweep Windows are not sequential. Errors possible ...');
end

% Update numSwpWindows
numSwpWindows = numel( iWin );


% ACQ
ktAcqSwp = reshape( squeeze( ACQ.Data ), nX, nY, nC, nT*nZ ); % squeeze as more efficient than swap_dim_reconframe and need order xycdl for reshape
ktAcqWin = single( zeros ( nX, nY, nC, numSwpWindows, swpWinFullWidth ) );
for iW = 1:numel(iWin)
    ktAcqWin(:,:,:,iW,:) = ktAcqSwp( :,:,:,swpWindows(:,iWin(iW) ));
end
ktAcqWin = permute( ktAcqWin, [1,2,5,3,4] );
clear ktAcqSwp

% TRN
ktTrnSwp = reshape( squeeze( TRN.Data ), nX, nY, nC, nT*nZ );
ktTrnWin = single( zeros ( nX, nY, nC, numSwpWindows, swpWinFullWidth ) );
for iW = 1:numel(iWin)
    ktTrnWin(:,:,:,iW,:) = ktTrnSwp( :,:,:,swpWindows(:,iWin(iW) ));
end
ktTrnWin = permute( ktTrnWin, [1,2,5,3,4] );
clear ktTrnSwp


%% View Sweep Window ACQ K-Space
if isVerbose
    
    figure;
    if numSwpWindows > 12
        numWinToPlot = 12;
    else
        numWinToPlot = numSwpWindows;
    end
    winToPlot = round(linspace(1,numSwpWindows,numWinToPlot));
    for ii = 1:numWinToPlot
        subplot(4,3,ii);
        imagesc(squeeze(abs( ...
            ktAcqWin(ceil(nX/2),:,:,ceil(nC/2),winToPlot(ii)) )), [0, 500] );
        title(['Window No. ' num2str(winToPlot(ii))]);
        colormap('gray');
    end
    
    % Save
    hFig = gcf; hFig.Name = strcat( outFilePrefix, '_swp_window_kspace' );
    saveas( hFig, [outputDirPath '/' hFig.Name, '.fig' ] );
    saveas( hFig, [outputDirPath '/' hFig.Name, '.png' ] ); clear hFig;
    
end


%% Create Weighted CSM
%	fprintf( 'Creating Weighted CSM ... \n' ),

csmSwp = swap_dim_reconframe_to_xydcl( ACQ.Parameter.Recon.Sensitivities.Sensitivity );

% CSM Original Locations
csmLoca = 1:nZnT;
csmLoca = reshape(csmLoca, [], nZ);

% CSM Identifier Matrix
csmIDs = 1:nZ;
csmIDs = repmat(csmIDs,nT,1);

% Calculate CSM Weighting Matrices
csmOverlapFractions = [];
uWgt = {}; % Note: need to be arrays to account for windows overlapping singular or multiple csmLoca
W    = {}; % Big matrix array of weights for simple .* with original CSM

for iW = iWin
    
    currWin  = swpWindows(:,iW);
    currWgt  = csmIDs(csmLoca(currWin));
    uWgt{iW} = unique(currWgt);
    
    % Overlap Between Windows and Original CSM
    for iWgt = 1:numel(uWgt{iW})
        csmOverlapFractions(iWgt,iW) = sum(currWgt == uWgt{iW}(iWgt));
        Wgts(iWgt,iW)          = csmOverlapFractions(iWgt,iW) / swpWinFullWidth;
        
        % Create Weights Matrices
        W{iW}(:,:,:,:,iWgt) = Wgts(iWgt,iW) * single( ones( [nX, nY, 1, nC] ) );
    end
    
end

% Create Weighted CSM by Linear Combination
csmWgt = single( zeros( [nX, nY, 1, nC, numSwpWindows] ) );

for iW = 1:numel(iWin)
    csmWgt(:,:,:,:,iW) = sum( W{ iWin(iW) } .* csmSwp(:,:,:,:,uWgt{ iWin(iW) }), 5 );
end


%% Collect Sweep Windowing Parameters

swpParams                 = struct;

swpParams.swpWindows      = swpWindows;
swpParams.numSwpWindows   = numSwpWindows;
swpParams.swpWinFullWidth = swpWinFullWidth;
swpParams.swpWinLoca      = swpWinLoca;
swpParams.swpWinSpacing   = swpWinSpacing;



end % ()