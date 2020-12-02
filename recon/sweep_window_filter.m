function [ xtRcn_Filtered, swpParams ] = sweep_window_filter( xtRcn, swpParams, rangeBadFrames )

%SWEEP_WINDOW_FILTER  filter k-t reconstructed Sweep data
%   
%   Filter frames at the beginning and end of k-t Sweep windows
% 
%   INPUT:
%       xtRcn - 4-D data
%       swpParams - Sweep window parameters
%
%
%   OUTPUT:
%       xtRcn_Filtered - 3-D unfolded Sweep data
%       swpParams - adds swpWindowsFilter which records frames used in
%           final unfolded 3-D Sweep data
%

%   TAR (t.roberts@kcl.ac.uk)


isVerbose = false;

%% Data Information
nX   = size( xtRcn, 1 );
nY   = size( xtRcn, 2 );
nZ   = size( xtRcn, 3 );
nT   = size( xtRcn, 4 );
nZnT = nZ*nT;


%% Determine Overlapping Frames

% Indices Window n
iWin1(:,1) = 1:swpParams.swpWinFullWidth;

% Indices Window n+1
iWin2(:,1) = 1+swpParams.swpWinStride : swpParams.swpWinStride+swpParams.swpWinFullWidth;

% Overlapping Frames
[~, iOF1, iOF2] = intersect(iWin1,iWin2);


%% View Overlapping Slices

if isVerbose
   
    iSlice1 = round( size(xtRcn,3)/2 ); iSlice2 = iSlice1 + 1;
    
    xtRcn_norm = xtRcn./max(xtRcn(:));

    xtRcn_diff = abs(squeeze(xtRcn_norm(:,:,iSlice1,iOF1))-...
                     squeeze(xtRcn_norm(:,:,iSlice2,iOF2)));

    implay_RR( [squeeze(xtRcn_norm(:,:,iSlice1,iOF1)),...
                squeeze(xtRcn_norm(:,:,iSlice2,iOF2)),...
                xtRcn_diff ],...
                'gray',[0,0.5]);
    
end


%% Assess Similarity Of Intersecting Frames

if nargin < 3
    rangeBadFrames = round( swpParams.swpWinFullWidth * (1/8) ); % e.g.: 96*1/8 * 2 = omit quarter of frames from each window 
    % rangeBadFrames = 0; % Do not perform window filtering
end
    
if isVerbose
    diffScore = squeeze( mean( mean (xtRcn_diff, 1), 2) );

    rangeSimilarFrames = 1+rangeBadFrames:(swpParams.swpWinFullWidth-swpParams.swpWinStride)-rangeBadFrames;

    figure; hold on;
    plot(diffScore,'.-k','Markersize',10);
    plot(rangeSimilarFrames,diffScore(rangeSimilarFrames),'.-r','Markersize',10);
    xlabel('Image number');
    ylabel('Difference Score (a.u.)');
    legend('All Frames','Frames After Filter','Location','NorthEast');
end


%% Filter Sweep Windows

% Remove Frames at Edges of Sweep Windows
% e.g.: remove first ten frames and last ten frames from each Sweep window
% swpWindows = swpParams.swpWindows;
% swpWindows(1:rangeBadFrames,:)         = [];
% swpWindows(end-rangeBadFrames+1:end,:) = [];

swpWindows = swpParams.swpWindows';
swpWindows(:,1:rangeBadFrames)         = [];
swpWindows(:,end-rangeBadFrames+1:end) = [];


xtRcn(:,:,:,1:rangeBadFrames)         = [];
xtRcn(:,:,:,end-rangeBadFrames+1:end) = [];

% View Intersecting Frames
if isVerbose
    iFrame = 200;
    [iDyn,iWin] = find(swpWindows'==iFrame);

    if numel(iWin) > 1
        imtar( [xtRcn_norm(:,:,iWin(1),iDyn(1)), ...
                xtRcn_norm(:,:,iWin(2),iDyn(2))], 0, 0.5 );
    elseif numel(iWin) == 1
        imtar( xtRcn_norm(:,:,iWin(1),iDyn(1)), 0, 0.5 );
    end
end


% Retain Single Instance of Each Frame
% swpWindows = reshape( swpWindows', size(swpWindows,1)*size(swpWindows,2), [] );
xtRcn      = reshape( xtRcn, nX, nY, size(swpWindows,1)*size(swpWindows,2) );

[~,iLast ,~] = unique(swpWindows,'last'); % last instance of frames
[~,iFirst,~] = unique(swpWindows,'first'); % first instance of frames

xtRcn_Filtered = xtRcn(:,:,iLast);


% View Unfolded xtRcn
if isVerbose
    implay_RR( xtRcn_Filtered );
end


%% Create Record of Filtered Frames
% nb: not the cleanest code - a little reverse-engineered
swpWindowsFilter = zeros( size(swpWindows) );
swpWindowsFilter( iLast  ) = 1;
swpWindowsPadding = zeros( size(swpWindows,1), rangeBadFrames );
swpWindowsFilter = cat(2, swpWindowsPadding, swpWindowsFilter, swpWindowsPadding );

swpParams.swpWindowsFilter = swpWindowsFilter';

% View Filters
if isVerbose
    swpWindowsFilter = zeros( size(swpWindows) );
    swpWindowsFilter( iLast  ) = 2;
    swpWindowsFilter( iFirst ) = swpWindowsFilter( iFirst ) + 1;
    swpWindowsPadding = zeros( size(swpWindows,1), rangeBadFrames );
    swpWindowsFilter = cat(2, swpWindowsPadding, swpWindowsFilter, swpWindowsPadding );
    
    imtar( swpWindowsFilter ); % not a great plot...
    colormap('hot');
    xlabel( 'Frame Number' );
    ylabel( 'Sweep Window Index' );
    title ( '0 = Not selected  //  1 = First  //  2 = Last  //  3 = Both' );
end

% sweep_window_filter(...)
end