function [ k, ytReferenceDyn ] = kt_sliding_window( kAcq, span, method )
% KT_SLIDING_WINDOW
%
%   k = kt_sliding_window( kAcq, span, method )
%
%   k space arrays are 4-d: freq, phase, dynamic, channel  
% 
%   span measured in frames/dynamics
% 
%   method default is nearest
% 

%   jfpva (joshua.vanamerom@kcl.ac.uk)


%% Initialise

if ~exist( 'span', 'var' )
    span = size(kAcq,3);
end

if ~exist( 'method', 'var' )
    method = 'nearest';
end


%% Setup

[ nRow, nCol, nDyn, nChan ] = size( kAcq );

if isa(kAcq,'single')
    k = zeros( size(kAcq), 'single' );
else
    k = zeros( size(kAcq ), 'double' );
end


%% Sampling Pattern

ytSamplingPattern = squeeze( sum( sum( kAcq, 4 ), 1 ) ~= 0 );
ytReferenceDyn    = nan( size( ytSamplingPattern ) );


%% Sliding Window

for iDyn = 1:nDyn
  
  absDistToAcqLines = abs( bsxfun( @minus, ytSamplingPattern, (1:nDyn)-iDyn+1 ) );
  absDistToAcqLines( absDistToAcqLines > span | ~ytSamplingPattern ) = NaN;
  
  switch method
      
      case 'nearest'

          [ ~, iDynOfClosestAcqLine ] = min( absDistToAcqLines, [], 2 );

          for iCol = 1:nCol

            k(:,iCol,iDyn,:) = kAcq(:,iCol,iDynOfClosestAcqLine(iCol),:); 

          end

          ytReferenceDyn( :, iDyn ) = iDynOfClosestAcqLine;
  
      otherwise
          
          error( 'Method %s invalid', method )
          
  end
  
end


end  % kt_sliding_window(...)
