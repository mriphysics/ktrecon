function frameDuration = mrecon_calc_frame_duration( MR )
%MRECON_CALC_FRAME_DURATION  calculate temporal resolution of MRecon object
%   
%   frameDuration = MRECON_CALC_FRAME_DURATION( MR ) returns the
%   temporal resolution, in seconds, of the ReconFrame object MR.
% 
%   If MR is not a dynamic scan, frameDuration = 1.

%   JFPvA (joshua.vanamerom@kcl.ac.uk)

nDyn = size( MR.Data, 5 );

if nDyn > 1,
    
    t = zeros(size(MR.Parameter.ImageInformation)); 
    
    for ii=1:numel(t), 
        t(ii) = MR.Parameter.ImageInformation(ii).DynamicScanTime; 
    end
    
    dt = diff( t, 1, 3 );
    
    frameDuration = median( dt(:) );

else
    
    frameDuration = 1;

end

end  % mrecon_calc_frame_duration(...)