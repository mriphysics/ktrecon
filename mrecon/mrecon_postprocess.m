function mrecon_postprocess( MR )
%MRECON_POSTPROCESS  post-process MRecon object image-space data
%
%   MRECON_POSTPROCESS( MR ) 
%
%   MR              - ReconFrame MRecon object

% NOTE: MRecon object is passed by reference

MR.SENSEUnfold;
MR.PartialFourier;
MR.ConcomitantFieldCorrection;
MR.DivideFlowSegments;
MR.CombineCoils;
MR.Average;
% MR.GeometryCorrection;
warning('MR.GeometryCorrection = switched off')
MR.RemoveOversampling;
MR.FlowPhaseCorrection;
MR.ReconTKE;
MR.ZeroFill;
MR.RotateImage;


end  % mrecon_postprocess(...)