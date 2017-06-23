function mrecon_k2i( MR )
%MRECON_K2I   transform MRecon object k-space data to image-space
%
%   MRECON_K2I( MR )
%
%   MR              - ReconFrame MRecon object

% NOTE: MRecon object is passed by reference

MR.GridData;
MR.RingingFilter;
MR.ZeroFill;
MR.K2IM;
MR.EPIPhaseCorrection;
MR.K2IP;
MR.GridderNormalization;

end  % mrecon_k2i(...)