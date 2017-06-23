function mrecon_preprocess( MR )
%MRECON_PREPROCESS   pre-process MRecon object k-space data
%
%   MRECON_PREPROCESS( MR )
%
%   MR              - ReconFrame MRecon object

% NOTE: MRecon object is passed by reference

MR.RandomPhaseCorrection;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.SortData;

end  % mrecon_preprocess(...)