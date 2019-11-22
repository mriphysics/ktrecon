function [fileName,pathModa]=generateNIIFileName(rec)

%GENERATENIIINFORMATION   Generates the geometric information necessary to
%write NII files
%   [MS,MT]=WRITEDATA(REC)
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   ** FILENAME is the file name
%   ** PATHMODA is the path to the modality
%

pathModa=fullfile(rec.Names.pathOu,numbe2Modal(rec.Par.Mine.Modal));
if ~exist(pathModa,'dir');mkdir(pathModa);end
fileName=fullfile(pathModa,rec.Names.Name);
if isfield(rec,'addRec')%RECONSTRUCTION OF REPEATED ACQUISITIONS
    for n=1:length(rec.addRec)
        fileName=strcat(fileName,rec.addRec{n}.Names.Name);
    end
end
