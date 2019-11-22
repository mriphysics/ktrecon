function pathIn=rawFolderDetection(path,rootIn,infer)

%RAWFOLDERDETECTION   Detects the location of a given raw data path on
%pnraw
%   PATHIN=RAWFOLDERDETECTION(PATH,{ROOTIN})
%   * PATH is the relative path to raw data
%   * {ROOTIN} is the folder where the raw data has been mounted. It
%   defaults to Data/rawSource (relative to the user's home folder)
%   * INFER tries to infer the existance of a given path by using
%   validation and date
%   ** PATHIN is the absolute path to the raw data
%

%BASIC PATH TO RAW DATA
if nargin<2;rootIn=fullfile(filesep,'home','lcg13','Data','rawSource');end
if nargin<3 || isempty(infer);infer=0;end

if infer
    pathSpl=strsplit(path,'/');
    path=pathSpl{1};
end

%SEVERAL OPTIONS FOR RAW DATA PLACEMENT IN ORDER OF PRECEDENCE
pathInV{1}=fullfile(rootIn,'pnraw','raw-nnu',path);
pathInV{2}=fullfile(rootIn,'archive-dhcp-rawdata','archive-nnu',path);
pathInV{3}=fullfile(rootIn,'archive-rawdata','archive-nnu',path);
pathInV{4}=fullfile(rootIn,'..','DISORDER/PHANTOM',path);
pathInV{5}=fullfile(rootIn,'..','DISORDER/VIVO',path);
pathInV{6}=fullfile(rootIn,path);
pathInV{7}=fullfile(rootIn,'pnraw','raw-ingenia',path);
pathInV{8}=fullfile(rootIn,'raw-nnu',path);

if ~infer
    pathIn=[];
    for n=1:length(pathInV)
        if exist(pathInV{n},'dir');pathIn=pathInV{n};break;end
    end
    assert(~isempty(pathIn),'Unable to find raw data folder %s',path);
else
    pathIn=0;
    for n=1:length(pathInV)
        if exist(pathInV{n},'dir')
            a=dir(pathInV{n});
            a={a.name};
            for m=1:length(a)
                if contains(a{m},pathSpl{2})
                    pathIn=1;
                    break
                end
            end
        end
        if pathIn;break;end
    end
end
