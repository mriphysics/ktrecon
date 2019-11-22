%MRECON_SC_SLW  Self-calibrated sliding window k-t SENSE reconstruction using ReconFrame
%
%   User input:
%       edit the Setup cell with details of patient to reconstruct
%
%   NOTE: currently designed to be used on beastie01 and gpubeastie03-pc
%   because MRecon is only licensed on beastie01. This script incorporates
%   code by jfpva and lcg. Requires MRecon 515 for lcg code.
%
%   see also: mrecon_kt.m, recon_kt_sc_slw.m, kt_sc_slw_preproc.m,
%   arrangeKT_wrapper.m, solverKT.m

% tar (t.roberts@kcl.ac.uk)


%% Setup
% e.g.:
% fcmrNum = 191;
% ingeniaPatDirPath = 'YYYY_MM_DD/ID_123456';
% seriesNos = [16 17 18 19 20];

fcmrNum = [];
ingeniaPatDirPath = '/';
seriesNos = [];


%% Host / Remote Computers

b01Path =     '/scratch/tr17/data/fcmr_4dflow';
g03Path = '/fastscratch/tr17/data/fcmr_4dflow';

[ ~, computerName ] = system('hostname');

if strcmp( computerName(1:end-1), 'beastie01' )
    hostComputerStudyPath   = b01Path;
    remoteComputerStudyPath = g03Path;
elseif strcmp( computerName(1:end-1), 'gpubeastie03-pc' )
    hostComputerStudyPath   = g03Path;
    remoteComputerStudyPath = b01Path;
else
    error('Unknown computer hostname.');
end
    

%% Directories

% Local
reconDir      = fullfile( hostComputerStudyPath, strcat( 'fcmr', num2str(fcmrNum) ) ) ;
outputDirPath = fullfile( reconDir, 'ktrecon');
maskDirPath   = fullfile( reconDir, 'mask');

% Remote
ktreconDirRemote = fullfile( remoteComputerStudyPath, strcat( 'fcmr', num2str(fcmrNum) ), 'ktrecon' ) ;
maskDirRemote    = fullfile( remoteComputerStudyPath, strcat( 'fcmr', num2str(fcmrNum) ), 'mask' ) ;

% Raw files
ingeniaRawDirPath = '/pnraw01/raw-ingenia';
rawDataDirPath    = fullfile( ingeniaRawDirPath, ingeniaPatDirPath);


%% Make Directories

mkdir(outputDirPath);
mkdir(maskDirPath);

% Create Folders on external computer for processing
if strcmp( computerName(1:end-1), 'beastie01' )
    eval(['!ssh tr17@gpubeastie03-pc "mkdir -p ' ktreconDirRemote ' && mkdir -p ' maskDirRemote '"' ]);
elseif strcmp( computerName(1:end-1), 'gpubeastie03-pc' )
    eval(['!ssh tr17@beastie01 "mkdir -p ' ktreconDirRemote ' && mkdir -p ' maskDirRemote '"' ]);
else
    error('Unknown computer hostname.');
end

fprintf('\nFolders created on Remote Computer.\n');


%% Dependencies

addpath( genpath('/home/tr17/MATLAB/General Useful Scripts') );
addpath( genpath('/home/tr17/MATLAB/ktrecon-dev') );


%% MRecon: need version 515 for arrangeKT

if strcmp( computerName(1:end-1), 'beastie01' )
    warning('off');
    rmpath ( genpath('/home/tr17/MATLAB/MRecon-3.0.557') );
    addpath( genpath('/home/tr17/MATLAB/MRecon-3.0.515') );
    warning('on');
else
    error('Must perform arrangeKT on computer with MRecon 515.');
end


%% Load Data with arrangeKT

if strcmp( computerName(1:end-1), 'beastie01' )

    for seriesNo = seriesNos

        [ rawDataFilePath, coilSurveyFilePath, senseRefFilePath ] = ...
            id_pnraw_data( rawDataDirPath, seriesNo );

        outFilePrefix = strcat( 's' , num2str(seriesNo) );

        mrecon_arrangeKT( rawDataFilePath, senseRefFilePath, coilSurveyFilePath, outputDirPath, outFilePrefix );

    end
    
else
    
    error('Must perform arrangeKT on computer with MRecon 515.');

end


%% Send k-Space Data to Remote Computer

% TODO: automate sending this data to gpubeastie03-pc, running
% recon_kt_sc_slw.m and then returning data to beastie01

fprintf('Sending to Remote Computer ...\n');

cd(outputDirPath);

filesToSend = '';

for seriesNo = seriesNos
    
    outFilePrefix = ['s' num2str(seriesNo)];
    fileStr = strcat( outFilePrefix, '_kspace_sc_slw.mat' );
    
    filesToSend = strcat( filesToSend, {' '}, fileStr );

end

if strcmp( computerName(1:end-1), 'beastie01' )
    eval(['!scp ' char(filesToSend) ' tr17@gpubeastie03-pc:' ktreconDirRemote]);
elseif strcmp( computerName(1:end-1), 'gpubeastie03-pc' )
    eval(['!scp ' char(filesToSend) ' tr17@beastie01:' ktreconDirRemote]);
else
    error('Unknown computer hostname.');
end

fprintf('\nData sent to Remote Computer.\n');


%% GPU warning

if strcmp( computerName(1:end-1), 'beastie01' )
    warning('Running solverKT on beastie01 will be slow compared to gpubeastie03-pc');
end


%% Perform Self-calibrated Sliding Window k-t Reconstruction

for seriesNo = seriesNos

        outFilePrefix = strcat( 's' , num2str(seriesNo) );
    
        recon_kt_sc_slw( outputDirPath, outFilePrefix );

end


%% Send Self-calibrated Sliding Window recon Data to Remote Computer

fprintf('Sending to Remote Computer ...\n');

cd(outputDirPath);

filesToSend = '';

for seriesNo = seriesNos
    
    outFilePrefix = ['s' num2str(seriesNo)];
    fileStr = strcat( outFilePrefix, '_sc_slw_recon.mat' );
    
    filesToSend = strcat( filesToSend, {' '}, fileStr );

end

if strcmp( computerName(1:end-1), 'beastie01' )
    eval(['!scp ' char(filesToSend) ' tr17@gpubeastie03-pc:' ktreconDirRemote]);
elseif strcmp( computerName(1:end-1), 'gpubeastie03-pc' )
    eval(['!scp ' char(filesToSend) ' tr17@beastie01:' ktreconDirRemote]);
else
    error('Unknown computer hostname.');
end

fprintf('\nData sent to Remote Computer.\n');


%% Check for Heart Masks

maskFileNames = dir( [ maskDirPath , '/s*_mask_heart.nii.gz' ] );

if isempty( maskFileNames )
    error('Heart masks not found in mask directory.');
end


%% Perform k-t SENSE Reconstruction using Self-calibrated Training Data 

% Need mRecon 557 for fcmr_ktrecon
if strcmp( computerName(1:end-1), 'beastie01' )
    warning('off');
    addpath( genpath('/home/tr17/MATLAB/MRecon-3.0.557') );
    rmpath ( genpath('/home/tr17/MATLAB/MRecon-3.0.515') );
    warning('on');
else
    error('Must perform mrecon_kt on computer with MRecon 557.');
end


% recon_exam

patchVersion       = 'PIH1';
isGeoCorrn         = true;
cusTrnDirPath      = fullfile( reconDir, 'ktrecon' );
makeHarmonicFilter = true;
isSelfCaliPreProc  = true;

recon_exam( fcmrNum, seriesNos, rawDataDirPath, patchVersion, ...
    isGeoCorrn, maskDirPath, cusTrnDirPath, makeHarmonicFilter, isSelfCaliPreProc, outputDirPath );

fprintf('\n____________ ALL STACKS RECONSTRUCTED ____________\n\n');
















