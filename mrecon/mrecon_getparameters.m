function PARAM = mrecon_getparameters( MR )
%MRECON_GETPARAMETERS  Get MRecon object parameters as Matlab struct
%
%   PARAM = MRECON_GETPARAMETERS( MR )

%   JFPvA (joshua.vanamerom@kcl.ac.uk)


%% Parameter Groups

paramGroup = fieldnames( MR.Parameter );

paramGroup2Read = { 'Scan', 'Encoding', 'Labels', 'Filename' };
    

%% Extract Parameters from MRecon Object

for iParamGroup = 1:numel(paramGroup2Read)
    
    if any( ismember( paramGroup, paramGroup2Read{iParamGroup} ) )
        fieldname = fieldnames( MR.Parameter.(paramGroup2Read{iParamGroup}) );
        for iField = 1:numel(fieldname)
            PARAM.(paramGroup2Read{iParamGroup}).(fieldname{iField})  = MR.Parameter.(paramGroup2Read{iParamGroup}).(fieldname{iField});
        end
    end
    
end


end  % mrecon_getparameters(...)
