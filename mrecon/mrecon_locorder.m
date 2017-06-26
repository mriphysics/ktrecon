function locOrder = mrecon_locorder( MR )
%MRECON_LOCORDER  Determine acquired loc order of MRecon Object
%
%   locOrder = MRECON_LOCORDER( MR )

%   JFPvA (joshua.vanamerom@kcl.ac.uk)


%% Get Unsorted Labels

P = MRecon( MR.Parameter.Filename.Parameter ); 


%% Determine Label Indices to Use

indTyp1Mix0 = P.Parameter.Labels.Index.typ==1 & P.Parameter.Labels.Index.mix==0;

loca = unique( P.Parameter.Labels.Index.loca( indTyp1Mix0 ) );


%% Find First Occurence of Each Loc

locFirstOccurenceIndex = nan( size( loca ) );

for iLoc = 1:numel(loca)
    locFirstOccurenceIndex( iLoc ) = find( P.Parameter.Labels.Index.loca( indTyp1Mix0 ) == loca( iLoc ), 1, 'first' );
end


%% Derive Order of Locs

[ ~, indLocOrder ] = sort( locFirstOccurenceIndex );

locOrder = loca( indLocOrder );


end  % mrecon_locorder(...)