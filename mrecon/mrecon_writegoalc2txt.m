function mrecon_writegoalc2txt( MR, filePath, varargin )
%MRECON_WRITEGOALC2TXT  write MRecon object GoalC information to text file
%
%   MRECON_WRITEGOALC2TXT( MR, filePath )
%
%   MRECON_WRITENIFTI2D( ... , 'param', val) 
% 
%   Input:
%       MR                  - ReconFrame MRecon object
%       filePath            - file path to output text file
%
%   Parameter-Value Options:
%       objname             - object names to write to file, 
%                             as cell array of character vectors

% JFPvA (joshua.vanamerom@kcl.ac.uk)


%% Default Values for Unspecified Input Arguments

default.objname        = {  'SQ`prebase', ...
                            'GR`s_prebase1', 'GR`r_prebase1', ...
                            'RF`prebase1', ...
                            'SQ`base', ...
                            'GR`mc[0]', 'GR`m[0]', 'GR`md', 'GR`m[3]', ...
                            'AQ`base', ...
                            'GR`py', 'GR`pyr', ...
                            'GR`s_ex', 'GR`r[0]', 'GR`d', ...
                            'RF`ex' };


%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(    p, 'MR', ...
    @(x) validateattributes( x, {'MRecon'}, {'scalar'}, mfilename) ); 
addRequired(    p, 'filePath', ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) ); 
add_param_fn(   p, 'objname', default.objname, ...
    @(x) validateattributes( x, {'cell'}, {}, mfilename) );

parse( p, MR, filePath, varargin{:} );

objname         = p.Results.objname;


%% Write Objects to File

fid = fopen( filePath, 'w' );

for iObj = 1:numel( objname )

    if ( MR.Parameter.IsObject(objname{iObj}) )
    
        T = evalc( 'MR.Parameter.DisplayObject( objname{iObj} )' );
        fprintf( fid, '%s\n\n', T );
        
    else
        
        warning( 'object ''%s'' not found', objname{iObj} )
    
    end

end

fclose( fid );


end  % mrecon_writegoalc2txt(...)
