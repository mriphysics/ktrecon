function fname = filename(file)
%FILENAME File name including extension of file.
%   FNAME = FILENAME(FILE) returns the file name including extenstion
%   extension for the specified FILE. The FILE input is the name of a file,
%   and can include a path and file name extension. The function interprets all
%   characters following the right-most path delimiter as a file name plus extension.
%
%   See also FILEPARTS, FULLFILE, PATHSEP, FILESEP.

% JFPvA (joshua.vanamerom@kcl.ac.uk)

[ ~, fname, fext ] = fileparts( file );
fname = strcat( fname, fext );

end
