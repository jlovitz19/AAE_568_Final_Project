thisScriptPath = matlab.desktop.editor.getActiveFilename;
thisFolder = fileparts(thisScriptPath);
targetPath = fullfile(thisFolder, 'attitude_fxns');
addpath(targetPath);

%{
to use this script, paste the following lines in your primary (runner) script 
(this only works if your script is in the same folder as add_to_path.m, you
may need to adjust)

script_path = matlab.desktop.editor.getActiveFilename;
run(fullfile(fileparts(mfilename('fullpath')), 'add_to_path.m'));
%}