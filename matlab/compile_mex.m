% COMPILE_MEX  Compile the MEX files for iboss_od_mex.
%
%   Run this script once before using iboss_od_mex.
%   The compiled binaries will be placed in the current directory.

fprintf('Compiling getIdx_mex.cpp ...\n');
mex(fullfile(fileparts(mfilename('fullpath')), 'getIdx_mex.cpp'), ...
    '-outdir', fileparts(mfilename('fullpath')));

fprintf('Compiling getIdxR_mex.cpp ...\n');
mex(fullfile(fileparts(mfilename('fullpath')), 'getIdxR_mex.cpp'), ...
    '-outdir', fileparts(mfilename('fullpath')));

fprintf('Done. MEX files ready.\n');
