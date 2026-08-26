% Absolute path to a file in the repo data directory, independent of the working
% directory. Figures run from the repo root, the slope computer runs from codes/,
% so a relative 'data/...' is only correct for one of them.
%
% You provide:
% *name - file name within data/
%
% Returns:
% *p - absolute path
%
% Nathan Laxague 2026
%
function p = repo_data_path(name)

here = fileparts(mfilename('fullpath'));
p = fullfile(fileparts(here),'data',name);
