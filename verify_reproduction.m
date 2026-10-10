% Small independent calculations covering the shared numerical library.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'tests'));
verify_compact_numerics(root,fullfile(root,'verification','quick'));
