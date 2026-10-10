function tf = parallel_loop_enabled(opts)
tf = isfield(opts, 'parallelActive') && logical(opts.parallelActive);
end

%% ============================================================
