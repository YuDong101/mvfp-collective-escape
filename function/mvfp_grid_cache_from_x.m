function grid = mvfp_grid_cache_from_x(x)
x = x(:);
n = numel(x);
if n < 2
    error('MVFP grid must contain at least two points.');
end
dx = x(2) - x(1);
w = dx .* ones(n, 1);
w([1, end]) = 0.5 .* dx;

grid = struct('x', x, 'dx', dx, 'w', w, ...
    'wx', w .* x, 'wx2', w .* x.^2, 'x3', x.^3);
end

%% ============================================================
