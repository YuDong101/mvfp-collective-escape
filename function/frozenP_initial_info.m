function info = frozenP_initial_info(P0, par)
info = struct('valid', 0, 'xLstar', nan, 'xRstar', nan);
F = initial_field_guess(P0, par);
[ok, xL, ~, xR] = extrema_from_F(F, par);
if ok
    info.valid = 1;
    info.xLstar = xL;
    info.xRstar = xR;
end
end

%% ============================================================
