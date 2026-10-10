function grid = mvfp_grid_cache(cfg)
x = linspace(-cfg.xMax, cfg.xMax, cfg.xGridN).';
grid = mvfp_grid_cache_from_x(x);
if isfield(cfg,'fluxScheme') && ~isempty(cfg.fluxScheme)
    fluxScheme = char(cfg.fluxScheme);
    switch fluxScheme
        case {'upwind','exponential-fit'}
            grid.fluxScheme = fluxScheme;
        otherwise
            error('MVFP:Grid:UnknownFluxScheme', ...
                'Unknown fluxScheme: %s',fluxScheme);
    end
end
if isfield(cfg,'massDiscretization') && ~isempty(cfg.massDiscretization)
    switch char(cfg.massDiscretization)
        case 'nodal-control-volume'
            % Endpoints represent half cells at the physical no-flux boundary.
            grid.volumeWeights = grid.w;
        case 'legacy-uniform-divergence'
            % Preserve the historical arithmetic, including endpoint scaling.
        otherwise
            error('MVFP:Grid:UnknownMassDiscretization', ...
                'Unknown massDiscretization: %s',char(cfg.massDiscretization));
    end
end
end

%% ============================================================
