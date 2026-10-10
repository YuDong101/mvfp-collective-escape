function base = merge_struct_fields(base, override)
if ~isstruct(override)
    error('Override must be a struct.');
end
fields = fieldnames(override);
for k = 1:numel(fields)
    name = fields{k};
    if isfield(base, name) && isstruct(base.(name)) && isstruct(override.(name))
        base.(name) = merge_struct_fields(base.(name), override.(name));
    else
        base.(name) = override.(name);
    end
end
end

%% ============================================================
