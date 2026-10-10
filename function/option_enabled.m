function tf = option_enabled(s, fieldName, defaultValue)
if isfield(s, fieldName)
    tf = logical(s.(fieldName));
else
    tf = defaultValue;
end
end

%% ============================================================
