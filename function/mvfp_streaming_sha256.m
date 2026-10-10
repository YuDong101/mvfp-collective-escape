function digestHex = mvfp_streaming_sha256(filename)
if ~isfile(filename)
    error('mvfp:DynamicFoldPilot:MissingHashSource', ...
        'Cannot hash a missing file: %s', filename);
end

digest = java.security.MessageDigest.getInstance('SHA-256');
stream = java.security.DigestInputStream( ...
    java.io.FileInputStream(char(filename)), digest);
cleanup = onCleanup(@() stream.close()); %#ok<NASGU>
buffer = zeros(1, 1024 * 1024, 'int8');
while true
    count = stream.read(buffer);
    if count < 0
        break;
    end
end
bytes = typecast(digest.digest(), 'uint8');
digestHex = lower(reshape(dec2hex(bytes, 2).', 1, []));
end
