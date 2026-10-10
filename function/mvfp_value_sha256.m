function digestHex = mvfp_value_sha256(value)
digest = java.security.MessageDigest.getInstance('SHA-256');
[~,~,machineEndian] = computer;
update_value(value);
raw = typecast(digest.digest(),'uint8');
digestHex = lower(reshape(dec2hex(raw,2).',1,[]));

    function update_value(item)
        update_text(class(item));
        update_size(size(item));
        if isnumeric(item)
            update_text('numeric');
            update_text(ternary(isreal(item),'real','complex'));
            update_numeric(real(item));
            if ~isreal(item)
                update_numeric(imag(item));
            end
        elseif islogical(item)
            update_text('logical');
            update_bytes(uint8(item(:)));
        elseif ischar(item)
            update_text('char');
            update_bytes(uint8(unicode2native(item(:).','UTF-8')));
        elseif isstring(item)
            update_text('string');
            for itemIndex = 1:numel(item)
                if ismissing(item(itemIndex))
                    update_text('missing');
                else
                    update_text(char(item(itemIndex)));
                end
            end
        elseif iscell(item)
            update_text('cell');
            for itemIndex = 1:numel(item)
                update_value(item{itemIndex});
            end
        elseif isstruct(item)
            update_text('struct');
            names = sort(fieldnames(item));
            update_size([numel(names) 1]);
            for itemIndex = 1:numel(item)
                for fieldIndex = 1:numel(names)
                    update_text(names{fieldIndex});
                    update_value(item(itemIndex).(names{fieldIndex}));
                end
            end
        else
            update_text('json-fallback');
            update_text(jsonencode(item));
        end
    end

    function update_numeric(item)
        if machineEndian == 'B' && ~isa(item,'int8') && ...
                ~isa(item,'uint8')
            item = swapbytes(item);
        end
        update_bytes(typecast(item(:),'uint8'));
    end

    function update_size(dimensions)
        update_text(sprintf('%d,',dimensions));
    end

    function update_text(textValue)
        bytes = uint8(unicode2native(char(textValue),'UTF-8'));
        lengthBytes = typecast(uint64(numel(bytes)),'uint8');
        if machineEndian == 'B'
            lengthBytes = fliplr(lengthBytes);
        end
        update_bytes(lengthBytes);
        update_bytes(bytes);
    end

    function update_bytes(bytes)
        if isempty(bytes)
            return;
        end
        digest.update(typecast(uint8(bytes(:)),'int8'));
    end
end

function value = ternary(condition,trueValue,falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end
