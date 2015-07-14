    function sReduced = selectStructEntries(s, entries)
    % Selects entries from structural array S where each field in S is a
    % (possibly multidimensional) array with the same size in the first
    % dimension. ENTRIES must be a binary array of elements to select.
    names = fieldnames(s);
    for ii = 1:length(names)
        iname = names{ii};
        isize = size(s.(iname));
        if isize(1) == length(entries)
            temparray = s.(iname);
            temparray = temparray(repmat(entries, [1, isize(2:end)]));
            sReduced.(iname) = reshape(temparray, ...
                [sum(entries), isize(2:end)]);
        else
            sReduced.(iname) = s.(iname);
        end
    end
    
    
    end