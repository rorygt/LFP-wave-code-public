function croppedS = cropStructure(S, rowIndices)
% CROPSTRUCTURE crops the structure S to only select only rows specified in
% ROWINDICES. ROWINDICES can be an array of indices, a binary array or a
% single string which will be found in S.type.
% Each element in S must be a matrix with m rows and a variable
% amount of columns.

fields = fieldnames(S);

% If input is a binary array, convert it to an array of indices
if length(rowIndices) == size(S.(fields{1}), 1)
    rowIndices = find(rowIndices);
end

% If input is a string, convert it to an array of indices by finding it in
% S.type
if ischar(rowIndices)
    rowIndices = find(strcmp(rowIndices, S.type));
end

% If input is a cell array of strings, convert it to an array of indices by
% finding all in S.type
if iscell(rowIndices)
    binaryArray = zeros(size(S.type));
    for ii = 1:length(rowIndices)
        istring = rowIndices{ii};
        binaryArray(strcmp(istring, S.type)) = 1;
    end
    rowIndices = find(binaryArray);   
end

% Iterate over every field in S and keep only required entries
for i = 1:length(fields)
    ifield = fields{i};
    % Remove field if it is the wrong size
    if size(S.(ifield), 1) >= max(rowIndices)
        croppedS.(ifield) = S.(ifield)(rowIndices, :);
    else
        croppedS.(ifield) = [];
        sprintf('WARNING: Field %s deleted', ifield)
    end
end

