function s_out = sortStruct(s, fname)
% SORTSTRUCT sorts every vector in the structure S by the values in field
% FNAME. Every field in S must be a vector with first dimension of length
% N, which is the same for each field.

names = fieldnames(s);
[~, order] = sort(s.(fname));

for iname = 1:length(names)
    name = names{iname};
    dims = repmat({':'}, 1, ndims(s.(name)) - 1);
    s_out.(name) = s.(name)(order, dims{:});
end