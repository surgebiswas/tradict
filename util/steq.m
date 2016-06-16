function mask = steq( cellstr, queries )
% Performs string comparison allowing for multiple queries.
% queries should be a row or column cell array of query strings.

mask = false(size(cellstr));
for i = 1 : length(queries)
    mask = mask | strcmpi(cellstr, queries{i});
end


end

