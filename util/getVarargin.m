function [ argument ] = getVarargin( vararginArr, query )
% For name-value argument pairs

% Useful when varargin is passed down through multiple function calls. 
% Unwrap the varargin like an onion.
if ~isempty(vararginArr)
    while ~ischar(vararginArr{1})
        vararginArr = vararginArr{1};
        if isempty(vararginArr)
            argument = [];
            return
        end
    end
end



ind = find(strcmpi(vararginArr, query)) + 1;
if isempty(ind)
    argument = [];
else
    argument = vararginArr{ind};
end

end

