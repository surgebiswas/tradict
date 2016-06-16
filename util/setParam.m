function [ val ] = setParam( vararg, param, default )

p = getVarargin(vararg, param);
if isempty(p)
    val = default;
else
    val = p;
end

end

