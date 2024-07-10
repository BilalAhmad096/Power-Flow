function [out] = randomvalues(lower,upper)
out = (upper-lower).*rand(500,1) + lower;
r_range = [min(out) max(out)];
end