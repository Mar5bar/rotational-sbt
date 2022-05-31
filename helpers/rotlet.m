function res = rotlet(x,y)
%% Returns the vector corresponding to a rotlet at y evaluated at x.
    r = x - y;
    norm_r_squared = sum(r.^2,1);
    res = r ./ (norm_r_squared).^(3/2);
end