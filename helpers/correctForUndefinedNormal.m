function r = correctForUndefinedNormal(n,t)
%% If the tangent is locally constant, then n will be NaN. We make a choice to
% fix n in terms of the local tangent t.
    if any(isnan(n))
        % Try the cross product with [1;0;0].
        r = cross(t,[1;0;0]);
        r = r / norm(r);
        % If this is still NaN, then t = [1;0;0] and we instead cross it with
        % [0;0;-1] to give a suitable n.
        if any(isnan(r))
            r = cross(t,[0;0;-1]);
            r = r / norm(r);
        end
    else
        r = n;
    end
end