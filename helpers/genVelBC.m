function [uFront, uSide, uBack, uOtherSide, Omega] = genVelBC(linVel, angVel, segmentMidpoints, eta, ers, params, localBasis)
    V = linVel(segmentMidpoints);
    N = params.N;
    ep = params.epsilon;
    Omega = zeros(3,N);
    uFront = zeros(3,N);
    uSide = zeros(3,N);
    uBack = zeros(3,N);
    uOtherSide = zeros(3,N);
    for pointInd = 1 : N
        s = segmentMidpoints(pointInd);
        Omega(:,pointInd) = localBasis(s) * angVel(s);
        uFront(:,pointInd) = V(:,pointInd) + ep * eta(s) * cross(Omega(:,pointInd), ers.front(:, pointInd));
        uSide(:,pointInd) =  V(:,pointInd) + ep * eta(s) * cross(Omega(:,pointInd), ers.side(:, pointInd));
        uBack(:,pointInd) =  V(:,pointInd) + ep * eta(s) * cross(Omega(:,pointInd), ers.back(:, pointInd));
        uOtherSide(:,pointInd) =  V(:,pointInd) + ep * eta(s) * cross(Omega(:,pointInd), ers.otherSide(:, pointInd));
    end
end