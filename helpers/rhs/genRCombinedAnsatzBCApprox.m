function R = genRCombinedAnsatzBCApprox(in)
% Form the RHS of the linear system, which is 
% various dot products of the velocities with et and eb.
    N = in.N;
    uFront = in.uFront;
    uSide = in.uSide;
    uBack = in.uBack;
    uOtherSide = in.uOtherSide;
    et = in.et;
    eb = in.eb;
    en = in.en;
    segmentMidpoints = in.segmentMidpoints;

    R = zeros(6*N,1);
    for pointInd = 1 : N
        s = segmentMidpoints(pointInd);

        R((pointInd-1)*6 + 1) = dot(uFront(:,pointInd) - uBack(:,pointInd), et(s));
        R((pointInd-1)*6 + 2) = dot(uFront(:,pointInd) - uBack(:,pointInd), eb(s));
        R((pointInd-1)*6 + 3) = dot(uSide(:,pointInd) - uOtherSide(:,pointInd), et(s));

        R((pointInd-1)*6 + 4) = dot(uFront(:,pointInd) + uBack(:,pointInd), et(s));
        R((pointInd-1)*6 + 5) = dot(uFront(:,pointInd) + uBack(:,pointInd), en(s));
        R((pointInd-1)*6 + 6) = dot(uFront(:,pointInd) + uBack(:,pointInd), eb(s));
    end
end