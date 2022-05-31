function R = genRRotletAnsatz(in)
% Form the RHS of the linear system for the rotlet ansatz, which is 
% various dot products of uFront and uSide with et and eb.
    
    % Unpack the struct.
    N = in.N;
    uFront = in.uFront;
    uSide = in.uSide;
    et = in.et;
    eb = in.eb;
    segmentMidpoints = in.segmentMidpoints;

    R = zeros(3*N,1);
    for pointInd = 1 : N
        s = segmentMidpoints(pointInd);
        R((pointInd-1)*3 + 1) = dot(uFront(:,pointInd), et(s));
        R((pointInd-1)*3 + 2) = dot(uFront(:,pointInd), eb(s));
        R((pointInd-1)*3 + 3) = dot(uSide(:,pointInd), et(s));
    end
end