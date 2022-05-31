function R = genRCombinedAnsatzRTTBCApprox(in)
% Form the RHS of the linear system for the combined ansatz, which is 
% the velocity on the surface at XFront and XSide. 
    uFront = in.uFront;
    Omega = in.Omega;

    R = [uFront(:); Omega(:)];
end