function A = genARotletAnsatzRTT(in, verbose);
% We form an 3Nx3N system and enforce the boundary condition at N 
% points along the length of the slender body. The RTT directly relates
% the local torque to the local angular velocity, and we need only to
% evaluate the integral I1(s) and form diagonal blocks.
    if nargin < 2
        verbose = false;
    end
    if verbose
        disp('Rotlet ansatz using RTT...')
    end
    A = kron(diag(in.I1(in.segmentMidpoints)),eye(3));
end