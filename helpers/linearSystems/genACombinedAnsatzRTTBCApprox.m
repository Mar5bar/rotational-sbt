function A = genACombinedAnsatzRTTBCApprox(in, verbose)
% We form an 6Nx6N system and enforce the boundary condition at N 
% points along the length of the slender body, using the leading-order
% boundary condition at XFront, and the RTT result at these N points to
% give a full rank system. The variables will be ordered as 
% [f1x,f1y,...,fNz,m1x,m1y,...mNz], whilst the equations will be 
% [uFront1,...uFrontN,RTT1,...,RTTN].
    if nargin < 2
        verbose = false;
    end    
    if verbose
        disp('Combined ansatz with RTT and leading-order rotlet integral...')    
        textprogressbar('Building matrix: ')
    end

    % Unpack the struct.
    N = in.N;
    ep = in.ep;
    segmentMidpoints = in.segmentMidpoints;
    segmentEndpoints = in.segmentEndpoints;
    XFront = in.XFront;
    XSide = in.XSide;
    XBack = in.XBack;
    XOtherSide = in.XOtherSide;
    xi = in.xi;
    eta = in.eta;
    et = in.et;
    en = in.en;
    eb = in.eb;
    opts = in.opts;
    dipoleWeightingFactor = in.dipoleWeightingFactor;
    regularisationParam = in.regularisationParam;
    e = in.e;
    I1 = in.I1;
    I2 = in.I2;
    ers = in.ers;

    A = zeros(6*N);
    if verbose
        textprogressbar('Building matrix: ')
    end
    for pointInd = 1 : N
    
        %----
        % Evaluate the translational part.
        %----
        % We'll integrate the translational kernel over each segment, then
        % assign the relevant components to the matrix.
        integrandIntermediate = @(s,X) regularisedStokeslet(X,xi(s),regularisationParam(s)) + dipoleWeightingFactor * (e^2 - s.^2).*regularisedPotentialDipole(X,xi(s),regularisationParam(s));
        integrand = @(s) integrandIntermediate(s,XFront(:,pointInd));

        % We'll perform the integrals over different segments in parallel,
        % as they are quite expensive.
        % integrals = zeros(3,3,N);
        % parfor segInd = 1 : N
        %     integrals(:,:,segInd) = integral(integrand, segmentEndpoints(segInd), segmentEndpoints(segInd+1), 'ArrayValued', true, 'AbsTol', tol);
        % end
        [~, sol] = ode15s(@(t,y) reshape(integrand(t),[],1), segmentEndpoints, zeros(9,1), opts);

        % Assign the computed integrals to the linear system.
        % A((pointInd-1)*3+1 : (pointInd-1)*3+3, 1:3*N) = reshape(integrals,3,3*N);
        A((pointInd-1)*3+1 : (pointInd-1)*3+3, 1:3*N) = reshape(diff(sol)',3,3*N);

        %----
        % Evaluate the rotlet part.
        %----
        % We'll use the leading-order result for the rotlet integrals,
        % which include only the contribution of the local torque.
        s = segmentMidpoints(pointInd);
        A((pointInd-1)*3+1 : (pointInd-1)*3+3, 3*N + (pointInd-1)*3+1 : 3*N + (pointInd-1)*3+3) = crossProductMatrix(ep * eta(s) * I1(s) * ers.front(:,pointInd) - I2(s) * et(s));

        if verbose
            textprogressbar(100 * pointInd / N)
        end
    end

    %----
    % Evaluate the RTT part all at once.
    %----
    % We'll use the full RTT relation at the N points.
    A(3*N+1:6*N,3*N+1:6*N) = kron(diag(I1(segmentMidpoints)),eye(3));

    if verbose
        textprogressbar('Done!')
    end
end