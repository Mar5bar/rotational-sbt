function A = genACombinedAnsatzBCApprox(in, verbose)
% We form an 6Nx6N system and enforce the boundary condition at 2N 
% points along the length of the slender body, using the surface 
% velocity at XFront and XSide to give a full rank system. The
% variables will be ordered as [f1x,f1y,...,fNz,m1x,m1y,...mNz], whilst
% the equations will be [uFront1,uSide1,uFront2,...,uSideN].
    if nargin < 2
        verbose = false;
    end    
    if verbose
        disp('Combined ansatz with leading-order rotlet integral...')    
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
    for pointInd = 1 : N
        s = segmentMidpoints(pointInd);

        %----
        % Evaluate the translational part.
        %----
        % We'll integrate the translational kernel over each segment, then
        % assign the relevant components to the matrix. We do this for
        % XFront, XBack, XSide, XOtherSide at the same time.
        integrandIntermediate = @(s,X) regularisedStokeslet(X,xi(s),regularisationParam(s)) + dipoleWeightingFactor * (e^2 - s.^2).*regularisedPotentialDipole(X,xi(s),regularisationParam(s));
        integrand = @(s) [integrandIntermediate(s,XFront(:,pointInd)); integrandIntermediate(s,XBack(:,pointInd)); integrandIntermediate(s,XSide(:,pointInd)); integrandIntermediate(s,XOtherSide(:,pointInd))];

        % Compute the integrals. ODE15S is very good at this.
        [~, sol] = ode15s(@(t,y) reshape(integrand(t),[],1), segmentEndpoints, zeros(36,1), opts);
        integrals = reshape(diff(sol)',12,3*N);

        uTransFront = integrals(1:3,:);
        uTransBack = integrals(4:6,:);
        uTransSide= integrals(7:9,:);
        uTransOtherSide= integrals(10:12,:);
        
        % These equations will isolate any angular dynamics using the difference between two velocities.
        A((pointInd-1)*6+1,1:3*N) = et(s)' * (uTransFront - uTransBack);
        A((pointInd-1)*6+2,1:3*N) = eb(s)' * (uTransFront - uTransBack);
        A((pointInd-1)*6+3,1:3*N) = et(s)' * (uTransSide - uTransOtherSide);
        
        % These equations will isolate any translational dynamics using the sum of two velocities.
        A((pointInd-1)*6+4,1:3*N) = et(s)' * (uTransFront + uTransBack);
        A((pointInd-1)*6+5,1:3*N) = en(s)' * (uTransFront + uTransBack);
        A((pointInd-1)*6+6,1:3*N) = eb(s)' * (uTransFront + uTransBack);

        %----
        % Evaluate the rotlet part.
        %----
        % We'll use the leading-order result for the rotlet integrals,
        % which include only the contribution of the local torque.
        uAngFront = zeros(3,N);
        uAngBack= zeros(3,N);
        uAngSide = zeros(3,N);
        uAngOtherSide = zeros(3,N);

        uAngFront(:,pointInd) = ep * eta(s) * I1(s) * ers.front(:,pointInd) - I2(s) * et(s);
        uAngBack(:,pointInd) = ep * eta(s) * I1(s) * ers.back(:,pointInd) - I2(s) * et(s);
        uAngSide(:,pointInd) = ep * eta(s) * I1(s) * ers.side(:,pointInd) - I2(s) * et(s);
        uAngOtherSide(:,pointInd) = ep * eta(s) * I1(s) * ers.otherSide(:,pointInd) - I2(s) * et(s);

        % These equations will isolate any angular dynamics using the difference between two velocities.
        A((pointInd-1)*6+1,3*N+1:end) = reshape(cross(uAngFront - uAngBack, repmat(et(s),1,N)),1,[]);
        A((pointInd-1)*6+2,3*N+1:end) = reshape(cross(uAngFront - uAngBack, repmat(eb(s),1,N)),1,[]);
        A((pointInd-1)*6+3,3*N+1:end) = reshape(cross(uAngSide - uAngOtherSide, repmat(et(s),1,N)),1,[]);

        % These equations will isolate any translational dynamics using the sum of two velocities.
        A((pointInd-1)*6+4,3*N+1:end) = reshape(cross(uAngFront + uAngBack, repmat(et(s),1,N)),1,[]);
        A((pointInd-1)*6+5,3*N+1:end) = reshape(cross(uAngFront + uAngBack, repmat(en(s),1,N)),1,[]);
        A((pointInd-1)*6+6,3*N+1:end) = reshape(cross(uAngFront + uAngBack, repmat(eb(s),1,N)),1,[]);

        if verbose
            textprogressbar(100 * pointInd / N)
        end
    end
    if verbose
        textprogressbar('Done!')
    end
end