function output = evaluateSBT(params, methods)
    tic
    methodsToDo = struct();
    if nargin < 2 
        methodsToDo.rotletAnsatz = true;
        methodsToDo.rotletAnsatzRTT = true;
        methodsToDo.rotletAnsatzRTTI1Approx = true;
        methodsToDo.combinedAnsatz = true;
        methodsToDo.combinedAnsatzBCApprox = true;
        methodsToDo.combinedAnsatzRTT = true;
        methodsToDo.combinedAnsatzRTTBCApprox = true;
    else
        methodsToDo.rotletAnsatz = false;
        methodsToDo.rotletAnsatzRTT = false;
        methodsToDo.rotletAnsatzRTTI1Approx = false;
        methodsToDo.combinedAnsatz = false;
        methodsToDo.combinedAnsatzBCApprox = false;
        methodsToDo.combinedAnsatzRTT = false;
        methodsToDo.combinedAnsatzRTTBCApprox = false;
        for method = methods
            methodsToDo.(method{1}) = true;
        end
    end

    %% Load in parameters and compute the slender-body setup.
    % Slenderness parameter.
    ep = params.epsilon;
    % Discretisation level.
    N = params.N;
        
    % Absolute tolerance for integration.
    if isfield(params,'tol')
        tol = params.tol;
    else
        tol = 1e-12;
    end

    % Select a parameterised body centreline function (x,y,z)(t) and derivative
    % (dx,dy,dz)(t) from centrelines.m.
    [x, y, z, dx, dy, dz, tRange, centrelineLab] = centrelines(params.xiSelector);
    disp(['xi(t) = ', centrelineLab, ', t in [',num2str(tRange),']'])

    % Select a shape function eta(s) from etas.m.
    [eta, etaLab] = etas(params.etaSelector);
    disp(['eta(s) = ', etaLab])
    % Redefine eta so that it has max 1.
    m = max(eta(linspace(-1,1,1e5)));
    eta = @(s) eta(s) / m;

    % Select a linear velocity linVel(s) from linVels.m.
    [linVel, linVelLab] = linVels(params.linVelSelector);
    disp(['V(s) = ', linVelLab])

    % Select an angular velocity angVel(s) from angVels.m.
    % angVel has components in the et, en, and eb directions.
    [angVel, angVelLab] = angVels(params.angVelSelector);
    disp(['Omega(s) = ', angVelLab])

    % We will want to have all of our quantities as functions of arclength,
    % not the parameter t. Hence, compute a map from arclength to t,
    % recalling that arclength ranges from -1 to 1.
    sampleParams = linspace(tRange(1), tRange(2), 1e4);
    arclengthsOriginal = [0,cumsum(sqrt(diff(x(sampleParams)).^2 + diff(y(sampleParams)).^2 + diff(z(sampleParams)).^2))];
    arclengthToParam = griddedInterpolant(2*arclengthsOriginal/arclengthsOriginal(end) - 1, sampleParams, 'linear', 'linear');

    % Now redefine x,y,z and dx,dy,dz in terms of arclength. Note that dx
    % etc are still derivatives wrt the parameter t, not s, but this will
    % not matter in what follows.
    x = @(s) 2*x(arclengthToParam(s)) / arclengthsOriginal(end);
    y = @(s) 2*y(arclengthToParam(s)) / arclengthsOriginal(end);
    z = @(s) 2*z(arclengthToParam(s)) / arclengthsOriginal(end);
    dx = @(s) 2*dx(arclengthToParam(s)) / arclengthsOriginal(end);
    dy = @(s) 2*dy(arclengthToParam(s)) / arclengthsOriginal(end);
    dz = @(s) 2*dz(arclengthToParam(s)) / arclengthsOriginal(end);

    % Define the centreline, xi(s).
    xi = @(s) [x(s); y(s); z(s)];

    % Define the forms of the Frenet-Serret basis vectors et, en, eb.
    et = @(s) [dx(s);dy(s);dz(s)]; et = @(s) et(s) / norm(et(s));
    ds = 1e-4; en = @(s) et(s+ds) - et(s-ds); en = @(s) en(s) / norm(en(s));
    % en will be NaN if et was locally constant, as en is formally
    % undefined in this case. We flag this occurrence to the user and set
    % en to be et x constant vector, the latter being either [1;0;0] or [0;0;-1].
    en = @(s) correctForUndefinedNormal(en(s),et(s));
    eb = @(s) cross(et(s), en(s));
    
    % Define er in terms of en and eb.
    er = @(s, phi) cos(phi) .* en(s) + sin(phi) .* eb(s);
    
    %% Compute derived quantities for use in integration.
    % Compute the effective eccentricity, used in the limits of integration
    % and the composite ansatz.
    e = sqrt(1 - ep^2);

    % Precompute the weighting factor applied to the potential dipole in the translational ansatz.
    dipoleWeightingFactor = - ep^2 ./ (2*e^2);
    
    % Define the arclength-dependent regularisation parameter used in the translational ansatz.
    regularisationParam = @(s) ep^2 * (1-s.^2 - eta(s).^2);

    % Decompose the integration domain into N elements, with endpoints segmentEndpoints.
    segmentEndpoints = linspace(-e,e,N+1);
    segmentMidpoints = movmean(segmentEndpoints,2,'Endpoints','discard'); % These points are at the middle of each segment.
    
    % Define the integrals used for constructing the leading-order rotlet
    % ansatz and the RTT.
    I1 = @(s) ((e - s)./((s-e).^2 + ep^2*eta(s).^2).^(0.5) + (s + e)./((s+e).^2 + ep^2*eta(s).^2).^(0.5)) ./ (ep^2 * eta(s).^2);
    I2 = @(s) (-1./((s-e).^2 + ep^2*eta(s).^2).^(0.5) + 1./((s+e).^2 + ep^2*eta(s).^2).^(0.5));
    
    % Define the integral used in the approximate form of the RTT.
    I1Approx = @(s) 2 ./ (ep^2 * eta(s).^2);
    
    % We'll evaluate the boundary condition at 2N points for the combined
    % ansatz, two per segment, and at N points for the pure rotlet ansatz.
    % We'll make use of 
    % XFront = xi + epsilon * eta * er(s,0), 
    % XSide = xi + epsilon * eta * er(s,pi/2), and
    % XBack = xi - epsilon * eta * er(s,0).
    erAtSegmentMidpoints = zeros(3,N);
    for pointInd = 1 : N
        erAtSegmentMidpoints(:,pointInd) = er(segmentMidpoints(pointInd), 0);
    end
    erFrontAtSegmentMidpoints = cell2mat(arrayfun(@(s) er(s,0), segmentMidpoints, 'UniformOutput', false));
    erSideAtSegmentMidpoints = cell2mat(arrayfun(@(s) er(s,pi/2), segmentMidpoints, 'UniformOutput', false));
    erBackAtSegmentMidpoints = cell2mat(arrayfun(@(s) er(s,pi), segmentMidpoints, 'UniformOutput', false));
    erOtherSideAtSegmentMidpoints = cell2mat(arrayfun(@(s) er(s,3*pi/2), segmentMidpoints, 'UniformOutput', false));
    
    XFront = xi(segmentMidpoints) + ep * eta(segmentMidpoints) .* erFrontAtSegmentMidpoints;
    XSide= xi(segmentMidpoints) + ep * eta(segmentMidpoints) .* erSideAtSegmentMidpoints;
    XBack = xi(segmentMidpoints) + ep * eta(segmentMidpoints) .* erBackAtSegmentMidpoints;
    XOtherSide= xi(segmentMidpoints) + ep * eta(segmentMidpoints) .* erOtherSideAtSegmentMidpoints;

    % Now evaluate the velocity at these points.
    V = linVel(segmentMidpoints);
    Omega = zeros(3,N);
    uFront = zeros(3,N);
    uSide = zeros(3,N);
    uBack = zeros(3,N);
    uOtherSide = zeros(3,N);
    for pointInd = 1 : N
        s = segmentMidpoints(pointInd);
        localBasis = [et(s), en(s), eb(s)];
        Omega(:,pointInd) = localBasis * angVel(s);
        uFront(:,pointInd) = V(:,pointInd) + ep * eta(s) * cross(Omega(:,pointInd), erFrontAtSegmentMidpoints(:, pointInd));
        uSide(:,pointInd) =  V(:,pointInd) + ep * eta(s) * cross(Omega(:,pointInd), erSideAtSegmentMidpoints(:, pointInd));
        uBack(:,pointInd) =  V(:,pointInd) + ep * eta(s) * cross(Omega(:,pointInd), erBackAtSegmentMidpoints(:, pointInd));
        uOtherSide(:,pointInd) =  V(:,pointInd) + ep * eta(s) * cross(Omega(:,pointInd), erOtherSideAtSegmentMidpoints(:, pointInd));
    end

    % Generate the final evaluation points and evaluate the prescribed velocity 
    % at these points. We won't sample directly at the tips of the slender 
    % body.
    arclengthsEvaluation = linspace(-1, 1, params.numArclengthEvaluationPoints + 2);
    arclengthsEvaluation = arclengthsEvaluation(2 : end-1);
    anglesEvaluation = linspace(0, 2*pi, params.numCircumferentialEvaluationPoints + 1);
    anglesEvaluation = anglesEvaluation(1 : end-1);
    numEvalPoints = length(arclengthsEvaluation)*length(anglesEvaluation);
    evaluationPoints = zeros(3,numEvalPoints);
    velPrescribed = zeros(3,numEvalPoints);
    for i = 1 : length(arclengthsEvaluation)
        s = arclengthsEvaluation(i);
        ers = cell2mat(arrayfun(@(phi) er(s, phi), anglesEvaluation, 'UniformOutput', false));
        evaluationPoints(:,(i-1)*length(anglesEvaluation)+1 : i*length(anglesEvaluation)) = xi(s) + ep * eta(s) * ers;
        localBasis = [et(s), en(s), eb(s)];
        velPrescribed(:,(i-1)*length(anglesEvaluation)+1 : i*length(anglesEvaluation)) = linVel(s) + ep * eta(s) * cross(repmat(localBasis*angVel(s),1,length(anglesEvaluation)), ers);
    end

    %% Package output.
    output = struct();
    output.epsilon = ep;
    output.N = N;
    output.e = e;
    output.centrelineVals = xi(segmentMidpoints);
    output.evaluationPoints = evaluationPoints;
    output.boundaryBCAppliedPoints = [XFront, XSide];
    output.discreteSegmentArclengthMidpoints = segmentMidpoints;
    output.velPrescribed = velPrescribed;
    output.evaluationArclengths = arclengthsEvaluation;
    output.evaluationAngles = anglesEvaluation;

    % We'll form interpolants of all computed forces and torques.
    interps = struct();

    %% Form and solve the linear systems.
    opts = odeset('AbsTol',tol,'RelTol',tol);
    if methodsToDo.rotletAnsatz
        name = 'rotletAnsatz';
        %% Full rotlet ansatz.
        % We form an 3Nx3N system and enforce the boundary condition at N 
        % points along the length of the slender body, using the surface 
        % velocity at XFront and XSide to give a full rank system.
        disp('Rotlet ansatz...')
        output.(name) = struct();
        ARotletAnsatz = zeros(3*N);
        textprogressbar('Building matrix: ')
        for pointInd = 1 : N
            s = segmentMidpoints(pointInd);

            % We'll want to integrate m x rotlet over each segment of the
            % centreline, where m is constant on each segment. To do this,
            % we'll integrate the rotlet kernel over each segment, then assign
            % the relevant components to the matrix to perform the cross
            % product. We do this for XFront and XSide at the same time.
            integrand = @(s) [rotlet(XFront(:,pointInd), xi(s)); rotlet(XSide(:,pointInd), xi(s))];

            % Compute the integrals. ODE15S is very good at this.
            [~, sol] = ode15s(@(t,y) integrand(t), segmentEndpoints, zeros(6,1), opts);
            integrals = diff(sol)';
            

            % Now assign the components of the computed integrals to the
            % matrix serially. For each segment with midpoint s we are solving
            % uFront(s) dot et(s) = sum over segments j( 
            %                   torque(j) dot (integral over segment j of 
            %                   R(XFront(s),xi(s')) ds')
            %                   cross et(s).
            % uFront(s) dot eb(s) = sum over segments j( 
            %                   torque(j) dot (integral over segment j of 
            %                   R(XFront(s),xi(s')) ds')
            %                   cross eb(s).
            % uSide(s) dot et(s) = sum over segments j( 
            %                   torque(j) dot (integral over segment j of 
            %                   R(XSide(s),xi(s')) ds')
            %                   cross et(s).
            
            % First equation, computing the cross products of the XFront integrals
            % with et(segmendMidpoints(pointInd)).
            ARotletAnsatz((pointInd-1)*3+1,:) = reshape(cross(integrals(1:3,:), repmat(et(s), 1, N)),[],1);

            % Second equation, computing the cross products of the XFront integrals
            % with eb(segmendMidpoints(pointInd)).
            ARotletAnsatz((pointInd-1)*3+2,:) = reshape(cross(integrals(1:3,:), repmat(eb(s), 1, N)),[],1);

            % Third equation, computing the cross products of the XSide integrals
            % with et(segmendMidpoints(pointInd)).
            ARotletAnsatz((pointInd-1)*3+3,:) = reshape(cross(integrals(4:6,:), repmat(et(s), 1, N)),[],1);
            textprogressbar(100 * pointInd / N)
        end
        textprogressbar('Done!')

        % Form the RHS of the linear system for the rotlet ansatz, which is 
        % various dot products of uFront and uSide with et and eb.
        RRotletAnsatz = zeros(3*N,1);
        for pointInd = 1 : N
            RRotletAnsatz((pointInd-1)*3 + 1) = dot(uFront(:,pointInd), et(s));
            RRotletAnsatz((pointInd-1)*3 + 2) = dot(uFront(:,pointInd), eb(s));
            RRotletAnsatz((pointInd-1)*3 + 3) = dot(uSide(:,pointInd), et(s));
        end

        % Solve the linear system for the torque density.
        sol = ARotletAnsatz \ RRotletAnsatz;
        output.(name).torque = reshape(sol,3,[]);
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end

    if methodsToDo.rotletAnsatzRTT
        name = 'rotletAnsatzRTT';
        output.(name) = struct();
        %% Rotlet ansatz using RTT.
        % We form an 3Nx3N system and enforce the boundary condition at N 
        % points along the length of the slender body. The RTT directly relates
        % the local torque to the local angular velocity, and we need only to
        % evaluate the integral I1(s) and form diagonal blocks.
        disp('Rotlet ansatz using RTT...')
        ARotletAnsatzRTT = kron(diag(I1(segmentMidpoints)),eye(3));

        % Form the RHS of the linear system for the rotlet ansatz with RTT, 
        % which is simply Omega(s).
        RRotletAnsatzRTT = Omega(:);

        % Solve the linear system for the torque density.
        sol = ARotletAnsatzRTT \ RRotletAnsatzRTT;
        output.(name).torque = reshape(sol,3,[]);
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end

    if methodsToDo.rotletAnsatzRTTI1Approx
        name = 'rotletAnsatzRTTI1Approx';
        %% Rotlet ansatz using approximate RTT.
        % We form an 3Nx3N system and enforce the boundary condition at N 
        % points along the length of the slender body. The approximate RTT 
        % directly relates the local torque to the local angular velocity, and 
        % we need only to evaluate the approximation to the integral I1(s) and 
        % form diagonal blocks.
        disp('Rotlet ansatz using RTT with I1 approximated...')
        output.(name) = struct();
        ARotletAnsatzRTTI1Approx = kron(diag(I1Approx(segmentMidpoints)),eye(3));

        % Form the RHS of the linear system for the rotlet ansatz with approximate RTT, 
        % which is simply Omega(s).
        RRotletAnsatzRTTApprox = Omega(:);

        % Solve the linear system for the torque density.
        sol = ARotletAnsatzRTTI1Approx \ RRotletAnsatzRTTApprox;
        output.(name).torque = reshape(sol,3,[]);
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end

    if methodsToDo.combinedAnsatz
        name = 'combinedAnsatz';
        %% Full combined ansatz.
        % We form an 6Nx6N system and enforce the boundary condition at 2N 
        % points along the length of the slender body, using the surface 
        % velocity at XFront and XSide to give a full rank system. The
        % variables will be ordered as [f1x,f1y,...,fNz,m1x,m1y,...mNz], whilst
        % the equations will be [uFront1,uSide1,uFront2,...,uSideN].
        disp('Combined ansatz...')
        output.(name) = struct();
        ACombinedAnsatz = zeros(6*N);
        textprogressbar('Building matrix: ')
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
            ACombinedAnsatz((pointInd-1)*6+1,1:3*N) = et(s)' * (uTransFront - uTransBack);
            ACombinedAnsatz((pointInd-1)*6+2,1:3*N) = eb(s)' * (uTransFront - uTransBack);
            ACombinedAnsatz((pointInd-1)*6+3,1:3*N) = et(s)' * (uTransSide - uTransOtherSide);
            
            % These equations will isolate any translational dynamics using the sum of two velocities.
            ACombinedAnsatz((pointInd-1)*6+4,1:3*N) = et(s)' * (uTransFront + uTransBack);
            ACombinedAnsatz((pointInd-1)*6+5,1:3*N) = en(s)' * (uTransFront + uTransBack);
            ACombinedAnsatz((pointInd-1)*6+6,1:3*N) = eb(s)' * (uTransFront + uTransBack);

            %----
            % Evaluate the rotlet part.
            %----
            % We'll want to integrate m x rotlet over each segment of the
            % centreline, where m is constant on each segment. To do this,
            % we'll integrate the rotlet kernel over each segment, then assign
            % the relevant components to the matrix to perform the cross
            % product. We do this for XFront and XSide at the same time.
            integrand = @(s) [rotlet(XFront(:,pointInd), xi(s)); rotlet(XBack(:,pointInd), xi(s)); rotlet(XSide(:,pointInd), xi(s)); rotlet(XOtherSide(:,pointInd), xi(s))];

            % Compute the integrals. ODE15S is very good at this.
            [~, sol] = ode15s(@(t,y) integrand(t), segmentEndpoints, zeros(12,1), opts);
            integrals = diff(sol)';

            uAngFront = integrals(1:3,:);
            uAngBack= integrals(4:6,:);
            uAngSide = integrals(7:9,:);
            uAngOtherSide = integrals(10:12,:);

            % These equations will isolate any angular dynamics using the difference between two velocities.
            ACombinedAnsatz((pointInd-1)*6+1,3*N+1:end) = reshape(cross(uAngFront - uAngBack, repmat(et(s),1,N)),1,[]);
            ACombinedAnsatz((pointInd-1)*6+2,3*N+1:end) = reshape(cross(uAngFront - uAngBack, repmat(eb(s),1,N)),1,[]);
            ACombinedAnsatz((pointInd-1)*6+3,3*N+1:end) = reshape(cross(uAngSide - uAngOtherSide, repmat(et(s),1,N)),1,[]);

            % These equations will isolate any translational dynamics using the sum of two velocities.
            ACombinedAnsatz((pointInd-1)*6+4,3*N+1:end) = reshape(cross(uAngFront + uAngBack, repmat(et(s),1,N)),1,[]);
            ACombinedAnsatz((pointInd-1)*6+5,3*N+1:end) = reshape(cross(uAngFront + uAngBack, repmat(en(s),1,N)),1,[]);
            ACombinedAnsatz((pointInd-1)*6+6,3*N+1:end) = reshape(cross(uAngFront + uAngBack, repmat(eb(s),1,N)),1,[]);
            
            textprogressbar(100 * pointInd / N);
        end
        textprogressbar('Done!')

        % Form the RHS of the linear system, which is 
        % various dot products of the velocities with et and eb.
        RCombinedAnsatz= zeros(6*N,1);
        for pointInd = 1 : N
            s = segmentMidpoints(pointInd);

            RCombinedAnsatz((pointInd-1)*6 + 1) = dot(uFront(:,pointInd) - uBack(:,pointInd), et(s));
            RCombinedAnsatz((pointInd-1)*6 + 2) = dot(uFront(:,pointInd) - uBack(:,pointInd), eb(s));
            RCombinedAnsatz((pointInd-1)*6 + 3) = dot(uSide(:,pointInd) - uOtherSide(:,pointInd), et(s));

            RCombinedAnsatz((pointInd-1)*6 + 4) = dot(uFront(:,pointInd) + uBack(:,pointInd), et(s));
            RCombinedAnsatz((pointInd-1)*6 + 5) = dot(uFront(:,pointInd) + uBack(:,pointInd), en(s));
            RCombinedAnsatz((pointInd-1)*6 + 6) = dot(uFront(:,pointInd) + uBack(:,pointInd), eb(s));
        end

        % Solve the linear system for the torque density.
       sol = ACombinedAnsatz \ RCombinedAnsatz;
        output.(name).force = reshape(sol(1:3*N),3,[]);
        output.(name).torque = reshape(sol(3*N+1:end),3,[]);
        interps.(name).force = griddedInterpolant(segmentMidpoints, transpose(output.(name).force), 'nearest', 'nearest');
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end

    if methodsToDo.combinedAnsatzBCApprox
        name = 'combinedAnsatzBCApprox';
        %% Combined ansatz with leading-order rotlet integral.
        % We form an 6Nx6N system and enforce the boundary condition at 2N 
        % points along the length of the slender body, using the surface 
        % velocity at XFront and XSide to give a full rank system. The
        % variables will be ordered as [f1x,f1y,...,fNz,m1x,m1y,...mNz], whilst
        % the equations will be [uFront1,uSide1,uFront2,...,uSideN].
        disp('Combined ansatz with leading-order rotlet integral...')    
        output.(name) = struct();
        ACombinedAnsatzBCApprox = zeros(6*N);
        textprogressbar('Building matrix: ')
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
            ACombinedAnsatzBCApprox((pointInd-1)*6+1,1:3*N) = et(s)' * (uTransFront - uTransBack);
            ACombinedAnsatzBCApprox((pointInd-1)*6+2,1:3*N) = eb(s)' * (uTransFront - uTransBack);
            ACombinedAnsatzBCApprox((pointInd-1)*6+3,1:3*N) = et(s)' * (uTransSide - uTransOtherSide);
            
            % These equations will isolate any translational dynamics using the sum of two velocities.
            ACombinedAnsatzBCApprox((pointInd-1)*6+4,1:3*N) = et(s)' * (uTransFront + uTransBack);
            ACombinedAnsatzBCApprox((pointInd-1)*6+5,1:3*N) = en(s)' * (uTransFront + uTransBack);
            ACombinedAnsatzBCApprox((pointInd-1)*6+6,1:3*N) = eb(s)' * (uTransFront + uTransBack);

            %----
            % Evaluate the rotlet part.
            %----
            % We'll use the leading-order result for the rotlet integrals,
            % which include only the contribution of the local torque.
            uAngFront = zeros(3,N);
            uAngBack= zeros(3,N);
            uAngSide = zeros(3,N);
            uAngOtherSide = zeros(3,N);

            uAngFront(:,pointInd) = ep * eta(s) * I1(s) * erFrontAtSegmentMidpoints(:,pointInd) - I2(s) * et(s);
            uAngBack(:,pointInd) = ep * eta(s) * I1(s) * erBackAtSegmentMidpoints(:,pointInd) - I2(s) * et(s);
            uAngSide(:,pointInd) = ep * eta(s) * I1(s) * erSideAtSegmentMidpoints(:,pointInd) - I2(s) * et(s);
            uAngOtherSide(:,pointInd) = ep * eta(s) * I1(s) * erOtherSideAtSegmentMidpoints(:,pointInd) - I2(s) * et(s);

            % These equations will isolate any angular dynamics using the difference between two velocities.
            ACombinedAnsatzBCApprox((pointInd-1)*6+1,3*N+1:end) = reshape(cross(uAngFront - uAngBack, repmat(et(s),1,N)),1,[]);
            ACombinedAnsatzBCApprox((pointInd-1)*6+2,3*N+1:end) = reshape(cross(uAngFront - uAngBack, repmat(eb(s),1,N)),1,[]);
            ACombinedAnsatzBCApprox((pointInd-1)*6+3,3*N+1:end) = reshape(cross(uAngSide - uAngOtherSide, repmat(et(s),1,N)),1,[]);

            % These equations will isolate any translational dynamics using the sum of two velocities.
            ACombinedAnsatzBCApprox((pointInd-1)*6+4,3*N+1:end) = reshape(cross(uAngFront + uAngBack, repmat(et(s),1,N)),1,[]);
            ACombinedAnsatzBCApprox((pointInd-1)*6+5,3*N+1:end) = reshape(cross(uAngFront + uAngBack, repmat(en(s),1,N)),1,[]);
            ACombinedAnsatzBCApprox((pointInd-1)*6+6,3*N+1:end) = reshape(cross(uAngFront + uAngBack, repmat(eb(s),1,N)),1,[]);

            textprogressbar(100 * pointInd / N)
        end
        textprogressbar('Done!')

        % Form the RHS of the linear system, which is 
        % various dot products of the velocities with et and eb.
        RCombinedAnsatzBCApprox = zeros(6*N,1);
        for pointInd = 1 : N
            s = segmentMidpoints(pointInd);

            RCombinedAnsatzBCApprox((pointInd-1)*6 + 1) = dot(uFront(:,pointInd) - uBack(:,pointInd), et(s));
            RCombinedAnsatzBCApprox((pointInd-1)*6 + 2) = dot(uFront(:,pointInd) - uBack(:,pointInd), eb(s));
            RCombinedAnsatzBCApprox((pointInd-1)*6 + 3) = dot(uSide(:,pointInd) - uOtherSide(:,pointInd), et(s));

            RCombinedAnsatzBCApprox((pointInd-1)*6 + 4) = dot(uFront(:,pointInd) + uBack(:,pointInd), et(s));
            RCombinedAnsatzBCApprox((pointInd-1)*6 + 5) = dot(uFront(:,pointInd) + uBack(:,pointInd), en(s));
            RCombinedAnsatzBCApprox((pointInd-1)*6 + 6) = dot(uFront(:,pointInd) + uBack(:,pointInd), eb(s));
        end

        % Solve the linear system for the torque density.
        sol = ACombinedAnsatzBCApprox \ RCombinedAnsatzBCApprox;
        output.(name).force = reshape(sol(1:3*N),3,[]);
        output.(name).torque = reshape(sol(3*N+1:end),3,[]);
        interps.(name).force = griddedInterpolant(segmentMidpoints, transpose(output.(name).force), 'nearest', 'nearest');
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end
    

    if methodsToDo.combinedAnsatzRTT
        name = 'combinedAnsatzRTT';
        %% Combined ansatz with RTT.
        % We form an 6Nx6N system and enforce the boundary condition at N 
        % points along the length of the slender body, using the full boundary 
        % condition at XFront, and the RTT relations at these N point to give a 
        % full rank system. The variables will be ordered as 
        % [f1x,f1y,...,fNz,m1x,m1y,...mNz], whilst the equations will be 
        % [uFront1,...uFrontN,RTT1,...,RTTN].
        disp('Combined ansatz with RTT...')    
        output.(name) = struct();
        ACombinedAnsatzRTT = zeros(6*N);
        textprogressbar('Building matrix: ')
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
            % ACombinedAnsatzRTT((pointInd-1)*3+1 : (pointInd-1)*3+3, 1:3*N) = reshape(integrals,3,3*N);
            ACombinedAnsatzRTT((pointInd-1)*3+1 : (pointInd-1)*3+3, 1:3*N) = reshape(diff(sol)',3,3*N);

            %----
            % Evaluate the rotlet part.
            %----
            % We'll want to integrate m x rotlet over each segment of the
            % centreline, where m is constant on each segment. To do this,
            % we'll integrate the rotlet kernel over each segment, then assign
            % the relevant components to the matrix to perform the cross
            % product.
            integrand = @(s) rotlet(XFront(:,pointInd), xi(s));

            % We'll perform the integrals over different segments in parallel,
            % as they are quite expensive.
            % integrals = zeros(3,N);
            % parfor segInd = 1 : N
            %     integrals(:,segInd) = integral(integrand, segmentEndpoints(segInd), segmentEndpoints(segInd+1), 'ArrayValued', true, 'AbsTol', tol);
            % end
            [~, sol] = ode15s(@(t,y) integrand(t), segmentEndpoints, zeros(3,1), opts);
            integrals = diff(sol)';

            % Assign the computed integrals to the linear system.
            ACombinedAnsatzRTT((pointInd-1)*3+1 : (pointInd-1)*3+3, 3*N+1:end) = cell2mat(arrayfun(@(i) crossProductMatrix(integrals(:,i)), 1:N, 'UniformOutput',false));
            textprogressbar(100 * pointInd / N)
        end
        textprogressbar('Done!')

        %----
        % Evaluate the RTT part all at once.
        %----
        % We'll use the full RTT relation at the N points.
        ACombinedAnsatzRTT(3*N+1:6*N,3*N+1:6*N) = kron(diag(I1(segmentMidpoints)),eye(3));

        % Form the RHS of the linear system for the combined ansatz, which is 
        % the velocity on the surface at XFront and XSide. 
        RCombinedAnsatzRTT = [uFront(:); Omega(:)];

        % Solve the linear system for the torque density.
        sol = ACombinedAnsatzRTT\ RCombinedAnsatzRTT;
        output.(name).force = reshape(sol(1:3*N),3,[]);
        output.(name).torque = reshape(sol(3*N+1:end),3,[]);
        interps.(name).force = griddedInterpolant(segmentMidpoints, transpose(output.(name).force), 'nearest', 'nearest');
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end


    if methodsToDo.combinedAnsatzRTTBCApprox
        name = 'combinedAnsatzRTTBCApprox';
        %% Combined ansatz with RTT, using leading-order rotlet integrals.
        % We form an 6Nx6N system and enforce the boundary condition at N 
        % points along the length of the slender body, using the leading-order
        % boundary condition at XFront, and the RTT result at these N points to
        % give a full rank system. The variables will be ordered as 
        % [f1x,f1y,...,fNz,m1x,m1y,...mNz], whilst the equations will be 
        % [uFront1,...uFrontN,RTT1,...,RTTN].
        disp('Combined ansatz with RTT and leading-order rotlet integral...')    
        output.(name) = struct();
        ACombinedAnsatzRTTBCApprox = zeros(6*N);
        textprogressbar('Building matrix: ')
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
            % ACombinedAnsatzRTTBCApprox((pointInd-1)*3+1 : (pointInd-1)*3+3, 1:3*N) = reshape(integrals,3,3*N);
            ACombinedAnsatzRTTBCApprox((pointInd-1)*3+1 : (pointInd-1)*3+3, 1:3*N) = reshape(diff(sol)',3,3*N);


            %----
            % Evaluate the rotlet part.
            %----
            % We'll use the leading-order result for the rotlet integrals,
            % which include only the contribution of the local torque.
            s = segmentMidpoints(pointInd);
            ACombinedAnsatzRTTBCApprox((pointInd-1)*3+1 : (pointInd-1)*3+3, 3*N + (pointInd-1)*3+1 : 3*N + (pointInd-1)*3+3) = crossProductMatrix(ep * eta(s) * I1(s) * erFrontAtSegmentMidpoints(:,pointInd) - I2(s) * et(s));
            textprogressbar(100 * pointInd / N)
        end

        %----
        % Evaluate the RTT part all at once.
        %----
        % We'll use the full RTT relation at the N points.
        ACombinedAnsatzRTTBCApprox(3*N+1:6*N,3*N+1:6*N) = kron(diag(I1(segmentMidpoints)),eye(3));

        % Form the RHS of the linear system for the combined ansatz, which is 
        % the velocity on the surface at XFront and XSide. 
        RCombinedAnsatzRTTBCApprox = [uFront(:); Omega(:)];

        % Solve the linear system for the torque density.
        sol = ACombinedAnsatzRTTBCApprox\ RCombinedAnsatzRTTBCApprox;
        output.(name).force = reshape(sol(1:3*N),3,[]);
        output.(name).torque = reshape(sol(3*N+1:end),3,[]);
        interps.(name).force = griddedInterpolant(segmentMidpoints, transpose(output.(name).force), 'nearest', 'nearest');
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end

    %% Evaluate the velocity generated by the different ansaetze at sample points on the surface.
    % We've now solved all of the problems and have the forces and torques
    % from each of the solution methods. We will now use these computed
    % values in the full forms of the ansaetze to evaluate the velocity on
    % the surface of the slender body at many points.
    disp('Testing boundary velocity...')    

    % For each set of forces and/or torques we have found, evaluate the
    % full ansaetz, using a piecewise-constant interpolant and numerically
    % evaluating the integrals in full.
    
    % Define the integrands.
    rotletAnsatzIntegrand = @(x,s,torque) cross(reshape(torque(s),3,1), rotlet(x, xi(s)));
    combinedAnsatzIntegrand = @(x,s,force,torque) (regularisedStokeslet(x,xi(s),regularisationParam(s)) + dipoleWeightingFactor * (e^2 - s.^2).*regularisedPotentialDipole(x,xi(s),regularisationParam(s))) * reshape(force(s),3,1) + rotletAnsatzIntegrand(x,s,torque);
    
    fnames = fieldnames(methodsToDo);
    for methodInd = 1 : numel(fnames)
        name = fnames{methodInd};
        if methodsToDo.(name)
            disp(name)
            temp = output.(name).vel;
            if contains(name, 'rotlet')
                fun = @(pointInd, s) rotletAnsatzIntegrand(evaluationPoints(:,pointInd), s, interps.(name).torque);
            else
                fun = @(pointInd, s) combinedAnsatzIntegrand(evaluationPoints(:,pointInd), s, interps.(name).force, interps.(name).torque);
            end
            parfor pointInd = 1 : numEvalPoints
                [~, sol] = ode15s(@(t,y) fun(pointInd, t), [-e,e], zeros(3,1), opts);
                temp(:,pointInd) = sol(end,:);
            end
            output.(name).vel = temp;
        end
    end

    toc

end