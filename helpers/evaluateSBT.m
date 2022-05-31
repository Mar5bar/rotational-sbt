function [output, As] = evaluateSBT(params, methods, As)
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

    if nargin < 3
        buildAs = true;
        As = struct();
    elseif isempty(As)
        buildAs = true;
        As = struct();
    else
        buildAs = false;
        disp("Using provided linear systems.")
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
    localBasis = @(s) [et(s), en(s), eb(s)];
    
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
    % XOtherSide = xi - epsilon * eta * er(s,pi/2), and
    ers = struct();
    ers.front = cell2mat(arrayfun(@(s) er(s,0), segmentMidpoints, 'UniformOutput', false));
    ers.side = cell2mat(arrayfun(@(s) er(s,pi/2), segmentMidpoints, 'UniformOutput', false));
    ers.back = cell2mat(arrayfun(@(s) er(s,pi), segmentMidpoints, 'UniformOutput', false));
    ers.otherSide = cell2mat(arrayfun(@(s) er(s,3*pi/2), segmentMidpoints, 'UniformOutput', false));
    
    XFront = xi(segmentMidpoints) + ep * eta(segmentMidpoints) .* ers.front;
    XSide = xi(segmentMidpoints) + ep * eta(segmentMidpoints) .* ers.side;
    XBack = xi(segmentMidpoints) + ep * eta(segmentMidpoints) .* ers.back;
    XOtherSide = xi(segmentMidpoints) + ep * eta(segmentMidpoints) .* ers.otherSide;

    % Now evaluate the velocity at these points.
    [uFront, uSide, uBack, uOtherSide, Omega] = genVelBC(linVel, angVel, segmentMidpoints, eta, ers, params, localBasis);

    % Generate the final evaluation points and evaluate the prescribed velocity 
    % at these points. We won't sample directly at the tips of the slender 
    % body.
    [evaluationPoints, velPrescribed, evaluationArclengths, evaluationAngles] = genEvalPoints(linVel, angVel, params, xi, eta, er, localBasis);
    numEvalPoints = size(evaluationPoints,2);

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
    output.evaluationArclengths = evaluationArclengths;
    output.evaluationAngles = evaluationAngles;

    %% Options for computing integrals.
    opts = odeset('AbsTol',tol,'RelTol',tol);

    %% Save intermediates in a struct for easy passing to external functions.
    intermediates = struct();
    intermediates.N = N;
    intermediates.ep = ep;
    intermediates.segmentEndpoints = segmentEndpoints;
    intermediates.segmentMidpoints = segmentMidpoints;
    intermediates.XFront = XFront;
    intermediates.XSide = XSide;
    intermediates.XBack = XBack;
    intermediates.XOtherSide = XOtherSide;
    intermediates.uFront = uFront;
    intermediates.uSide = uSide;
    intermediates.uBack = uBack;
    intermediates.uOtherSide = uOtherSide;
    intermediates.opts = opts;
    intermediates.xi = xi;
    intermediates.eta = eta;
    intermediates.et = et;
    intermediates.en = en;
    intermediates.eb = eb;
    intermediates.ers = ers;
    intermediates.Omega = Omega;
    intermediates.I1 = I1;
    intermediates.I1Approx = I1Approx;
    intermediates.I2 = I2;
    intermediates.dipoleWeightingFactor = dipoleWeightingFactor;
    intermediates.regularisationParam = regularisationParam;
    intermediates.e = e;

    % We'll form interpolants of all computed forces and torques.
    interps = struct();

    % We'll store the computed linear systems and RHSs.
    % Only rebuild the linear systems if not passed.
    % We'll always rebuild the RHS.
    Rs = struct();

    %% Form and solve the linear systems.
    if methodsToDo.rotletAnsatz
        %% Full rotlet ansatz.
        name = 'rotletAnsatz';
        output.(name) = struct();
        if buildAs
            As.(name) = genARotletAnsatz(intermediates);
        end
        Rs.(name) = genRRotletAnsatz(intermediates);

        % Solve the linear system for the torque density.
        sol = As.(name) \ Rs.(name);
        output.(name).torque = reshape(sol,3,[]);
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end

    if methodsToDo.rotletAnsatzRTT
        %% Rotlet ansatz using RTT.
        name = 'rotletAnsatzRTT';
        output.(name) = struct();
        
        if buildAs
            As.(name) = genARotletAnsatzRTT(intermediates);
        end
        Rs.(name) = genRRotletAnsatzRTT(intermediates);

        % Solve the linear system for the torque density.
        sol = As.(name) \ Rs.(name);
        output.(name).torque = reshape(sol,3,[]);
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end

    if methodsToDo.rotletAnsatzRTTI1Approx
        %% Rotlet ansatz using approximate RTT.
        name = 'rotletAnsatzRTTI1Approx';
        output.(name) = struct();
        if buildAs
            As.(name) = genARotletAnsatzRTTI1Approx(intermediates);
        end
        Rs.(name) = genRRotletAnsatzRTTI1Approx(intermediates);

        % Solve the linear system for the torque density.
        sol = As.(name) \ Rs.(name);
        output.(name).torque = reshape(sol,3,[]);
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end

    if methodsToDo.combinedAnsatz
        %% Full combined ansatz.
        name = 'combinedAnsatz';
        output.(name) = struct();
        if buildAs
            As.(name) = genACombinedAnsatz(intermediates);
        end
        Rs.(name) = genRCombinedAnsatz(intermediates);

        % Solve the linear system for the torque density.
        sol = As.(name) \ Rs.(name);
        output.(name).force = reshape(sol(1:3*N),3,[]);
        output.(name).torque = reshape(sol(3*N+1:end),3,[]);
        interps.(name).force = griddedInterpolant(segmentMidpoints, transpose(output.(name).force), 'nearest', 'nearest');
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end

    if methodsToDo.combinedAnsatzBCApprox
        %% Combined ansatz with leading-order rotlet integral.
        name = 'combinedAnsatzBCApprox';
        output.(name) = struct();
        if buildAs
            As.(name) = genACombinedAnsatzBCApprox(intermediates);
        end
        Rs.(name) = genRCombinedAnsatzBCApprox(intermediates);

        % Solve the linear system for the torque density.
        sol = As.(name) \ Rs.(name);
        output.(name).force = reshape(sol(1:3*N),3,[]);
        output.(name).torque = reshape(sol(3*N+1:end),3,[]);
        interps.(name).force = griddedInterpolant(segmentMidpoints, transpose(output.(name).force), 'nearest', 'nearest');
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end
    

    if methodsToDo.combinedAnsatzRTT
        %% Combined ansatz with RTT.
        name = 'combinedAnsatzRTT';
        output.(name) = struct();
        if buildAs
            As.(name) = genACombinedAnsatzRTT(intermediates);
        end
        Rs.(name) = genRCombinedAnsatzRTT(intermediates);

        % Solve the linear system for the torque density.
        sol = As.(name)\ Rs.(name);
        output.(name).force = reshape(sol(1:3*N),3,[]);
        output.(name).torque = reshape(sol(3*N+1:end),3,[]);
        interps.(name).force = griddedInterpolant(segmentMidpoints, transpose(output.(name).force), 'nearest', 'nearest');
        interps.(name).torque = griddedInterpolant(segmentMidpoints, transpose(output.(name).torque), 'nearest', 'nearest');
        output.(name).vel = zeros(3,numEvalPoints);
    end


    if methodsToDo.combinedAnsatzRTTBCApprox
        %% Combined ansatz with RTT, using leading-order rotlet integrals.
        name = 'combinedAnsatzRTTBCApprox';
        output.(name) = struct();
        if buildAs
            As.(name) = genACombinedAnsatzRTTBCApprox(intermediates);
        end
        Rs.(name) = genRCombinedAnsatzRTTBCApprox(intermediates);

        % Solve the linear system for the torque density.
        sol = As.(name) \ Rs.(name);
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