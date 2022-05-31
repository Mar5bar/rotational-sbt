function A = genARotletAnsatz(in,verbose)
% We form an 3Nx3N system and enforce the boundary condition at N 
% points along the length of the slender body, using the surface 
% velocity at XFront and XSide to give a full rank system. 
    if nargin < 2
        verbose = false;
    end
    if verbose
        disp('Rotlet ansatz...')
    end
    if verbose
        textprogressbar('Building matrix: ')
    end
    % Unpack the struct.
    N = in.N;
    segmentMidpoints = in.segmentMidpoints;
    segmentEndpoints = in.segmentEndpoints;
    XFront = in.XFront;
    XSide = in.XSide;
    xi = in.xi;
    et = in.et;
    eb = in.eb;
    opts = in.opts;
    
    A = zeros(3*N);
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
        A((pointInd-1)*3+1,:) = reshape(cross(integrals(1:3,:), repmat(et(s), 1, N)),[],1);

        % Second equation, computing the cross products of the XFront integrals
        % with eb(segmendMidpoints(pointInd)).
        A((pointInd-1)*3+2,:) = reshape(cross(integrals(1:3,:), repmat(eb(s), 1, N)),[],1);

        % Third equation, computing the cross products of the XSide integrals
        % with et(segmendMidpoints(pointInd)).
        A((pointInd-1)*3+3,:) = reshape(cross(integrals(4:6,:), repmat(et(s), 1, N)),[],1);
        if verbose
            textprogressbar(100 * pointInd / N)
        end
    end
    if verbose
        textprogressbar('Done!')
    end
end