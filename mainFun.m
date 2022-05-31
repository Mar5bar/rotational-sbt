function [output, errors] = mainFun(params)
    arguments
        params.epsilon = 1e-2   % Slenderness parameter.
        params.N = 40    % Level of numerical discretisation.
        params.xiSelector = 1 % Centreline
        params.etaSelector = 1 % Shape function
        params.linVelSelector = 1 % Linear velocity
        params.angVelSelector = 1 % Angular velocity
        params.numArclengthEvaluationPoints = 11 % Number of eval cross sections
        params.numCircumferentialEvaluationPoints = 10 % Number of eval radii
    end

    % Evaluate the SBT ansaetze.
    output = evaluateSBT(params);

   % Evaluate the error in the velocity fields.
    errors = struct();
    errors.rotletAnsatz = evalVelErrors('rotletAnsatz', output);
    errors.rotletAnsatzRTT = evalVelErrors('rotletAnsatzRTT', output);
    errors.rotletAnsatzRTTI1Approx = evalVelErrors('rotletAnsatzRTTI1Approx', output);
    errors.combinedAnsatz = evalVelErrors('combinedAnsatz',output);
    errors.combinedAnsatzBCApprox = evalVelErrors('combinedAnsatzBCApprox',output);
    errors.combinedAnsatzRTT = evalVelErrors('combinedAnsatzRTT',output);
    errors.combinedAnsatzRTTBCApprox = evalVelErrors('combinedAnsatzRTTBCApprox',output);

end