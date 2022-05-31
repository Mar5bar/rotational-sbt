function [evaluationPoints, velPrescribed, evaluationArclengths, evaluationAngles] = genEvalPoints(linVel, angVel, params, xi, eta, er, localBasis)
    ep = params.epsilon;
    evaluationArclengths = linspace(-1, 1, params.numArclengthEvaluationPoints + 2);
    evaluationArclengths = evaluationArclengths(2 : end-1);
    evaluationAngles = linspace(0, 2*pi, params.numCircumferentialEvaluationPoints + 1);
    evaluationAngles = evaluationAngles(1 : end-1);
    numEvalPoints = length(evaluationArclengths)*length(evaluationAngles);
    evaluationPoints = zeros(3,numEvalPoints);
    velPrescribed = zeros(3,numEvalPoints);
    for i = 1 : length(evaluationArclengths)
        s = evaluationArclengths(i);
        ers = cell2mat(arrayfun(@(phi) er(s, phi), evaluationAngles, 'UniformOutput', false));
        evaluationPoints(:,(i-1)*length(evaluationAngles)+1 : i*length(evaluationAngles)) = xi(s) + ep * eta(s) * ers;
        velPrescribed(:,(i-1)*length(evaluationAngles)+1 : i*length(evaluationAngles)) = linVel(s) + ep * eta(s) * cross(repmat(localBasis(s)*angVel(s),1,length(evaluationAngles)), ers);
    end
end