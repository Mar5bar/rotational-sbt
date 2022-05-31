function errors = evalForceTorqueErrors(name, output, NameValueArgs)
    arguments
        name char
        output struct
        NameValueArgs.trueTorque double
        NameValueArgs.trueForce double
        NameValueArgs.structure struct
    end

    if isfield(NameValueArgs, 'structure')
        errors = NameValueArgs.structure;
    else
        errors = struct();
    end

    torque = output.(name).torque;
    trueTorque = NameValueArgs.trueTorque;
    errors.absTorque = abs(torque - trueTorque);
    errors.relTorque =  abs((torque - trueTorque) ./ trueTorque);
    errors.maxAbsTorque = max(max(errors.absTorque));
    errors.maxRelTorque = max(max(errors.relTorque(~isinf(errors.relTorque) & abs(trueTorque) > 1e-12)));
    
    if contains(name, 'combined')
        force = output.(name).force;
        trueForce = NameValueArgs.trueForce;
        errors.absForce = abs(force - trueForce);
        errors.relForce =  abs((force - trueForce) ./ trueForce);
        errors.maxAbsForce = max(max(errors.absForce));
        errors.maxRelForce = max(max(errors.relForce(~isinf(errors.relForce) & abs(trueForce) > 1e-12)));
    end
    
end