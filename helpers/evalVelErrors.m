function errors = evalVelErrors(name, output, NameValueArgs)
    arguments
        name char
        output struct
        NameValueArgs.structure struct
    end

    if isfield(NameValueArgs, 'structure')
        errors = NameValueArgs.structure;
    else
        errors = struct();
    end

    vel = output.(name).vel;

    errors.absVel = abs(vel - output.velPrescribed);
    errors.relVel =  abs((vel - output.velPrescribed) ./ output.velPrescribed);
    errors.maxAbsVel = max(max(errors.absVel));
    errors.maxRelVel = max(max(errors.relVel(~isinf(errors.relVel) & abs(output.velPrescribed) > 1e-12)));
    
end