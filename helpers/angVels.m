function [angVel, angVelLab] = angVels(selection)
%% angVels returns a function handle angVel and a label angVelLab corresponding to the
% integer selection made. The components of angVel are written in the et, en, eb basis.
    switch selection
        case 1
            angVel = @(s) [1;0;0] + 0.*s; angVelLab = 'et';
        case 2
            angVel = @(s) [0;1;0] + 0.*s; angVelLab = 'en';
        case 3
            angVel = @(s) [1;0;0].*reshape(s,1,[]); angVelLab = 's*et';
        otherwise
            angVel = @(s) selection*[1;0;0] + 0.*s; angVelLab = [num2str(selection),' et'];
    end
end