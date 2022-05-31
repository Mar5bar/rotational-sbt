function [linVel, linVelLab] = linVels(selection)
%% linVels returns a function handle linVel and a label linVelLab corresponding to the
% integer selection made.
    switch selection
        case 1
            linVel = @(s) [0;0;0] + 0.*s; linVelLab = '0';
        case 2
            linVel = @(s) [1;1;1]/sqrt(3) + 0.*s; linVelLab = '[1;1;1] / sqrt(3)';
        case 3
            linVel = @(s) 10*[1;1;1]/sqrt(3) + 0.*s; linVelLab = '10[1;1;1] / sqrt(3)';
    end
end