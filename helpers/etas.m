function [eta, etaLab] = etas(selection)
%% ETAS returns a function handle eta and a label etaLab corresponding to the
% integer selection made.
    switch selection
        case 1
            eta = @(s) sqrt(1-s.^2); etaLab = 'sqrt(1-s^2)';
        case 2
            eta = @(s) (1-s).*(s+1).*(s.^3+1).*(s.^2 + 0.2); etaLab = '(1-s)(s+1)(s^3+1)(s^2+0.2)';
        case 3
            eta = @(s) sqrt(1-s.^2).*(1.1+sin(9*pi*s)); etaLab = 'sqrt(1-s^2)(1.1+sin(9pis)';
        case 4
            eta = @(s) sqrt(1-s.^2).*(1-0.1*cos(2*pi*s)); etaLab = 'sqrt(1-s^2)(1-0.1cos(2pis)';
        case 5
            eta = @(s) (1+s).*(1-s).^2; etaLab = '(1+s)(1-s)^2';
        case 6
            eta = @(s) (1-abs(s).^0.25) .* (1 + (s+0.1).^2); etaLab = '(1-abs(s)^0.25)(1+(s+0.1)^2)';
        case 7
            eta = @(s) 1 + 0*s; etaLab = '1';
    end
end