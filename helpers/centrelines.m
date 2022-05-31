function [x, y, z, dx, dy, dz, tRange, centrelineLab] = centrelines(selection)
%% CENTRELINES returns function handles x,y,z, parameter range tRange, and a
% label centrelineLab corresponding to the integer selection made.
    switch selection
        case 1
            x = @(t) t; dx = @(t) 1 + 0*t;
            y = @(t) 0*t; dy = @(t) 0*t;
            z = @(t) 0*t; dz = @(t) 0*t;
            tRange = [-1,1];
            centrelineLab = '(t,0,0)';
        case 2
            x = @(t) t; dx = @(t) 1 + 0*t;
            y = @(t) t.^2; dy = @(t) 2*t;
            z = @(t) 0*t; dz = @(t) 0*t;
            tRange = [-1,1];
            centrelineLab = '(t,t^2,0)';
        case 3
            x = @(t) sin(t); dx = @(t) cos(t);
            y = @(t) sin(t).*cos(t); dy = @(t) cos(2*t);
            z = @(t) 0*t; dz = @(t) 0*t;
            tRange = [pi/2,3*pi/2];
            centrelineLab = '(sin(t),sin(t)cos(t),0)';
        case 4
            x = @(t) t.^0.25.*sin(t); dx = @(t) t.^0.25.*cos(t) + 0.25*t.^(-0.75).*sin(t);
            y = @(t) t.^0.25.*cos(t); dy = @(t) -t.^0.25.*sin(t) + 0.25*t.^(-0.75).*cos(t);
            z = @(t) 0*t; dz = @(t) 0*t;
            tRange = [0.1,2*pi];
            centrelineLab = '(t^0.25sin(t),t^0.25cos(t),0)';
        case 5
            x = @(t) t; dx = @(t) 1 + 0*t;
            y = @(t) tanh(100*t); dy = @(t) 100*sech(100*t).^2;
            z = @(t) 0*t; dz = @(t) 0*t;
            tRange = [-0.5,0.5];
            centrelineLab = '(t,tanh(100t),0)';
        case 6
            x = @(t) 16*sin(t).^3; dx = @(t) 48*sin(t).^2.*cos(t);
            y = @(t) 13*cos(t)-5*cos(2*t)-2*cos(3*t)-cos(4*t); dy = @(t) -13*sin(t) + 10*sin(2*t) + 6*sin(3*t) +4*sin(4*t);
            z = @(t) 0*t; dz = @(t) 0*t;
            tRange = [pi/2,3*pi/2]; 
            centrelineLab = '(16sin(t)^3,13cos(t)-5cos(2t)-2cos(3t)-cos(4t),0)';
        case 7
            x = @(t) sin(t); dx = @(t) cos(t);
            y = @(t) cos(t); dy = @(t) -sin(t);
            z = @(t) t; dz = @(t) 1 + 0*t;
            tRange = [0,2*pi];
            centrelineLab = '(sin(t),cos(t),t)';
        case 8
            x = @(t) t; dx = @(t) 1 + 0*t;
            y = @(t) 0.3*cos(pi*t-0.2); dy = @(t) -0.3*pi*sin(pi*t-0.2);
            z = @(t) 0*t; dz = @(t) 0*t;
            tRange = [0,2];
            centrelineLab = '(t,0.3cos(tpi-0.2),0)';
    end
end