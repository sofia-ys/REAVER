% HOW TO USE:
% tldr; DO NOT TOUCH THIS
% why? this is a class file, in matlab you can ONLY define the class here (to
% run it you need a different script)

% IN A SEPARATE FILE
% 1)    create an instance of the class
%       dvCalc = deltaV();

% 2)    define the input parameters you need
%       r_1 = 7000;
%       r_2 = 42164;
%       i_1 = deg2rad(28.5);
%       i_2 = deg2rad(0);
%       RAAN_1 = deg2rad(0);
%       RAAN_2 = deg2rad(20);

% 3)    run the function you want
%       dV_total = dvCalc.delV(r_1, r_2, i_1, i_2, RAAN_1, RAAN_2)

classdef deltaV
    properties
        mu_Earth = 398600.4415;  % [km^2 / s^2]
        r_Earth = 6378.137;  % [km]
    end
    methods        
        % orbital velocity at a given radius
        function V = velocity(obj, r)
            V = sqrt(obj.mu_Earth / r);
        end
        
        % cos(theta) value from the inclination and RAAN changes
        function theta = cosTheta(obj, i_1, i_2, RAAN_1, RAAN_2)
            cos_theta = cos(i_1) * cos(i_2) + sin(i_1) * sin(i_2) * cos(RAAN_2 - RAAN_1);
            theta = acos(cos_theta);
        end

        % transfer orbit delta-v stuff
        function [vt1, vt2] = transferV(obj, r_1, r_2)
            vt1 = sqrt(obj.mu_Earth * (2/r_1 - 2/(r_1 + r_2)));
            vt2 = sqrt(obj.mu_Earth * (2/r_2 - 2/(r_1 + r_2)));
        end

        % find the theta_1 value for the delta v calc
        function theta1 = findTheta1(obj, r_1, r_2, i_1, i_2, RAAN_1, RAAN_2)
            v1 = obj.velocity(r_1);
            v2 = obj.velocity(r_2);

            [vt1, vt2] = obj.transferV(r_1, r_2);

            theta = obj.cosTheta(i_1, i_2, RAAN_1, RAAN_2);

            f = @(theta1) ...
                v1 * vt1 * sin(theta1) ./ sqrt(v1^2 + vt1^2 - 2*v1*vt1*cos(theta1)) ...
                - ...
                v2 * vt2 * sin(theta - theta1) ./ sqrt(v2^2 + vt2^2 - 2*v2*vt2*cos(theta - theta1));

            theta1 = fzero(f, theta/2);
        end

        % calculate final delta v 
        function dV = delV(obj, r_1, r_2, i_1, i_2, RAAN_1, RAAN_2)
            v1 = obj.velocity(r_1);
            v2 = obj.velocity(r_2);
            [vt1, vt2] = obj.transferV(r_1, r_2);

            theta = obj.cosTheta(i_1, i_2, RAAN_1, RAAN_2);
            theta1 = obj.findTheta1(r_1, r_2, i_1, i_2, RAAN_1, RAAN_2);

            dV1 = sqrt(v1^2 + vt1^2 - 2*v1*vt1*cos(theta1));
            dV2 = sqrt(v2^2 + vt2^2 - 2*v2*vt2*cos(theta - theta1));
            dV = dV1 + dV2;
        end
    end
end

% VERIFICATION:
% tested with pure altitude change hohmann transfer 
% tested with pure inclination change 