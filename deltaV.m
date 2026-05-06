% inputs for each orbit: altitude (assuming circular orbit), inclination,
% RAAN

% r_1 = 35786;  % GEO [km]

% orbit1 = (r_1, i_1, RAAN_1);
% orbit2 = (r_2, i_2, RAAN_2);

% function V = deltaV(orbit1, orbit2)

% 
% classdef deltaV
%     properties
%         mu_Earth = 398600.4415;  % [km^2 / s^2]
%         r_Earth = 6378.137;  % [km]
%     end
%     methods
%         function V = velocity(r)
%             V = sqrt(obj.mu_Earth / r);
%         end
%         function cos_theta = cosTheta(i_1, i_2, RAAN_1, RAAN_2)
%             cos_theta = cos(i_1) * cos(i_2) + sin(i_1) * sin(i_2) * cos(RAAN_2 - RAAN_1);
%         end
% 
% 
%         function dV = transferV(r_1, r_2)
%             vt1 = sqrt(obj.mu_Earth * (2/r_1 - 2/(r_1 + r_2) ))
%             vt2 = sqrt(obj.mu_Earth * (2/r_2 - 2/(r_1 + r_2) ))
%         end
%     end
% end


mu_Earth = 398600.4415;  % [km^2 / s^2]
r_Earth = 6378.137;  % [km]

% input parameters
h_1 = 35700;
h_2 = 36500;
r_1 = h_1 + r_Earth;  % [km]
r_2 = h_2 + r_Earth;
i_1 = 20 * pi/180;  % [rad]
i_2 = 20 * pi/180;
RAAN_1 = 0 * pi/180;
RAAN_2 = 180 * pi/180;  % [rad]

cos_theta = cos(i_1) * cos(i_2) + sin(i_1) * sin(i_2) * cos(RAAN_2 - RAAN_1);
theta = acos(cos_theta);

v1 = sqrt(mu_Earth / r_1);
v2 = sqrt(mu_Earth / r_2);

vt1 = sqrt(mu_Earth * (2/r_1 - 2/(r_1 + r_2) ));
vt2 = sqrt(mu_Earth * (2/r_2 - 2/(r_1 + r_2) ));

syms theta1

eqn = v1 * vt1 * sin(theta1)/(sqrt(v1^2 + vt1^2 - 2 * v1 *vt1 * cos(theta1))) == v2 * vt2 * sin(theta - theta1)/(sqrt(v2^2 + vt2^2 - 2 * v2 *vt2 * cos(theta - theta1)));
solution = solve(eqn, theta1)

% deltaV = (sqrt(v1^2 + vt1^2 - 2 * v1 *vt1 * cos(theta1))