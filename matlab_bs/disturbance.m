% INPUT TLE HERE
tle_COMSATBW1 = "1 35943U 09054B   26159.14836492  .00000040  00000-0  00000+0 0  9995 2 35943   0.0583  75.8570 0003010   6.2458 290.8552  1.00271790 61211";

% tle to keplerian orbits 
function [a, e, i, RAAN, arg_p, true_anomaly] = TLEtoKeplerian(tle)
    tle = split(tle);
    i = str2double(tle(12)); 
    RAAN = str2double(tle(13));
    e = str2double("0." + tle(14));  % the decimal is implied (for no reason <3)
    arg_p = str2double(tle(15));
    
    % semi-major axis
    mu_Earth = 3.986004418e14;
    n = str2double(tle(17)) * (2*pi/86400);
    a = (mu_Earth/n^2)^(1/3);
    
    % true anomaly
    M = deg2rad(str2double(tle(16)));
    E = M;
    E_prev = 0;  % initialisation value
    while abs(E - E_prev) > 1e-6
        E = E + (M - E + e*sin(E)) / (1 - e*cos(E));
        E_prev = E;
    end
    true_anomaly = 2 * atan2(sqrt(1 + e) * sin(E/2), sqrt(1 - e) * cos(E/2));
end

% block parameters >> orbit >> initial state format >> orbital elements
[a, e, i, RAAN, arg_p, true_anomaly] = TLEtoKeplerian(tle_COMSATBW1);


% block parameters >> orbit >> initial state format >> icrf state vector
position_velocity_matrices;
ts_pos = timeseries(position*10^3, time_elapsed);
ts_vel = timeseries(velocity*10^3, time_elapsed);