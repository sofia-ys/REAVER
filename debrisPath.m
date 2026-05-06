% r_Earth = 6378.137;  % [km]
% 
% % recycling hub orbital parameters
% h_0 = 37586;  % [km]
% i_0 = 7;  % [deg]
% RAAN_0 = 0;  % [deg]
% rh_orbit = [h_0 + r_Earth, i_0 * 180/pi, RAAN_0 * 180/pi];
% 
% % debris 1 orbital parameters
% h_1 = 3570;  % [km]
% i_1 = 20;  % [deg]
% RAAN_1 = 0;  % [deg]
% d1_orbit = [h_1 + r_Earth, i_1 * 180/pi, RAAN_1 * 180/pi];
% 
% % debris 2 orbital parameters
% h_2 = 3650;  % [km]
% i_2 = 10;  % [deg]
% RAAN_2 = 0;  % [deg]
% d2_orbit = [h_2 + r_Earth, i_2 * 180/pi, RAAN_2 * 180/pi];
% 
% % debris 3 orbital parameters
% h_3 = 3630;  % [km]
% i_3 = 10;  % [deg]
% RAAN_3 = 10;  % [deg]
% d3_orbit = [h_3 + r_Earth, i_3 * 180/pi, RAAN_3 * 180/pi];
% 
% orbits = [
%     rh_orbit;
%     d1_orbit;
%     d2_orbit;
%     d3_orbit
% ];
% 
% n_targets = 3;
% debris_perms = perms((2:n_targets+1));
% 
% orbits_list = cell(size(debris_perms, 1), 1);
% 
% for k = 1:size(debris_perms, 1)
%     debris_perm = debris_perms(k, :);
% 
%     % start at RH, visit debris, return to RH
%     path_indices = [1, debris_perm, 1];
% 
%     % store actual orbital data for this path
%     orbits_list{k} = orbits(path_indices, :);
% end
% 
% dvCalc = deltaV();


dvCalc = deltaV();
r_1 = 7000;
r_2 = 42164;
i_1 = deg2rad(28.5);
i_2 = deg2rad(0);
RAAN_1 = deg2rad(0);
RAAN_2 = deg2rad(20);
dV_total = dvCalc.delV(r_1, r_2, i_1, i_2, RAAN_1, RAAN_2)

% can assume any RN for RH 
% circular orbits is fine 
% why would tugs be 500kg each? -> tugs have electric prop