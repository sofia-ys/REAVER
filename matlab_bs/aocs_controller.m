% these values will be taken from the CATIA model
Ixx = 1800;
Ixy = 50;
Ixz = 30;
Iyy = 2200;
Iyz = 40;
Izz = 2600;

I = [Ixx, -Ixy, -Ixz;
    -Ixy, Iyy, -Iyz;
    -Ixz, -Iyz, Izz]

% defining the state matrices
A = [0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1;
     zeros(3,6)];

B = [zeros(3,3);
     inv(I)];

C = eye(6, 6);

D = zeros(6, 3);

% state-space system
satellite = ss(A, B, C, D);

% controller
Q = eye(6);
R_c = eye(3);
K = lqr(A,B,Q,R_c);

C_i = C(1:3,:); % augmenting C to integrate attitude angle error (3 integrators)
n = size(A,1);       
m = size(B,2);       
p_i = size(C_i,1);
Ai = [A, zeros(n,p_i);  % augmented A and B matrices 
      -C_i, zeros(p_i,p_i)];
Bi = [B;
      zeros(p_i,m)];

% Qx and Qz needs to be tuned 
Qx = diag([100,100,100,1,1,1]);  % states that are penalised in controller logic (angles that are more importnat to control)
Qz = diag([1000,1000,1000]);  % high value reduces steady-state integral error faster, but too much =overshoot
Q_i = blkdiag(Qx, Qz);  % adds these two matrices together for the diagonal
R = eye(m); 

% gain matrices
K_aug = lqr(Ai, Bi, Q_i, R);  
Kx = K_aug(:,1:n);
Ki_z = K_aug(:, n+1:end);         



% modelling disturbances
r_geo = 42164e3;
n = sqrt(3.986004418e14 / r_geo^3);   % rad/s (~7.292e-5)
t = (0:60:24*3600)';                  % sample every 60 s for 1 day
phi0 = 0;                             % initial phase (rad)

pos = [ r_geo*cos(n*t + phi0), ...
    r_geo*sin(n*t + phi0), ...
    zeros(size(t)) ];             % Nx3 matrix, units meters

% To use in Simulink: create timeseries
pos_ts = timeseries(pos, t);

a = [0 1 2;
    0 2 2;
    3 1 1];

disp("-----------------")
for el = 1:length(a)
    disp(a(el, :))
end