% these values will be taken from the CATIA model
Ixx = 1800;
Ixy = 50;
Ixz = 30;
Iyy = 2200;
Iyz = 40;
Izz = 2600;

I = [Ixx, -Ixy, -Ixz;
    -Ixy, Iyy, -Iyz;
    -Ixz, -Iyz, Izz];

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
R = eye(3);
K = lqr(A,B,Q,R);
