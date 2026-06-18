omega_0 = 7.272 * 10^-5;  % orbital rate of satellite wrt earth (one rev per day)

I_x = 3770;
I_y = 730;
I_z = 4020;

%% Revised design: integrate only attitude outputs (first 3 outputs)
% Define numeric plant (from your data)
Ixx = 1800;
Ixy = 50;
Ixz = 30;
Iyy = 2200;
Iyz = 40;
Izz = 2600;

I = [Ixx, -Ixy, -Ixz;
    -Ixy, Iyy, -Iyz;
    -Ixz, -Iyz, Izz];

A = [0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1;
     zeros(3,6)];

B = [zeros(3,3);
     inv(I)];   % 6x3

C = eye(6);
D = zeros(6,3);

% Choose which outputs to integrate: here integrate first 3 (attitude angles)
int_idx = 1:3;  % change if attitude outputs are different rows
Ci = C(int_idx, :);   % 3x6

% Form augmented plant for 3 integrators
na = size(A,1);    % 6
p_int = size(Ci,1);% 3
A_aug = [A, zeros(na,p_int);
         -Ci, zeros(p_int,p_int)];  % (6+3)x(6+3) = 9x9
B_aug = [B;
         zeros(p_int,size(B,2))];  % 9x3

% Check controllability
Co = ctrb(A_aug, B_aug);
rankCo = rank(Co);
fprintf('Augmented controllability rank = %d (should equal %d for full controllability)\n', rankCo, size(A_aug,1));
if rankCo < size(A_aug,1)
    warning('Augmented system is not fully controllable. LQR or place may still produce stabilizing gains on controllable subspace.');
end

% LQR weights (tune these)
Qx = 100 * eye(na);          % weight on plant states
Qi = 300 * eye(p_int);       % weight on integrator states (increase to reduce steady-state error)
Q = blkdiag(Qx, Qi);         % 9x9
R = 0.05 * eye(size(B,2));   % 3x3 (tune; larger R = less aggressive)

% Attempt LQR
try
    [K_aug, S, e] = lqr(A_aug, B_aug, Q, R);   % K_aug is 3x9
    disp('LQR succeeded. Closed-loop eigenvalues:');
    disp(e);
catch ME
    disp('LQR failed on the 9x9 augmented system:');
    disp(ME.message);
    disp('Attempting pole placement on controllable subspace...');
    % If LQR fails, we will attempt pole placement on controllable subspace.
    % Choose desired poles (example values)
    desired_poles = [-1 -1.2 -1.5 -0.8 -1.0 -1.3 -0.9 -0.95 -1.1]; % 9 poles (adjust)
    if length(desired_poles) ~= size(A_aug,1)
        error('Adjust desired_poles length to match augmented system order.');
    end
    % Try place (may fail if not controllable)
    K_aug = place(A_aug, B_aug, desired_poles);
    disp('Pole placement succeeded.');
end

% Partition K_aug into Kx and Ki (note Ki now 3x3)
Kx = K_aug(:, 1:na);        % 3x6
Ki = K_aug(:, na+1:end);    % 3x3

disp('Kx (3x6) ='); disp(Kx);
disp('Ki (3x3) ='); disp(Ki);

% Prepare to write to Simulink model 'aocs'
model = 'aocs';
if ~bdIsLoaded(model)
    load_system(model);
end

% Block paths (adjust if block names differ)
blk_Kp = [model '/Gain'];        % Gain
blk_Ki = [model '/Gain2'];       % Gain2
blk_Int = [model '/Integrator']; % Integrator

% Convert to strings
Kx_str = mat2str(Kx);
Ki_str = mat2str(Ki);

% Apply gains
set_param(blk_Kp, 'Gain', Kx_str);
% Since Ki is 3x3 now, but your current Gain2 was 3x6, we will place Ki in the Gain2 block
% and you must ensure integrator input/outputs match shape: integrator should integrate 3 errors (angles)
set_param(blk_Ki, 'Gain', Ki_str);

% Set integrator initial condition to zeros(3,1) if integrator has dimension 3
set_param(blk_Int, 'InitialCondition', 'zeros(3,1)');

save_system(model);
disp('Gains written to model. Note: Ki is 3x3 (integrator on 3 attitude outputs).');

% Remind user to verify wiring
disp('IMPORTANT: Ensure the integrator block has dimension 3 (integrating 3 signals) and it is placed on the attitude error signals (rows 1:3). If your integrator currently integrates 6 signals, you must replace it with a 3-channel integrator or rewire so Ki multiplies the 3 integrator outputs.');
