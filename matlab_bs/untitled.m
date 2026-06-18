%% Apply Kx/Ki to your 'aocs' model (uses actual top-level paths)
model = 'aocs';

% --- Design weights (tune if needed) ---
Qx = 100 * eye(6);    % state weight
Qi = 300 * eye(3);    % integrator weight
R  = 0.05 * eye(3);   % input weight

% --- Numeric plant (user-provided) ---
Ixx = 1800; Ixy = 50; Ixz = 30; Iyy = 2200; Iyz = 40; Izz = 2600;
I = [Ixx, -Ixy, -Ixz; -Ixy, Iyy, -Iyz; -Ixz, -Iyz, Izz];

A = [0,0,0,1,0,0;
     0,0,0,0,1,0;
     0,0,0,0,0,1;
     zeros(3,6)];
B = [zeros(3,3); inv(I)];   % 6x3
C = eye(6);
D = zeros(6,3);

% --- Integrate first 3 outputs (assumed attitudes) ---
int_idx = 1:3;
Ci = C(int_idx, :);   % 3x6

% --- Augmented plant (9x9) ---
na = size(A,1);
p_int = size(Ci,1);
A_aug = [A, zeros(na,p_int); -Ci, zeros(p_int,p_int)];
B_aug = [B; zeros(p_int, size(B,2))];

% --- Check controllability ---
rankCo = rank(ctrb(A_aug, B_aug));
fprintf('Augmented controllability rank = %d of %d\n', rankCo, size(A_aug,1));
if rankCo < size(A_aug,1)
    warning('Augmented system not fully controllable. LQR may still work on controllable subspace or fallback to place will be used.');
end

% --- LQR design with fallback to place if needed ---
Q = blkdiag(Qx, Qi);
try
    [K_aug, ~, E] = lqr(A_aug, B_aug, Q, R);   % K_aug is 3x9
    disp('LQR succeeded.');
catch ME
    warning('LQR failed: %s', ME.message);
    % fallback: example pole placement (tune poles as needed)
    desired_poles = [-1,-1.2,-1.5,-0.8,-1.0,-1.3,-0.9,-0.95,-1.1];
    if length(desired_poles) ~= size(A_aug,1)
        error('Adjust desired_poles to length %d', size(A_aug,1));
    end
    K_aug = place(A_aug, B_aug, desired_poles);
    disp('Pole placement succeeded.');
end

% --- Partition gains ---
Kx = K_aug(:, 1:na);      % 3x6
Ki = K_aug(:, na+1:end);  % 3x3

disp('Kx ='); disp(Kx);
disp('Ki ='); disp(Ki);

% % --- Load and open model ---
% if ~bdIsLoaded(model)
%     load_system(model);
% end
% open_system(model);
% 
% % --- Block paths (as discovered) ---
% blk_Kp = [model '/Gain'];       % Gain — aocs/Gain
% blk_Ki = [model '/Gain2'];      % Gain2 — aocs/Gain2
% blk_Int = [model '/Integrator'];% Integrator — aocs/Integrator
% blk_Sum1 = [model '/Sum1'];     % Sum1
% blk_Sum2 = [model '/Sum2'];     % Sum2
% blk_Gain1 = [model '/Gain1'];   % Gain1
% 
% % --- Verify existence of critical blocks ---
% req = {blk_Kp, blk_Ki, blk_Int, blk_Sum1, blk_Sum2, blk_Gain1};
% for i=1:numel(req)
%     if isempty(find_system(model,'SearchDepth',1,'Path',req{i}))
%         warning('Block not found: %s. Update path if block is inside a subsystem.', req{i});
%     end
% end
% 
% % --- Write Kx and Ki into blocks ---
% set_param(blk_Kp, 'Gain', mat2str(Kx));
% set_param(blk_Ki, 'Gain', mat2str(Ki));
% 
% % --- Ensure Integrator is 3-channel and set initial condition ---
% if ~isempty(find_system(model,'SearchDepth',1,'Path',blk_Int))
%     set_param(blk_Int, 'InitialCondition', 'zeros(3,1)');
%     disp('Integrator initial condition set to zeros(3,1). Integrator is expected to be 3-channel.');
% else
%     % If integrator somehow missing, add one
%     add_block('simulink/Continuous/Integrator', blk_Int, 'Position', [200 300 250 360], 'InitialCondition', 'zeros(3,1)');
%     disp('Integrator block added at top level.');
% end
% 
% % --- Add Saturation block (ControllerSaturation) if not present ---
% satBlk = [model '/ControllerSaturation'];
% if isempty(find_system(model,'SearchDepth',1,'Name','ControllerSaturation'))
%     add_block('simulink/Commonly Used Blocks/Saturation', satBlk, 'Position', [450 300 500 340]);
%     set_param(satBlk, 'UpperLimit', '1', 'LowerLimit', '-1');
%     disp('ControllerSaturation block added at top level. Manual wiring may be needed.');
% else
%     disp('ControllerSaturation already present.');
% end
% 
% % --- Save model ---
% save_system(model);
% disp('Gains written to the model and integrator configured. IMPORTANT: verify wiring manually:');
% disp(' - Ensure attitude error signals (3 channels) feed aocs/Integrator output.');
% disp(' - Ensure aocs/Gain2 (Ki, 3x3) multiplies integrator output before summing into the actuator path.');
% disp(' - Ensure aocs/Gain (Kx, 3x6) multiplies the 6-state vector as intended.');
% disp('Run simulation and inspect scopes. If automatic wiring is incorrect, I can provide exact add_line/delete_line commands once you confirm where the attitude error signals originate.');
