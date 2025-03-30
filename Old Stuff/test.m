% === SETUP ===
clear; clc;

% Simulation settings
T = 50;              % Number of time steps
r = 30;              % Desired grip force (target)
kP = 1.2;            % Proportional gain

% Initial state and pulsewidth
x = zeros(4,1);      % Initial state
%PW = 0.0002;         % Initial pulse width
PW_f = 0;
PW_e = 0;

% === SYSTEM MATRICES FROM hammerstein CLASS ===
Phi = [0.82 0.008 0 0;
       0 0.82 0 0;
       0 0 0.78 0.008;
       0 0 0 0.78];

Gamma = [0 0;
         0.009 0;
         0 0;
         0 0.009];

% === LQR CONTROLLER DESIGN ===
C = [5436.56, 0, -6795.7, 0];
Q = 1000 * (C' * C);
R = eye(2);
K = dlqr(Phi, Gamma, Q, R);  % Compute optimal gain matrix

% === DATA LOGGING ===
force_log = zeros(1,T);
uf_log = zeros(1,T);
ue_log = zeros(1,T);

% === SIMULATION LOOP ===
for k = 1:T
    % Create Hammerstein model object
    model = hammerstein2(PW_f, PW_e, x);   % Use your class name exactly
    
    % Get output grip force
    y = model.yk();  % Output force
    
    % FES control input
    u = FESController(x, r, y, K, kP);  % Returns [uf; ue]
    uf = u(1);
    ue = u(2);
    
    % Store results
    uf_log(k) = uf;
    ue_log(k) = ue;
    force_log(k) = y;

    % Update system state
    x = model.xk1_bar();

    PW_f = uf;
    PW_e = ue;
end

% === PLOTTING ===
figure;
subplot(3,1,1);
plot(force_log, 'LineWidth', 1.5);
yline(r, 'r--', 'Target');
title('Grip Force'); ylabel('Force'); grid on;

subplot(3,1,2);
plot(uf_log, 'LineWidth', 1.5);
title('Flexor Stimulation (uf)'); ylabel('uf'); grid on;

subplot(3,1,3);
plot(ue_log, 'LineWidth', 1.5);
title('Extensor Stimulation (ue)'); ylabel('ue'); xlabel('Time'); grid on;
