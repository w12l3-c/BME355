%% GripSim.m
% This script simulates a grip experiment over 5 seconds.
% The simulation runs for different object weights (50 to 500 g in steps of 75 g)
% and different friction coefficients (e.g., 0.2 to 0.9).
%
% For each simulation, the reference grip force is set as:
%   Fn = m*g/mu   (with m in kg, g=9.81 m/s^2).
%
% The simulation updates the muscle model (Hammerstein), observer, and FESController,
% and then integrates a simple kinematic model to compute z displacement.
%
% Finally, a surface plot is created showing final z displacement vs. weight and friction.

clear; clc; close all;

%% Simulation Setup Parameters
T_total = 5;       % Total simulation time in seconds
dt = 0.01;         % Sampling period (s)
time_vector = (0:dt:T_total)';
t_delay = 0;

% Define weight and friction ranges
weights = 50:75:50;    % weights in grams
frictions = 0.2:0.1:0.2;  % friction coefficients

g_val = 9.81;  % gravitational acceleration
F_thresh = 0;    % force threshold

% Controller gains 
w_f = 19.8451;
w_e = 24.8461;
w_f = 20;
w_e = 25;

A = [-w_f 1 0 0
    0 -w_f 0 0
    0 0 -w_e 1
    0 0 0 -w_e];
B = [0 0
    1 0
    0 0
    0 1];

A_f = A(1:2, 1:2);
A_e = A(3:4, 3:4);
Phi_f = exp(1)^(A_f*0.01);
Phi_e = exp(1)^(A_e*0.01);
Phi = exp(1)^(A*0.01);
Gamma = inv(A)*(Phi - eye(4))*B;

C = [5436.56,  0,  -6795.7,  0];
% State matrix
Phi = [0.82  0.008  0     0;
       0     0.82   0     0;
       0     0      0.78  0.008;
       0     0      0     0.78];

% Input matrix
Gamma = [0       0;
         0.009   0;
         0       0;
         0       0.009];
Q = 1e3 * (C' * C);
R = eye(size(Gamma,2));  % R must be the same size as the number of inputs
[K, S, e] = dlqr(Phi, Gamma, Q, R);

kU = 1;
kP = 2;
kI = 110;
kD = 0;
K = [484.18, 15.4, -518.1, -15.75;
     -559,   -162.79, 605.4, 17.1];
%%
kU_arr = 0:1:0.05;
kP_arr = 0.25:2:0.25;
kI_arr = 10:150:10;
kD_arr = 0:1:0.05;

totalSims = length(weights) * length(frictions);
simCount = 1;

n_rows = length(weights);   
n_cols = length(frictions);

numW = length(weights);
numMu = length(frictions);
final_z = zeros(numW, numMu);

forceFigure = figure('Name', 'Reference vs. Output Force');
pwFigure = figure('Name', 'FES Pulse Widths');

for w_idx = 1:numW
    for mu_idx = 1:numMu
        % Convert weight from grams to kg
        simCount = simCount + 1;
        m = weights(w_idx) / 1000;  % mass in kg
        mu = frictions(mu_idx);
        
        % Compute desired grip force (in Newtons) as Fn = (m*g)/mu
        F_thresh = (m * g_val) / mu;
        F_desired = F_thresh * 1.10;
        
        t_break = [0, 1, 2, 3, 4, 5];   
        noise_level = 0.05 * F_desired;  
        
        F_break = [0, ...
                   F_desired , ... 
                   F_desired , ... 
                   F_desired , ... 
                   F_desired , ... 
                   0];

        F_thresh_break = [0, ...
                   F_thresh , ... 
                   F_thresh , ... 
                   F_thresh , ... 
                   F_thresh , ... 
                   0];

        F_ref = interp1(t_break, F_break, time_vector, 'pchip');
        F_ori = interp1(t_break, F_thresh_break, time_vector, "pchip");
        % F_ref = awgn(F_ref, 5);    % gaussian
        
        %% Initialize Models
        initial_state = zeros(4, 1);
        muscle_model = Hammerstein(initial_state);
        observer_model = StateObserver(initial_state);

        % observer_model.stability();

        numSteps = length(time_vector);
        output_forces = zeros(1, numSteps);
        z_disp = zeros(1, numSteps);  % z displacement in m
        z_vel = zeros(1, numSteps);   % z velocity in m/s

        PW_f_history = zeros(1, numSteps);  % to store computed PW_f at each step
        PW_e_history = zeros(1, numSteps);  % to store computed PW_e at each step
        stateHistory = zeros(4, numSteps);
        
        prev_y = 0;
        i_term = 0;  % Initialize the integral of error
        prev_error = 0;      % Initialize the previous error

        
        %% Simulation Loop for this combination
        for i = 1:numSteps
            current_ref = F_ref(i);
            [PW_f, PW_e, i_term, prev_error] = FESController(observer_model.xk_bar_hat, current_ref, prev_y, K, kU, kP, kI, kD, dt, i_term, prev_error);
            % PW_f = PW_f * 1e1;
            % PW_e = PW_e * 1e1;
            PW_f_history(i) = PW_f;
            PW_e_history(i) = PW_e;
            [muscle_model, y] = muscle_model.update(PW_f, PW_e);
            
            if i < 1 + t_delay
                % [muscle_model, y] = muscle_model.update(0, 0);
                observer_model = observer_model.update(0, 0, y);
            else
                % [muscle_model, y] = muscle_model.update(PW_f_history(i-t_delay), PW_e_history(i-t_delay));
                observer_model = observer_model.update(PW_f_history(i-t_delay), PW_e_history(i-t_delay), y);
            end

            output_forces(i) = y;
            stateHistory(:, i) = muscle_model.xk_bar;
            fprintf("Output: %d\n", y);
            
            % Compute z velocity
            if time_vector(i) < 1 | time_vector(i) > 4
                a_z = 0;
                z_vel(i) = 0;
            else
                if y < F_thresh
                    a_z = (y * mu - m * g_val) / m;
                else
                    a_z = 0;
                    z_vel(i) = 0;
                end
            end
            
            % Update z displacement:
            if i < numSteps
                z_vel(i+1) = z_vel(i) + a_z * dt;
                z_disp(i+1) = z_disp(i) + z_vel(i+1) * dt;
            end
            
            % Update previous measured force
            prev_y = output_forces(i);
        end
        
        final_z(w_idx, mu_idx) = z_disp(end);
        
        fprintf('Weight = %dg, mu = %.2f, Final z displacement = %.4f m\n', ...
                weights(w_idx), mu, z_disp(end));

        index = (w_idx - 1) * n_cols + mu_idx;

        % Plotting: Reference vs. Output Force in forceFigure
        figure(forceFigure);
        subplot(n_rows, n_cols, index);
        plot(time_vector, F_ref, 'k--', 'LineWidth', 1.5); hold on;
        plot(time_vector, F_ori, 'r-', 'LineWidth', 1.5); hold on;
        plot(time_vector, output_forces, 'b-', 'LineWidth', 1.5);
        title(sprintf('%dg, \\mu=%.2f\nz=%.3f m', weights(w_idx), mu, z_disp(end)), 'FontSize', 8);
        xlabel('Time (s)', 'FontSize', 7);
        ylabel('Force (N)', 'FontSize', 7);
        legend('Ref', 'Thres', 'Out', 'Location', 'Best', 'FontSize', 6);
        grid on;
        
        % Plotting: FES Pulse Widths in pwFigure
        figure(pwFigure);
        subplot(n_rows, n_cols, index);
        plot(time_vector, PW_f_history, 'm-', 'LineWidth', 1.5); hold on;
        plot(time_vector, PW_e_history, 'c-', 'LineWidth', 1.5);
        title(sprintf('%dg, \\mu=%.2f', weights(w_idx), mu), 'FontSize', 8);
        xlabel('Time (s)', 'FontSize', 7);
        ylabel('Pulse Width', 'FontSize', 7);
        legend('PW_f', 'PW_e', 'Location', 'Best', 'FontSize', 6);
        grid on;

        
    end
end

%% Plot
% State plot
figure('Name','State Variables over Time');
subplot(2,2,1);
plot(time_vector, stateHistory(1,:), 'b-', 'LineWidth',1.5);
title('State x_1');
xlabel('Time (s)');
ylabel('x_1');

subplot(2,2,2);
plot(time_vector, stateHistory(2,:), 'r-', 'LineWidth',1.5);
title('State x_2');
xlabel('Time (s)');
ylabel('x_2');

subplot(2,2,3);
plot(time_vector, stateHistory(3,:), 'g-', 'LineWidth',1.5);
title('State x_3');
xlabel('Time (s)');
ylabel('x_3');

subplot(2,2,4);
plot(time_vector, stateHistory(4,:), 'm-', 'LineWidth',1.5);
title('State x_4');
xlabel('Time (s)');
ylabel('x_4');

sgtitle('Time Evolution of Hammerstein Model States');

% Create grid for weights and frictions
[MU, W] = meshgrid(frictions, weights);

figure;
surf(MU, W, final_z);
xlabel('Friction Coefficient (\mu)');
ylabel('Weight (grams)');
zlabel('Final Z Displacement (m)');
title('Final Z Displacement for Various Weights and Friction Coefficients');
colorbar;
shading interp;

% contour plot:
figure;
contourf(MU, W, final_z, 20, 'LineColor', 'none');
xlabel('Friction Coefficient (\mu)');
ylabel('Weight (grams)');
title('Contour Plot of Final Z Displacement (m)');
colorbar;

