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

% Define weight and friction ranges
weights = 50:75:50;    % weights in grams
frictions = 0.2:0.1:0.2;  % friction coefficients

g_val = 9.81;  % gravitational acceleration
alpha_z = 1e-4;  % scaling factor to convert force (N) to upward velocity (m/s)
F_thresh = 0;    % force threshold

% Controller gains 
kP = 2;
kI = 0:0.5:0.1;
kD = 0:0.5:0.1;
K = [484.18, 15.4, -518.1, -15.75;
     -559,   -162.79, 605.4, 17.1];

totalSims = length(weights) * length(frictions);
simCount = 1;

n_rows = length(weights);   
n_cols = length(frictions);

numW = length(weights);
numMu = length(frictions);
final_z = zeros(numW, numMu);

figure;

for w_idx = 1:numW
    for mu_idx = 1:numMu
        % Convert weight from grams to kg
        m = weights(w_idx) / 1000;  % mass in kg
        mu = frictions(mu_idx);
        
        % Compute desired grip force (in Newtons) as Fn = (m*g)/mu
        F_desired = (m * g_val) / mu;
        
        t_break = [0, 1, 2, 3, 4, 5];   
        noise_level = 0.01 * F_desired;  
        
        F_break = [0, ...
                   F_desired + noise_level * randn, ... 
                   F_desired + noise_level * randn, ... 
                   F_desired + noise_level * randn, ... 
                   F_desired + noise_level * randn, ... 
                   0];

        F_ref = interp1(t_break, F_break, time_vector, 'pchip');
        % F_ref = awgn(F_ref, 5);    % gaussian
        
        %% Initialize Models
        initial_state = zeros(4, 1);
        muscle_model = Hammerstein(initial_state);
        observer_model = StateObserver(initial_state);

        numSteps = length(time_vector);
        output_forces = zeros(1, numSteps);
        z_disp = zeros(1, numSteps);  % z displacement in m
        z_vel = zeros(1, numSteps);
        
        prev_y = 0;
        
        %% Simulation Loop for this combination
        for i = 1:numSteps
            current_ref = F_ref(i);
            [PW_f, PW_e] = FESController(observer_model.xk_bar_hat, current_ref, prev_y, K, kP);
            [muscle_model, y] = muscle_model.update(PW_f, PW_e);
            observer_model = observer_model.update(PW_f, PW_e, y);
            output_forces(i) = y;
            
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
            prev_y = y;
        end
        
        final_z(w_idx, mu_idx) = z_disp(end);
        
        fprintf('Weight = %dg, mu = %.2f, Final z displacement = %.4f m\n', ...
                weights(w_idx), mu, z_disp(end));

        index = (w_idx - 1) * n_cols + mu_idx;

        subplot(n_rows, n_cols, index);
        plot(time_vector, F_ref, 'k--', 'LineWidth', 1.5); hold on;
        plot(time_vector, output_forces, 'b-', 'LineWidth', 1.5);
        title(sprintf('%dg, \\mu=%.2f\nz_{final}=%.3f m', weights(w_idx), mu, z_disp(end)), 'FontSize', 8);
        xlabel('t (s)', 'FontSize', 8); ylabel('Force (N)', 'FontSize', 8);
        grid on;
        
    end
end

%% Plot

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
