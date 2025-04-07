%% GripSim.m
% This script simulates a grip experiment over 3 seconds.
% With t=0 to 1 increasing the output force from 0 and t=4 to 5 decreasing
% the output force back to 0
%
% The simulation runs for different object weights (50 to 350 g in steps of 75 g)
% and different friction coefficients (e.g., 0.3 to 0.9 in steps of 0.2).
%
% For each simulation, the reference grip force is set as:
%   Fn = m*g/mu * 1.1   (with m in kg, g=9.81 m/s^2).
% Where the 1.1 is a factor to make sure the force output is larger than
% the force required
%
% The simulation updates the muscle model (Hammerstein) observer, and PID integrated FESController,
% and then integrates a kinematic model to compute z displacement.
%
% Finally, a surface plot is created showing final z displacement vs. weight and friction.

clear; clc; close all;

%% Simulation Setup Parameters
T_total = 5;       % Total simulation time in seconds
dt = 0.01;         % Sampling period (s)
time_vector = (0:dt:T_total)';
t_delay = 0;       % Time delay (in dt seconds)

% Define weight and friction ranges
weights = 50:75:275;    % weights in grams
frictions = 0.3:0.2:0.9;  % friction coefficients

g_val = 9.81;  % gravitational acceleration
F_thresh = 0;    % force threshold

% Controller gains 
C = [543.656,  0,  -679.57, 0];
kU = 0;
kP = 2;
kI = 100;
kD = 0.1;
K = [484.18, 15.4, -518.1, -15.75;
     -559,   -162.79, 605.4, 17.1];

%% Plotting Parameters
totalSims = length(weights) * length(frictions);
simCount = 1;

n_rows = length(weights);   
n_cols = length(frictions);

numW = length(weights);
numMu = length(frictions);
final_z = zeros(numW, numMu);

forceFigure = figure('Name', 'Reference vs. Output Force');
pwFigure = figure('Name', 'FES Pulse Widths');

%% Simulation
for w_idx = 1:numW
    for mu_idx = 1:numMu
        % Convert weight from grams to kg
        simCount = simCount + 1;
        m = weights(w_idx) / 1000;  % mass in kg
        mu = frictions(mu_idx);     % friction coefficient
        
        % Compute desired grip force (in Newtons) as Fn = (m*g)/mu
        F_thresh = (m * g_val) / mu;
        F_desired = F_thresh * 1.10;
        
        % Create reference force profile
        t_break = [0, 1, 2, 3, 4, 5];   
        noise_level = 0.05 * F_desired;  
        
        % F = 0 when t = 0 & 5, F = ref_force when 1 < t < 4
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

        % Interpolating both reference and threshold curve for smooth
        % transition
        F_ref = interp1(t_break, F_break, time_vector, 'pchip');
        F_ori = interp1(t_break, F_thresh_break, time_vector, "pchip");
        
        %% Initialize Models
        initial_state = zeros(4, 1);
        muscle_model = Hammerstein(initial_state);
        observer_model = StateObserver(initial_state);

        % observer_model.stability();   % Stability circle
        
        numSteps = length(time_vector);     % Number of time steps
        output_forces = zeros(1, numSteps); % Array to save the output forces
        z_disp = zeros(1, numSteps);  % z displacement in m
        z_vel = zeros(1, numSteps);   % z velocity in m/s

        PW_f_history = zeros(1, numSteps);  % to store computed PW_f at each step
        PW_e_history = zeros(1, numSteps);  % to store computed PW_e at each step
        stateHistory = zeros(4, numSteps);  % to store the state variables at each step
        
        prev_y = 0;  % Previous output
        i_term = 0;  % Integral of error
        prev_error = 0;      % Previous error

        
        %% Simulation Loop for this combination of mu and weight
        for i = 1:numSteps
            current_ref = F_ref(i);

            % FES Controller
            [PW_f, PW_e, i_term, prev_error] = FESController(observer_model.xk_bar_hat, current_ref, prev_y, K, kU, kP, kI, kD, dt, i_term, prev_error);
            
            PW_f_history(i) = PW_f; % Saving Pulse Width for flexion 
            PW_e_history(i) = PW_e; % Saving Pulse Width for extension  
            
            % Time delay for state observer and Hammerstein update
            if i < 1 + t_delay
                [muscle_model, y] = muscle_model.update(0, 0);
                observer_model = observer_model.update(0, 0, y);
            else
                [muscle_model, y] = muscle_model.update(PW_f_history(i-t_delay), PW_e_history(i-t_delay));
                observer_model = observer_model.update(PW_f_history(i-t_delay), PW_e_history(i-t_delay), y);
            end
            stateHistory(:, i) = muscle_model.xk_bar;   % Save current state variable

            output_forces(i) = y; % Output forces
            fprintf("Output: %d\n", y);
            
            % Compute z velocity
            % Doesn't compute velocity as the cup is not grasp
            if time_vector(i) < 1 | time_vector(i) > 4
                a_z = 0;
                z_vel(i) = 0;
            else
                % Output force larger than the threshold (not reference, reference is 1.1 * threshold)
                if y < F_ori(i)
                    % acceleration base on F = (N*mu - mg)
                    a_z = (y * mu - m * g_val) / m; 
                else
                    a_z = 0;
                    z_vel(i) = 0;   
                end
            end
            
            % Update z displacement base on kinematics
            if i < numSteps
                z_vel(i+1) = z_vel(i) + a_z * dt;   % Update velocity
                z_disp(i+1) = z_disp(i) + z_vel(i+1) * dt;  % Update displacement
            end
            
            % Update previous measured force
            prev_y = output_forces(i);
        end
        
        % Final z displacement
        final_z(w_idx, mu_idx) = z_disp(end);
        fprintf('Weight = %dg, mu = %.2f, Final z displacement = %.4f m\n', ...
                weights(w_idx), mu, z_disp(end));

        % Index for plotting
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

%% More Plots
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

% contour plot for weight * mu over z_displacement
figure;
contourf(MU, W, final_z, 20, 'LineColor', 'none');
xlabel('Friction Coefficient (\mu)');
ylabel('Weight (grams)');
title('Contour Plot of Final Z Displacement (m)');
colorbar;

