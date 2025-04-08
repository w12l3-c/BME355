%% PIDSweep.m
% This script sweeps through a range of kI and kD values (with kU and kP fixed)
% and for each combination it runs the simulation (via runGripSim) and then plots
% the output force time series. The subplots are arranged in a 5x5 grid.

clear; clc; close all;

% Define parameter ranges:
kU_range = linspace(0, 0, 1);   % Only one value (0)
kP_range = linspace(1, 5, 1);   
kI_range = [0, 3, 10, 100, 500];
kD_range = [0, 0.01, 0.05, 0.1, 0.5];

% Define simulation settings (if not defined elsewhere)
T_total = 5;      % Total simulation time in seconds
dt = 0.01;        % Time step (s)
time_vector = (0:dt:T_total)';

% kU and kP are constants
kU = kU_range(1);
kP = kP_range(2);

numI = length(kI_range);
numD = length(kD_range);

% Create a figure for the 5x5 subplots
figure('Name','Force Output Time Series (kI vs. kD)');
for iI = 1:numI
    for iD = 1:numD
        current_kI = kI_range(iI);
        current_kD = kD_range(iD);
        
        % Run the simulation for this combination.
        % We assume runGripSim returns [metric, isStable, force_output].
        [~, force_output, F_ref, F_ori] = runGripSim(kU, kP, current_kI, current_kD, force_threshold);
        
        % Determine the subplot index (row-wise ordering)
        subplot(numI, numD, (iI-1)*numD + iD);
        plot(time_vector, F_ref, 'k--', 'LineWidth', 1.5); hold on;
        plot(time_vector, F_ori, 'r-', 'LineWidth', 1.5); hold on;
        plot(time_vector, force_output, 'b-', 'LineWidth',1.5);
        title(sprintf('kI=%.2f, kD=%.2f', current_kI, current_kD), 'FontSize',8);
        xlabel('Time (s)', 'FontSize',7);
        ylabel('Force (N)', 'FontSize',7);
        legend('Ref', 'Thres', 'Out', 'Location', 'Best', 'FontSize', 6);
        grid on;
    end
end

%% Function for the simulation
function [metric, output_forces, F_ref, F_ori] = runGripSim(kU, kP, kI, kD)
    % runGripSim runs the grip simulation with the specified PID gains.
    % It returns a performance metric (final z displacement)
    % and the full force output time series (force_output).
    %
    % Inputs:
    %   kU, kP, kI, kD : PID gains.
    %
    % Outputs:
    %   metric       : Final z displacement (m)
    %   force_output : Vector of output force values over time
    %   F_ref        : Reference Force
    %   F_ori        : Threshold Force
    
    % Simulation parameters
    T_total = 5;  % seconds
    dt = 0.01;    % time step (s)
    time_vector = (0:dt:T_total)';
    
    % Use a single weight and friction
    weight = 50;      % grams
    friction = 0.2;   % coefficient of friction
    m = weight / 1000;  % mass in kg
    g_val = 9.81;
    
    % Compute threshold force and desired force
    F_thresh = (m * g_val) / friction;
    F_desired = F_thresh * 1.10;
    
    % Create a reference force profile using pchip interpolation
    t_break = [0, 1, 2, 3, 4, 5];
    F_break = [0, F_desired, F_desired, F_desired, F_desired, 0];
    F_ref = interp1(t_break, F_break, time_vector, 'pchip');
    F_break = [0, F_thresh, F_thresh, F_thresh, F_thresh, 0];
    F_ori = interp1(t_break, F_break, time_vector, 'pchip');
    
    % Initialize models
    initial_state = zeros(4, 1);
    muscle_model = Hammerstein(initial_state);
    observer_model = StateObserver(initial_state);
    
    numSteps = length(time_vector);
    output_forces = zeros(1, numSteps);
    z_disp = zeros(1, numSteps);  % z displacement (m)
    z_vel = zeros(1, numSteps);  % z velocity (m/s)
    PW_f_history = zeros(1, numSteps);
    PW_e_history = zeros(1, numSteps);
    
    % Initialize PID state variables
    prev_y = 0;
    integral_error = 0;
    prev_error = 0;
    t_delay = 0;  
    
    % Fixed gain matrix K
    K = [484.18, 15.4, -518.1, -15.75;
         -559,   -162.79, 605.4, 17.1];
    
    % Simulation Loop
    for i = 1:numSteps
        current_ref = F_ref(i);
        error = current_ref - prev_y;
        integral_error = integral_error + error * dt;
        derivative_error = (error - prev_error) / dt;
        
        % Compute state feedback term using observer's state estimate:
        u_intermediate = [1, -1] * (-K * observer_model.xk_bar_hat);
        
        % Combine PID terms scaled by kU
        uc = kU * u_intermediate + kP * error + kI * integral_error + kD * derivative_error + prev_y;
        
        % Switching rule: if uc >= 0, stimulate flexor; if uc < 0, stimulate extensor.
        if uc >= 0
            uf = uc;
            ue = 0;
        else
            uf = 0;
            ue = abs(uc);
        end
        
        % Compute pulse widths using your InverseIRC function
        [PW_f, PW_e] = InverseIRC(uf, ue);
        PW_f_history(i) = PW_f;
        PW_e_history(i) = PW_e;
        
        % Update muscle model and observer
        [muscle_model, y] = muscle_model.update(PW_f, PW_e);
        if i < 1 + t_delay
            observer_model = observer_model.update(0, 0, y);
        else
            observer_model = observer_model.update(PW_f_history(i-t_delay), PW_e_history(i-t_delay), y);
        end
        
        output_forces(i) = y;
        
        % Update z displacement (active between t=1 and t=4 seconds)
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
        
        prev_y = output_forces(i);
        prev_error = error;
    end
    
    % Use final z displacement as the performance metric.
    metric = z_disp(end);
end
