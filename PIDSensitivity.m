%% PIDSensitivity.m
% This script performs a full parameter sweep over the four PID gains:
% kU, kP, kI, and kD. For each combination, it runs a grip simulation (via runGripSim)
% that returns a performance metric (final z displacement) and a stability flag.
%
% Because the parameter space is 4D, the script uses a parallel coordinates
% plot and pairwise scatter plots to help visualize the results.

clear; clc; close all;

kU_range = linspace(0, 0, 1);   
kP_range = linspace(2, 2, 1);      
kI_range = linspace(0, 3, 5);     
kD_range = linspace(0, 1, 5);      

% Define the stability criterion threshold 
force_threshold = 2; 

%% Preallocate a results array.
% Columns: [kU, kP, kI, kD, metric]
N = length(kU_range) * length(kP_range) * length(kI_range) * length(kD_range);
results = zeros(N, 5);
stability_flags = false(N, 1);  % Flag for stability

counter = 1;

%% Loop over every combination of kU, kP, kI, and kD
for iU = 1:length(kU_range)
    for iP = 1:length(kP_range)
        for iI = 1:length(kI_range)
            for iD = 1:length(kD_range)
                current_kU = kU_range(iU);
                current_kP = kP_range(iP);
                current_kI = kI_range(iI);
                current_kD = kD_range(iD);
                
                % Run Simulting Function
                [metric, isStable] = runGripSim(current_kU, current_kP, current_kI, current_kD, force_threshold);
                
                results(counter, :) = [current_kU, current_kP, current_kI, current_kD, metric];
                stability_flags(counter) = isStable;
                
                fprintf('kU=%.2f, kP=%.2f, kI=%.2f, kD=%.2f => Final z=%.4f m, Stable: %d\n',...
                    current_kU, current_kP, current_kI, current_kD, metric, isStable);
                counter = counter + 1;
            end
        end
    end
end

%% Visualization
% 3D Scatter Plot: kP vs. kU vs. Final z displacement
figure;
scatter3(results(:,2), results(:,3), results(:,5), 50, results(:,5), 'filled');
xlabel('kP');
ylabel('kI');
zlabel('Final z displacement (m)');
title('3D Scatter: kP vs. kI vs. Final z displacement');
colorbar;
grid on;

% 3D Scatter Plot: kP vs. kD vs. Final z displacement
figure;
scatter3(results(:,2), results(:,4), results(:,5), 50, results(:,5), 'filled');
xlabel('kP');
ylabel('kD');
zlabel('Final z displacement (m)');
title('3D Scatter: kP vs. kD vs. Final z displacement');
colorbar;
grid on;

% 3D Scatter Plot: kP vs. kU vs. Final z displacement
figure;
scatter3(results(:,2), results(:,1), results(:,5), 50, results(:,5), 'filled');
xlabel('kP');
ylabel('kU');
zlabel('Final z displacement (m)');
title('3D Scatter: kP vs. kU vs. Final z displacement');
colorbar;
grid on;

%% Print all combinations
% Filter out only the stable combinations
stable_results = results(stability_flags, :);

fprintf('\nStable PID Combinations (kU, kP, kI, kD, Metric):\n');
disp(stable_results);

% Compute the range for each gain that are stable
if ~isempty(stable_results)
    min_kU = min(stable_results(:,1));
    max_kU = max(stable_results(:,1));
    
    min_kP = min(stable_results(:,2));
    max_kP = max(stable_results(:,2));
    
    min_kI = min(stable_results(:,3));
    max_kI = max(stable_results(:,3));
    
    min_kD = min(stable_results(:,4));
    max_kD = max(stable_results(:,4));
    
    fprintf('Stable kU range: [%.3f, %.3f]\n', min_kU, max_kU);
    fprintf('Stable kP range: [%.3f, %.3f]\n', min_kP, max_kP);
    fprintf('Stable kI range: [%.3f, %.3f]\n', min_kI, max_kI);
    fprintf('Stable kD range: [%.3f, %.3f]\n', min_kD, max_kD);
else
    fprintf('No stable combinations found.\n');
end

%% Function for the simulation
function [metric, isStable] = runGripSim(kU, kP, kI, kD, force_threshold)
    % runGripSim runs the grip simulation with the specified PID gains.
    % It returns a performance metric (final z displacement), the
    % stability, and the full force output time series (force_output).
    %
    % Inputs:
    %   kU, kP, kI, kD : PID gains.
    %   force_threshold: stability constant
    %
    % Outputs:
    %   metric       : Final z displacement (m)
    %   isStable     : Boolean flag on whether the force graph is stable

    % Simulation parameters
    T_total = 5;  % seconds
    dt = 0.01;    % seconds
    time_vector = (0:dt:T_total)';
    
    % A single weight and friction
    weight = 50;    % grams
    mu = 0.2; % coefficient
    m = weight / 1000;  % mass in kg
    g_val = 9.81;
    
    % Compute threshold force and desired force
    F_thresh = (m * g_val) / mu;
    F_desired = F_thresh * 1.10;
    
    % Create reference force profile
    t_break = [0, 1, 2, 3, 4, 5];
    F_break = [0, F_desired, F_desired, F_desired, F_desired, 0];
    F_ref = interp1(t_break, F_break, time_vector, 'pchip');
    
    % Initialize models
    initial_state = zeros(4, 1);
    muscle_model = Hammerstein(initial_state);
    observer_model = StateObserver(initial_state);
    
    numSteps = length(time_vector);
    output_forces = zeros(1, numSteps);
    z_disp = zeros(1, numSteps);  % z displacement in m
    z_vel = zeros(1, numSteps);   % z velocity in m/s
    PW_f_history = zeros(1, numSteps);
    PW_e_history = zeros(1, numSteps);
    
    prev_y = 0;
    integral_error = 0;
    prev_error = 0;
    t_delay = 0; 
    
    % Main simulation loop
    for i = 1:numSteps
        current_ref = F_ref(i);
        error = current_ref - prev_y;
        integral_error = integral_error + error * dt;
        derivative_error = (error - prev_error) / dt;

        K = [484.18, 15.4, -518.1, -15.75;
             -559,   -162.79, 605.4, 17.1];
        u_intermediate = [1, -1] * (-K * observer_model.xk_bar_hat);
        
        % Combine PID terms with state feedback, scaled by kU
        uc = kU * u_intermediate + kP * error + kI * integral_error + kD * derivative_error + prev_y;
        
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
        
        % Update the muscle model and observer
        [muscle_model, y] = muscle_model.update(PW_f, PW_e);
        if i < 1 + t_delay
            observer_model = observer_model.update(0, 0, y);
        else
            observer_model = observer_model.update(PW_f_history(i-t_delay), PW_e_history(i-t_delay), y);
        end
        
        output_forces(i) = y;
        
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
    
    % Use final z displacement as metric
    metric = z_disp(end);
    
    % Determine stability:
    windowSize = round(0.5 / dt);

    % Compute the trend using a moving average
    trend = movmean(output_forces, windowSize);
    
    % Compute the detrended output
    detrended_output = output_forces - trend;
    
    % Compute the oscillation amplitude from the detrended output
    oscillation_amplitude = max(detrended_output) - min(detrended_output);
    if oscillation_amplitude > force_threshold
        isStable = false;
    else
        isStable = true;
    end
end