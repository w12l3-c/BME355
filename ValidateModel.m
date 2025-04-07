%% ValidateImpulse.m
% Clear workspace and command window
clear; clc;

%% Define the system matrices from the hammerstein2 class (constants)
initial_state = zeros(4, 1);
muscle_model = Hammerstein(initial_state);
observer_model = StateObserver(initial_state);

%%
Phi = [0.82  0.008  0     0;
       0     0.82   0     0;
       0     0      0.78  0.008;
       0     0      0     0.78];

% Input matrix
Gamma = [0        0;
         0.009    0;
         0        0;
         0        0.009];

% Output matrix
C = [54.3656,  0,  -67.957, 0];

Q = 1000*(C*C');
R = eye(size(Gamma, 2));
[K, s, p] = dlqr(Phi, Gamma, Q, R);

%% Simulation Parameters

% Get reference forces from Figure 6
data = readmatrix("resampled_data.csv");
time_data = data(:, 1); % Gets array of timepoints
reference_forces = data(:, 2); % Gets array of reference force data

numSteps = length(time_data);       % Number of discrete time steps
prev_y = 0;
kP = 1.25;
%   K = [484.18, 15.4, -518.1, -15.75;
%        -559, -162.79, 605.4, 17.1];

reference_force = 0;
output_forces = zeros(1, numSteps);
history_pwf = zeros(1, numSteps+2);
history_pwe = zeros(1, numSteps+2);

history_pwf(1) = 0;
history_pwf(1) = 0;
history_pwe(1) = 0;
history_pwe(2) = 0;

%% Simulation Loop
% We will simulate the discrete system over "numSteps" steps.
for i = 1:numSteps

    % Gets the reference force at the ith sample
    reference_force = reference_forces(i);

    % Calculates the PW for FES stimulation
    [PW_f, PW_e] = FESController(observer_model.xk_bar_hat, reference_force, prev_y, K, kP);

    history_pwf(i+2) = PW_f;
    history_pwe(i+2) = PW_e;
    
    % Updates Hammerstein muscle model
    [muscle_model, y] = muscle_model.update(history_pwf(i), history_pwe(i));

    % Updates observer model
    observer_model = observer_model.update(history_pwf(i), history_pwe(i), y);

    % Save force output for plotting
    output_forces(i) = y;

    % Set variable for next iteration for error calculation
    prev_y = y;
end

%% Plot the Results

% Create a time vector for discrete steps
time = 0:numSteps-1;

% plot(time, output_forces, 'b-', time, reference_forces, 'r--');

plot(time, history_pwf(1:numSteps), 'b-', time, history_pwe(1:numSteps), 'r--');

% Add labels and title
xlabel('x values');
ylabel('y values');
title('Plot of x vs y');
