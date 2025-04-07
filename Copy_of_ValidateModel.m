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
C = [543.656,  0,  -679.57, 0];

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
kP = 2;
%K = [484.18, 15.4, -518.1, -15.75;
%     -559, -162.79, 605.4, 17.1];

output_forces = zeros(1, numSteps);
history_pwf = zeros(1, numSteps+2);
history_pwe = zeros(1, numSteps+2);

% history_pwf(1) = 0;
% history_pwf(1) = 0;
% history_pwe(1) = 0;
% history_pwe(2) = 0;
desired_forces = zeros(1, numSteps);

%% Simulation Loop
% We will simulate the discrete system over "numSteps" steps.
for i = 1:numSteps

    % Gets the reference force at the ith 
    if i < 150
        reference_force = 10;
    else
        reference_force = 5;
    end

    desired_forces(i) = reference_force;
    % Calculates the PW for FES stimulation
    [PW_f, PW_e] = FESController(observer_model.xk_bar_hat, reference_force, prev_y, K, kP);

    history_pwf(i) = PW_f;
    history_pwe(i) = PW_e;
    
    % Updates Hammerstein muscle model
    [muscle_model, y] = muscle_model.update(PW_f, PW_e);
    % disp(muscle_model.xk_bar)

    % Updates observer model
    observer_model = observer_model.update(PW_f, PW_e, prev_y);

    % Save force output for plotting
    output_forces(i) = y;

    % Set variable for next iteration for error calculation
    prev_y = y;
end

%% Plot the Results

% Create a time vector for discrete steps
time = 0:numSteps-1;

plot(time, output_forces, 'b-', time, desired_forces, 'r--');

% plot(time, history_pwf(1:numSteps), 'b-', time, history_pwe(1:numSteps), 'r--');

% Add labels and title
xlabel('x values');
ylabel('y values');
title('Plot of x vs y');
