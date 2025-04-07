%% ValidateImpulse.m
% Clear workspace and command window
clear; clc;

%% Define the system matrices from the hammerstein2 class (constants)
initial_state = zeros(4, 1);
muscle_model = Hammerstein(initial_state);
observer_model = StateObserver(initial_state);

%% Simulation Parameters

% Get reference forces from Figure 6
data = readmatrix("resampled_data.csv");
time_data = data(:, 1); % Gets array of timepoints
reference_forces = data(:, 2); % Gets array of reference force data

numSteps = length(time_data);       % Number of discrete time steps
prev_y = 0;
kP = 0.125;
%K = [484.18, 15.4, -518.1, -15.75;
%     -559, -162.79, 605.4, 17.1];
reference_force = 0;
output_forces = zeros(1, numSteps);
history_pwf = zeros(1, numSteps+2);
history_pwe = zeros(1, numSteps+2);

%% Simulation Loop
% We will simulate the discrete system over "numSteps" steps.
for i = 1:numSteps
    
    % if i < 150
    %     desired_force = 10;
    % else
    %     desired_force = 0;
    % end
    % 
    % error = desired_force - prev_y;
    % 
    % 
    % if error > 0
    %     % Activate the flexors
    %     [PW_f, PW_e] = InverseIRC(desired_force, 0);
    %     disp("flexor");
    % else
    %     % Activate the extensors
    %     [PW_f, PW_e] = InverseIRC(0, desired_force);
    %     disp("extensor");
    % end


    if i < 150
        [PW_f, PW_e] = InverseIRC(10, 0);
    elseif i < 160
        [PW_f, PW_e] = InverseIRC(0, 1);
    else
        [PW_f, PW_e] = InverseIRC(prev_y, 0);
    end
    disp(PW_f);
    % Gets the reference force at the ith sample

    
    % Updates Hammerstein muscle model
    [muscle_model, y] = muscle_model.update(PW_f, PW_e);

    % Updates observer model
    observer_model = observer_model.update(PW_f, PW_e, y);

    % Save force output for plotting
    output_forces(i) = y;

    % Set variable for next iteration for error calculation
    prev_y = y;
end

%% Plot the Results

% Create a time vector for discrete steps
time = 0:numSteps-1;

plot(time, output_forces, 'b-');

% plot(time, history_pwf(1:numSteps), 'b-', time, history_pwe(1:numSteps), 'r--');

% Add labels and title
xlabel('x values');
ylabel('y values');
title('Plot of x vs y');
