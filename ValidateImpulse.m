%% ValidateModel.m
% Clear workspace and command window
clear; clc;

%% Define the system matrices from the hammerstein2 class (constants)
initial_state = zeros(4, 1);
muscle_model = Hammerstein(initial_state);
observer_model = StateObserver(initial_state);


%% Simulation Parameters
numSteps = 50;       % 50*0.01s = 500ms
kP = 0.25;
K = [484.18, 15.4, -518.1, -15.75;
     -559, -162.79, 605.4, 17.1];

% If you want to impulse just the flexors, set PW_f to 1 and PW_e to 0
% If you want to impulse just the extensors, set PW_f to 0 and PW_e to 1
PW_f = 1;
PW_e = 0;

% Array to store output forces for plotting
output_forces = zeros(1, numSteps);
%% Simulation Loop
for i = 1:numSteps

    % Ensure that after the initial impulse, the forces become 0
    if i ~= 1
        PW_f = 0;
        PW_e = 0;
    end
    
    % Stimulate muscle model and get output force
    [muscle_model, y] = muscle_model.update(PW_f, PW_e);
    
    % Store output force
    output_forces(i) = y;

end

%% Plot the Results

% Create a time vector for discrete steps
time = 0:numSteps-1;

plot(time, output_forces);

% Add labels and title
xlabel('x values');
ylabel('y values');
title('Plot of x vs y');
