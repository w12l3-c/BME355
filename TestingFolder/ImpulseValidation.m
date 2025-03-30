%% ValidateModel.m
% Clear workspace and command window
clear; clc;

%% Define the system matrices from the hammerstein2 class (constants)
initial_state = zeros(4, 1);
muscle_model = Hammerstein(initial_state);
observer_model = StateObserver(initial_state);


%% Simulation Parameters
numSteps = 50;       % Number of discrete time steps
prev_y = 0;
kP = 0.25;
K = [484.18, 15.4, -518.1, -15.75;
     -559, -162.79, 605.4, 17.1];


PW_f = 0;
PW_e = 1;

output_forces = zeros(1, numSteps);
%% Simulation Loop
% We will simulate the discrete system over "numSteps" steps.
for i = 1:numSteps

    if i ~= 1
        PW_f = 0;
        PW_e = 0
    end
    
    [muscle_model, y] = muscle_model.update(PW_f, PW_e);
    
    output_forces(i) = y;
    prev_y = y;
end

%% Plot the Results

% Create a time vector for discrete steps
time = 0:numSteps-1;

plot(time, output_forces);

% Add labels and title
xlabel('x values');
ylabel('y values');
title('Plot of x vs y');
