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

reference_force = 1;

output_forces = zeros(1, numSteps);
%% Simulation Loop
% We will simulate the discrete system over "numSteps" steps.
for i = 1:numSteps

    if i ~= 1
        reference_force = 0;
    end
    
    [PW_f, PW_e] = FESController(observer_model.xk_bar_hat, reference_force, prev_y, K, kP);
    

    [muscle_model, y] = muscle_model.update(PW_f, PW_e);

    observer_model = observer_model.update(PW_f, PW_e, y);
    
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
