%% ValidateModel.m
% Clear workspace and command window
clear; clc;

%% Define the system matrices from the hammerstein2 class (constants)
inital_state = zeros(4, 1);
muscle_model = Hammerstein(initial_state);
observer_model = StateObserver(initial_state);


%% Simulation Parameters
numSteps = 100;       % Number of discrete time steps
prev_y = 0;
kP = 0.25;
K = [484.18, 15.4, -518.1, -15.75;
     -559, -162.79, 605.4, 17.1];

reference_force = 0
%% Simulation Loop
% We will simulate the discrete system over "numSteps" steps.
for k = 1:numSteps
    reference_force = 0 %TODO: GET VALUE AT TIME K
    
    [PW_f, PW_e] = FESController(observer_model.xk_bar_hat, reference_force, prev_y, K, kP);

    [muscle_model, y] = muscle_model.update(PW_f, PW_e);

    observer_model = observer_model.update(PW_f, PW_e, y);

end

%% Plot the Results

% Create a time vector for discrete steps
time = 0:numSteps-1;

% Plot each state: compare true state and observer estimate.
figure;
for i = 1:n
    subplot(n,1,i);
    plot(time, x_true(i,:), 'b-', 'LineWidth', 2); hold on;
    plot(time, xhat(i,:), 'r--', 'LineWidth', 2);
    xlabel('Discrete Time Step');
    ylabel(['State x', num2str(i)]);
    legend('True', 'Observer Estimate');
    title(['State x', num2str(i), ' vs. Observer Estimate']);
    grid on;
end

% Plot the output error over time:
y_true = C * x_true;
y_hat = C * xhat;
figure;
plot(time, y_true - y_hat, 'k-', 'LineWidth', 2);
xlabel('Discrete Time Step');
ylabel('Output Error (y_{true} - y_{hat})');
title('Observer Output Error');
grid on;
