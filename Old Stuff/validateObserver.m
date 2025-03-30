%% validateObserver.m
% This script validates the full-order state observer for the hammerstein2 model.
% It computes the observer gain L via pole placement and then simulates the system
% and observer over a number of steps.

% Clear workspace and command window
clear; clc;

%% Define the system matrices from the hammerstein2 class (constants)
Phi = hammerstein2.Phi;
Gamma = hammerstein2.Gamma;
C = hammerstein2.C;

% Choose desired observer poles (they should be placed inside the unit circle for a discrete system).
% For example, for a 4-state system:
desiredObserverPoles = [0.5, 0.45, 0.4, 0.35];

% Compute observer gain L using the duality approach:
L = place(Phi', C', desiredObserverPoles)';
disp('Observer Gain L:');
disp(L);

%% Simulation Parameters
numSteps = 50;       % Number of discrete time steps
n = size(Phi,1);     % Number of states (should be 4)

% Preallocate arrays for true states and observer estimates
x_true = zeros(n, numSteps);
xhat = zeros(n, numSteps);

% Initial Conditions:
% Let's assume the true state starts at some nonzero value and the observer starts at zero.
x_true(:,1) = [1000; 0; 1000; 0];  % Example initial state (units as defined by your model)
xhat(:,1) = zeros(n,1);            % Observer's initial state estimate

% Define constant stimulation pulse widths for this simulation
PW_f = 200e-6;  % Flexor pulse width (example value)
PW_e = 150e-6;  % Extensor pulse width (example value)

%% Simulation Loop
% We will simulate the discrete system over "numSteps" steps.
for k = 1:numSteps-1
    % Create a hammerstein2 object for the current time step using the true state.
    % (We assume that the current pulse widths are maintained constant.)
    hModel = hammerstein2(PW_f, PW_e, x_true(:,k));
    
    % Compute the input vector using the recruitment curves:
    % This yields a 2x1 vector corresponding to the flexor and extensor contributions.
    u = hModel.u_bar();
    
    % True system update:
    % x_true(:,k+1) = Phi * x_true(:,k) + Gamma * u;
    x_true(:,k+1) = Phi * x_true(:,k) + Gamma * u;
    
    % Measured output from the true system:
    yk = C * x_true(:,k);
    
    % Observer update:
    % The observer uses the same input and measures the error (yk - C*xhat)
    xhat(:,k+1) = Phi * xhat(:,k) + Gamma * u + L * (yk - C * xhat(:,k));
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
