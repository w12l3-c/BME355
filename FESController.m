function [PW_f, PW_e, integral_term, prev_error] = FESController(x_hat, r, y, K, kU, kP, kI, kD, dt, i_term, prev_error)
% FESController - Compute the FES control inputs for an antagonist muscle pair
%               using a combination of state feedback and PID error correction.
%
% Syntax:
%   [PW_f, PW_e, integral_term, prev_error] = FESController(x_hat, r, y, K, kU, kP, kI, kD, dt, i_term, prev_error)
%
% Inputs:
%   x_hat       - Estimated state vector (4x1) from the Hammerstein muscle model.
%   r           - Reference force (desired grip force, scalar).
%   y           - Measured force (actual grip force, scalar).
%   K           - Feedback gain matrix (2x4) where K = [K_f; K_e] (designed via LQR).
%   kU          - Scaling gain for the overall control signal.
%   kP          - Proportional gain applied to the force error (r - y).
%   kI          - Integral gain applied to the accumulated force error.
%   kD          - Derivative gain applied to the change in error.
%   dt          - Time step in seconds (used for integration and differentiation).
%   i_term      - Current accumulated integral term (scalar), passed from the previous cycle.
%   prev_error  - Previous error value (scalar), used to compute the derivative of the error.
%
% Outputs:
%   PW_f        - Computed pulse width for the flexor muscle (scalar).
%   PW_e        - Computed pulse width for the extensor muscle (scalar).
%   integral_term - Updated integral term (accumulated error) for the next cycle.
%   prev_error  - Updated error value (current error) for the next cycle.

    % Compute the error between reference and measured force
    error = r - y;

    % Integral term updated using Euler integration
    integral_term = i_term + error * dt;
    I = kI * integral_term;
    
    % Derivative term calculated based on difference in error
    D = kD * (error - prev_error) / dt;
    
    % Compute the state feedback term
    % Then premultiply by [1 -1] to get a scalar intermediate control signal.
    disp(x_hat);
    u_intermediate = [1, -1] * (-K * x_hat); % (1,1)
    
    % Add the proportional term based on force error 
    uc = kU * u_intermediate + (kP * error) + I + D + y;
    prev_error = error; % Update previous error term for next iteration
    
    fprintf('Error (r-y): %.2f, u_int: %.2f, P: %.2f, I: %.2f, D: %.2f, uc: %.2f\n', ...
         error, kU*u_intermediate, kP*error, I, D, uc);
    
    % Apply switching rule to separate stimulation for flexor and extensor:
    if uc >= 0
        % flexor receives positive control signal
        uf = uc;  
        ue = 0;
    else
        % extensor receives the magnitude of the negative control signal
        uf = 0; 
        ue = abs(uc); 
    end
    ubar = [uf ue];
    ubar = transpose(ubar);
    fprintf("Generate Ubar: %d\n", ubar);

    [PW_f, PW_e] = InverseIRC(uf, ue);

end
