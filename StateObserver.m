function xhat_next = StateObserver(Phi, Gamma, C, L, xhat, uBuffer, k, tauSteps, y)
% StateObserver updates the state estimate for a delayed-input system.
%
% Syntax:
%   xhat_next = StateObserver(Phi, Gamma, C, L, xhat, uBuffer, k, tauSteps, y)
%
% Inputs:
%   Phi       - State transition matrix (n x n)
%   Gamma     - Input matrix (n x m)
%   C         - Output matrix (p x n)
%   L         - Observer gain matrix (n x p)
%   xhat      - Current state estimate (n x 1)
%   uBuffer   - Array storing past inputs, size (m x totalSteps)
%   k         - Current time index (integer)
%   tauSteps  - Number of discrete steps of delay (integer)
%   y         - Current measured output (p x 1)
%
% Output:
%   xhat_next - Updated state estimate (n x 1)

    % 1. Compute the delayed index for the input
    delayedIndex = k - tauSteps;
    
    % 2. If the delayed index is less than or equal to zero, assume zero input
    if delayedIndex <= 0
        uDelayed = zeros(size(uBuffer, 1), 1);
    else
        uDelayed = uBuffer(:, delayedIndex);
    end
    
    % 3. Compute the observer correction term based on output error
    Correction = L * (y - C * xhat);
    
    % 4. Update the observer state
    xhat_next = Phi * xhat + Gamma * uDelayed + Correction;
end
