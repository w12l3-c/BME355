function u = FESController(x, r, y, K, kP)
% FESController - Compute FES control input for an antagonist muscle pair.
%
% Syntax:  u = FESController(x, r, y, K, kP)
%
% Inputs:
%    x  - State vector (4x1) from the Hammerstein muscle model (or its estimate)
%    r  - Reference force (desired grip force, scalar)
%    y  - Measured force (actual grip force, scalar)
%    K  - Feedback gain matrix (2x4) where K = [Kf; Ke] (designed via LQR)
%    kP - Proportional gain for output error (scalar)
%
% Outputs:
%    u  - Input vector for FES [uf; ue] where:
%           uf is the stimulation command for the flexor muscle,
%           ue is the stimulation command for the extensor muscle.
%
% The controller computes an intermediate control signal:
%    uc = [1 -1] * (-K * x) + kP*(r - y)
%
% Then it applies a switching rule:
%    if uc >= 0, then uf = uc and ue = 0,
%    else,         uf = 0 and ue = abs(uc).
%
% This premultiplication by [1 -1] avoids coactivation of both muscles.
%
% Example:
%    % x: current state (4x1), r: desired force, y: measured force
%    % K: feedback gain matrix, kP: proportional gain
%    u = FESController(x, r, y, K, kP);
%

    % Compute the error between reference and measured force
    error = r - y;
    
    % Compute the state feedback term
    % Then premultiply by [1 -1] to get a scalar intermediate control signal.
    u_intermediate = [1, -1] * (-K * x);
    
    % Add the proportional term based on force error (The sum circle)
    uc = u_intermediate + kP * error;
    
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
    
    % Output the FES input vector [uf; ue]
    u = [uf; ue];
end
