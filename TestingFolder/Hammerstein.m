classdef Hammerstein
    %% Setting up constants for the Hammerstein model
    properties (Constant)
        % Recruitment curve constants
        c1_flex = 4.14;
        c2_flex = 2655.88;
        c1_ext  = 4.75;
        c2_ext  = 913.2;

        % State matrix
        Phi = [0.82  0.008  0     0;
               0     0.82  0     0;
               0     0     0.78  0.008;
               0     0     0     0.78];

        % Input matrix
        Gamma = [0    0;
                 0.009    0;
                 0    0;
                 0    0.009];

        % Output matrix
        C = [5436.56  0  -6795.7  0];

        % System delay
        tau = 0.02; %20 ms
    end

    properties
        % State vector
        xk_bar (4,1) {mustBeNumeric} = zeros(4,1);
    end

    methods
        %% Constructor
        function obj = Hammerstein(initial_state)
            obj.xk_bar = initial_state;
        end

        %% Generate one time step
        function [obj, output] = update(obj, PW_f, PW_e)

            % Get the input PWM forces
            u_bar = [
                obj.c1_flex * abs(tanh(obj.c2_flex * PW_f / 2));
                obj.c1_ext  * abs(tanh(obj.c2_ext  * PW_e / 2))
            ];

            % 4x4 * 4x1 + 4x2 * 2*1
            next_xk_bar = obj.Phi * obj.xk_bar + obj.Gamma * u_bar; % Calculate the next x_bar (4, 1)
    
            % 1x4 * 4*1
            output = obj.C * obj.xk_bar; % Get the force output (1, 1)

            obj.xk_bar = next_xk_bar; % Update u for next time step
        end
    end
end
