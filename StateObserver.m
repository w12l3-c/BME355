classdef StateObserver
    %% Setting up constants for the Hammerstein model
    properties (Constant)
        % Recruitment curve constants
        c1_flex = 4.14;
        c2_flex = 2655.88;
        c1_ext  = 4.75;
        c2_ext  = 913.2;

        % State matrix
        Phi = [0.82  0.008  0     0;
               0     0.82   0     0;
               0     0      0.78  0.008;
               0     0      0     0.78];

        % Input matrix
        Gamma = [0        0;
                 0.009    0;
                 0        0;
                 0        0.009];

        % Output matrix
        C = [5436.56,  0,  -6795.7, 0];
        % C = [54365.6,  0,  -67957,  0];

        % System delay
        tau = 0.02; %20 ms
    end

    properties
        % State vector
        xk_bar_hat (4,1) {mustBeNumeric} = zeros(4,1);
        PW_history
        k
    end

    methods
        %% Constructor
        function obj = StateObserver(initial_state)
            obj.xk_bar_hat = initial_state;
            obj.PW_history = [];
            empty_col_vector = zeros(2, 1);
            obj.PW_history = [obj.PW_history, empty_col_vector];
            obj.PW_history = [obj.PW_history, empty_col_vector];
            obj.k = 1;
        end

        %% Generate one time step
        function [obj] = update(obj, PW_f, PW_e, yk)

            % Get the input PWM forces
            u_bar = [
                obj.c1_flex * abs(tanh(obj.c2_flex * PW_f / 2));
                obj.c1_ext  * abs(tanh(obj.c2_ext  * PW_e / 2))
            ]; % (2, 1)
            % u_bar = [
            %     PW_f;
            %     PW_e
            % ]; % (2, 1)

            obj.PW_history = [obj.PW_history, u_bar];
            % 4x4 * 4x1       
            % poles = [19.8451, 19.85, 24.8461, 24.85];
            poles = [0.82, 0.81, 0.78, 0.77];
            % poles = [0.72, 0.71, 0.58, 0.57];

            L = place(obj.Phi', obj.C', poles)'; % (4x2)


            % 4x4 * 4x1 + 4x2 * 2x1 + (4x2 * 1x1) = 4x1
            % next_xk_bar_hat = obj.Phi * obj.xk_bar_hat + obj.Gamma * obj.PW_history(:, obj.k) + L * (yk - obj.C* obj.xk_bar_hat); % Calculate the next x_bar
            next_xk_bar_hat = obj.Phi * obj.xk_bar_hat + obj.Gamma * u_bar + L * (yk - obj.C * obj.xk_bar_hat);

            obj.xk_bar_hat = next_xk_bar_hat; % Update u for next time step
            obj.k = obj.k+1;
        end

        function obj = stability(obj)
            poles = [0.82, 0.81, 0.78, 0.77];
            L = place(obj.Phi', obj.C', poles)'; % (4x2)
            A_minus_LC = obj.Phi - L * obj.C;
            eigVals = eig(A_minus_LC);

            theta = linspace(0, 2*pi, 100);
            unit_circle = exp(1i*theta);
            figure;
            plot(real(unit_circle), imag(unit_circle), 'k--','LineWidth',1.5); hold on;
            plot(real(eigVals), imag(eigVals), 'ro','MarkerSize',8, 'LineWidth',2);
            xlabel('Real Part'); ylabel('Imaginary Part');
            title('Observer Error Dynamics: Eigenvalues vs. Unit Circle');
            axis equal;
            grid on;
        end
    end
end
