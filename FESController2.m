function u = FESController2(x, r, y, K, kP)
% Compute error between desired and actual grip force
e = r - y;

% Compute base control effort for each muscle
u_feedback = -K * x;         % 2x1: [uf_raw; ue_raw]
u_prop = [kP * e; -kP * e];  % Encourage flexor if e > 0, extensor if e < 0

% Combine feedback and proportional term
u_raw = u_feedback + u_prop;

% Clamp negative values to 0 (muscles can't be stimulated negatively)
uf = max(u_raw(1), 0);
ue = max(u_raw(2), 0);

u = [uf; ue];
end
